using OpenTK.Graphics.OpenGL.Compatibility;
using OpenTK.Mathematics;
using KimeraCS.Core;

namespace KimeraCS.Rendering
{
    using static Utils;

    /// <summary>
    /// Modern skeleton renderer using VAO/VBO instead of immediate mode GL.
    /// Computes bone transforms on CPU and renders using LineMesh/PointMesh.
    /// </summary>
    public static class SkeletonRenderer
    {
        // Cached meshes for bone visualization
        private static LineMesh? _boneLinesMesh;
        private static PointMesh? _boneJointsMesh;

        /// <summary>
        /// Compute world-space bone positions for a unified skeleton.
        /// Works with both field and battle skeletons through the unified format.
        /// </summary>
        public static void ComputeBonePositions(
            UnifiedSkeleton skeleton,
            UnifiedFrame frame,
            out List<Vector3> jointPositions,
            out List<Vector3> boneStarts,
            out List<Vector3> boneEnds)
        {
            jointPositions = new List<Vector3>();
            boneStarts = new List<Vector3>();
            boneEnds = new List<Vector3>();

            if (skeleton?.Bones == null || skeleton.Bones.Count == 0 || frame == null)
                return;

            // Stack for hierarchical transforms
            var matrixStack = new Stack<Matrix4>();
            var parentStack = new Stack<int>();

            // Bone direction multiplier (field = -1, battle = +1)
            float boneDir = skeleton.BoneDirection == BoneDirection.NegativeZ ? -1f : 1f;

            // Start with root transform (Y already negated in unified frame for field)
            Matrix4 rootTranslation = Matrix4.CreateTranslation(frame.RootTranslation);

            Matrix4 rootRotation = BuildRotationMatrixWithQuaternions(
                frame.RootRotation.Alpha,
                frame.RootRotation.Beta,
                frame.RootRotation.Gamma);

            Matrix4 currentMatrix = rootRotation * rootTranslation;
            matrixStack.Push(currentMatrix);
            parentStack.Push(-1); // Root has no parent

            for (int i = 0; i < skeleton.Bones.Count; i++)
            {
                var bone = skeleton.Bones[i];

                // Pop matrices until we find the parent
                while (parentStack.Count > 0 && parentStack.Peek() != bone.ParentIndex)
                {
                    if (matrixStack.Count > 1) // Keep at least root matrix
                    {
                        matrixStack.Pop();
                        parentStack.Pop();
                    }
                    else
                    {
                        break;
                    }
                }

                currentMatrix = matrixStack.Peek();

                // Push current matrix for children
                matrixStack.Push(currentMatrix);
                parentStack.Push(i);

                // Apply bone rotation
                if (i < frame.BoneRotations.Count)
                {
                    var rot = frame.BoneRotations[i];
                    Matrix4 boneRotation = BuildRotationMatrixWithQuaternions(
                        rot.Alpha, rot.Beta, rot.Gamma);
                    currentMatrix = boneRotation * currentMatrix;
                }

                // Get bone start position (origin in current transform)
                Vector3 boneStart = Vector3.TransformPosition(Vector3.Zero, currentMatrix);

                // Get bone end position (along Z axis with bone direction)
                Vector3 boneEnd = Vector3.TransformPosition(
                    new Vector3(0, 0, boneDir * bone.Length),
                    currentMatrix);

                jointPositions.Add(boneStart);
                jointPositions.Add(boneEnd);
                boneStarts.Add(boneStart);
                boneEnds.Add(boneEnd);

                // Update stack matrix with translation along bone
                Matrix4 boneTranslation = Matrix4.CreateTranslation(0, 0, boneDir * bone.Length);
                matrixStack.Pop();
                matrixStack.Push(boneTranslation * currentMatrix);
            }
        }

        /// <summary>
        /// Render unified skeleton bones using modern OpenGL.
        /// </summary>
        public static void RenderSkeletonBones(
            UnifiedSkeleton skeleton,
            UnifiedFrame frame,
            float lineR = 1f, float lineG = 1f, float lineB = 1f,
            float jointR = 1f, float jointG = 0f, float jointB = 0f)
        {
            if (!GLRenderer.IsReady || skeleton == null || frame == null) return;

            ComputeBonePositions(skeleton, frame,
                out var jointPositions, out var boneStarts, out var boneEnds);

            if (boneStarts.Count == 0) return;

            // Sync legacy GL matrices to GLRenderer for modern rendering
            double[] projMatrix = new double[16];
            GL.GetDouble(GetPName.ProjectionMatrix, projMatrix);
            var legacyProjection = ToMatrix4(projMatrix);

            double[] mvMatrix = new double[16];
            GL.GetDouble(GetPName.ModelviewMatrix, mvMatrix);
            var legacyModelView = ToMatrix4(mvMatrix);

            // Save original matrices
            var savedProjection = GLRenderer.ProjectionMatrix;
            var savedView = GLRenderer.ViewMatrix;
            var savedModel = GLRenderer.ModelMatrix;

            // Use legacy matrices directly
            GLRenderer.ProjectionMatrix = legacyProjection;
            GLRenderer.ViewMatrix = Matrix4.Identity;
            GLRenderer.ModelMatrix = legacyModelView;

            // Build line vertices for bones
            var lineVertices = new LineVertex[boneStarts.Count * 2];
            for (int i = 0; i < boneStarts.Count; i++)
            {
                lineVertices[i * 2] = new LineVertex
                {
                    Position = boneStarts[i],
                    Color = new Vector4(lineR, lineG, lineB, 1f)
                };
                lineVertices[i * 2 + 1] = new LineVertex
                {
                    Position = boneEnds[i],
                    Color = new Vector4(lineR, lineG, lineB, 1f)
                };
            }

            // Build point vertices for joints
            var pointVertices = new LineVertex[jointPositions.Count];
            for (int i = 0; i < jointPositions.Count; i++)
            {
                pointVertices[i] = new LineVertex
                {
                    Position = jointPositions[i],
                    Color = new Vector4(jointR, jointG, jointB, 1f)
                };
            }

            // Create or update meshes
            if (_boneLinesMesh == null)
                _boneLinesMesh = new LineMesh();
            _boneLinesMesh.Upload(lineVertices);

            if (_boneJointsMesh == null)
                _boneJointsMesh = new PointMesh();
            _boneJointsMesh.Upload(pointVertices);

            // Render using GLRenderer
            GLRenderer.DrawLinesModern(_boneLinesMesh);
            GLRenderer.DrawPointsModern(_boneJointsMesh, 5.0f);

            // Restore original matrices
            GLRenderer.ProjectionMatrix = savedProjection;
            GLRenderer.ViewMatrix = savedView;
            GLRenderer.ModelMatrix = savedModel;
        }

        /// <summary>
        /// Cleanup cached meshes.
        /// </summary>
        public static void Cleanup()
        {
            _boneLinesMesh?.Dispose();
            _boneLinesMesh = null;

            _boneJointsMesh?.Dispose();
            _boneJointsMesh = null;
        }
    }
}
