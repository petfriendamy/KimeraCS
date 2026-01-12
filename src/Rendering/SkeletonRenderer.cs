using System;
using System.Collections.Generic;
using OpenTK.Graphics.OpenGL.Compatibility;
using OpenTK.Mathematics;
using KimeraCS.Core;

namespace KimeraCS.Rendering
{
    using static FF7FieldSkeleton;
    using static FF7FieldAnimation;
    using static FF7BattleSkeleton;
    using static FF7BattleAnimation;
    using static Utils;

    /// <summary>
    /// Modern skeleton renderer using VAO/VBO instead of immediate mode GL.
    /// Computes bone transforms on CPU and renders using LineMesh/PointMesh.
    /// </summary>
    public static class SkeletonRenderer
    {
        // Cached meshes for bone visualization
        private static LineMesh _boneLinesFieldMesh;
        private static PointMesh _boneJointsFieldMesh;
        private static LineMesh _boneLinesBattleMesh;
        private static PointMesh _boneJointsBattleMesh;

        /// <summary>
        /// Compute world-space bone positions for a field skeleton.
        /// Returns arrays of joint start/end positions for rendering.
        /// </summary>
        public static void ComputeFieldBonePositions(
            FieldSkeleton fSkeleton,
            FieldFrame fFrame,
            out List<Vector3> jointPositions,
            out List<Vector3> boneStarts,
            out List<Vector3> boneEnds)
        {
            jointPositions = new List<Vector3>();
            boneStarts = new List<Vector3>();
            boneEnds = new List<Vector3>();

            if (fSkeleton.bones == null || fSkeleton.bones.Count == 0)
                return;

            // Stack for hierarchical transforms
            var matrixStack = new Stack<Matrix4>();
            var jointStack = new string[fSkeleton.bones.Count + 1];
            int jsp = 0;

            // Start with root transform
            Matrix4 rootTranslation = Matrix4.CreateTranslation(
                (float)fFrame.rootTranslationX,
                (float)-fFrame.rootTranslationY,
                (float)fFrame.rootTranslationZ);

            Matrix4 rootRotation = CreateRotationMatrixFromQuaternions(
                fFrame.rootRotationAlpha,
                fFrame.rootRotationBeta,
                fFrame.rootRotationGamma);

            Matrix4 currentMatrix = rootRotation * rootTranslation;
            matrixStack.Push(currentMatrix);

            jointStack[jsp] = fSkeleton.bones[0].joint_f;

            for (int iBoneIdx = 0; iBoneIdx < fSkeleton.bones.Count; iBoneIdx++)
            {
                // Pop matrices until we find matching joint
                while (fSkeleton.bones[iBoneIdx].joint_f != jointStack[jsp] && jsp > 0)
                {
                    if (matrixStack.Count > 0)
                        currentMatrix = matrixStack.Pop();
                    jsp--;
                }

                // Push current matrix for this bone
                matrixStack.Push(currentMatrix);

                // Apply bone rotation
                Matrix4 boneRotation = CreateRotationMatrixFromQuaternions(
                    fFrame.rotations[iBoneIdx].alpha,
                    fFrame.rotations[iBoneIdx].beta,
                    fFrame.rotations[iBoneIdx].gamma);

                currentMatrix = boneRotation * currentMatrix;

                // Get bone start position (origin in current transform)
                Vector3 boneStart = Vector3.TransformPosition(Vector3.Zero, currentMatrix);

                // Get bone end position
                Vector3 boneEnd = Vector3.TransformPosition(
                    new Vector3(0, 0, (float)-fSkeleton.bones[iBoneIdx].len),
                    currentMatrix);

                jointPositions.Add(boneStart);
                jointPositions.Add(boneEnd);
                boneStarts.Add(boneStart);
                boneEnds.Add(boneEnd);

                // Translate along bone for next iteration
                Matrix4 boneTranslation = Matrix4.CreateTranslation(0, 0, (float)-fSkeleton.bones[iBoneIdx].len);
                currentMatrix = boneTranslation * currentMatrix;

                jsp++;
                jointStack[jsp] = fSkeleton.bones[iBoneIdx].joint_i;
            }
        }

        /// <summary>
        /// Compute world-space bone positions for a battle skeleton.
        /// </summary>
        public static void ComputeBattleBonePositions(
            BattleSkeleton bSkeleton,
            BattleFrame bFrame,
            out List<Vector3> jointPositions,
            out List<Vector3> boneStarts,
            out List<Vector3> boneEnds)
        {
            jointPositions = new List<Vector3>();
            boneStarts = new List<Vector3>();
            boneEnds = new List<Vector3>();

            if (bSkeleton.bones == null || bSkeleton.bones.Count == 0)
                return;

            // Stack for hierarchical transforms
            var matrixStack = new Stack<Matrix4>();
            int jsp = 0;

            // Start with root transform
            Matrix4 rootTranslation = Matrix4.CreateTranslation(
                (float)bFrame.startX,
                (float)bFrame.startY,
                (float)bFrame.startZ);

            Matrix4 rootRotation = CreateRotationMatrixFromQuaternions(
                bFrame.bones[0].alpha,
                bFrame.bones[0].beta,
                bFrame.bones[0].gamma);

            Matrix4 currentMatrix = rootRotation * rootTranslation;
            matrixStack.Push(currentMatrix);

            for (int bi = 0; bi < bSkeleton.nBones; bi++)
            {
                // Handle parent relationships (simplified - battle skeletons may have different hierarchy)
                while (bSkeleton.bones[bi].parentBone < jsp && jsp > 0)
                {
                    if (matrixStack.Count > 0)
                        currentMatrix = matrixStack.Pop();
                    jsp--;
                }

                matrixStack.Push(currentMatrix);

                // Apply bone rotation (skip first bone as it's handled by root)
                if (bi > 0 && bi + 1 < bFrame.bones.Count)
                {
                    Matrix4 boneRotation = CreateRotationMatrixFromQuaternions(
                        bFrame.bones[bi + 1].alpha,
                        bFrame.bones[bi + 1].beta,
                        bFrame.bones[bi + 1].gamma);

                    currentMatrix = boneRotation * currentMatrix;
                }

                // Get bone positions
                Vector3 boneStart = Vector3.TransformPosition(Vector3.Zero, currentMatrix);
                Vector3 boneEnd = Vector3.TransformPosition(
                    new Vector3(0, 0, (float)bSkeleton.bones[bi].len),
                    currentMatrix);

                jointPositions.Add(boneStart);
                jointPositions.Add(boneEnd);
                boneStarts.Add(boneStart);
                boneEnds.Add(boneEnd);

                // Translate along bone
                Matrix4 boneTranslation = Matrix4.CreateTranslation(0, 0, (float)bSkeleton.bones[bi].len);
                currentMatrix = boneTranslation * currentMatrix;

                jsp++;
            }
        }

        /// <summary>
        /// Create a rotation matrix from Euler angles using quaternions (matches legacy BuildRotationMatrixWithQuaternions).
        /// </summary>
        private static Matrix4 CreateRotationMatrixFromQuaternions(double alpha, double beta, double gamma)
        {
            // Convert to radians and create rotation
            float alphaRad = (float)(alpha * Math.PI / 180.0);
            float betaRad = (float)(beta * Math.PI / 180.0);
            float gammaRad = (float)(gamma * Math.PI / 180.0);

            // Build rotation in YXZ order (matches legacy code)
            Matrix4 rotY = Matrix4.CreateRotationY(betaRad);
            Matrix4 rotX = Matrix4.CreateRotationX(alphaRad);
            Matrix4 rotZ = Matrix4.CreateRotationZ(gammaRad);

            return rotZ * rotX * rotY;
        }

        /// <summary>
        /// Render field skeleton bones using modern OpenGL.
        /// </summary>
        public static void RenderFieldSkeletonBones(
            FieldSkeleton fSkeleton,
            FieldFrame fFrame,
            float lineR = 1f, float lineG = 1f, float lineB = 1f,
            float jointR = 1f, float jointG = 0f, float jointB = 0f)
        {
            if (!GLRenderer.IsReady) return;

            ComputeFieldBonePositions(fSkeleton, fFrame,
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
            if (_boneLinesFieldMesh == null)
                _boneLinesFieldMesh = new LineMesh();
            _boneLinesFieldMesh.Upload(lineVertices);

            if (_boneJointsFieldMesh == null)
                _boneJointsFieldMesh = new PointMesh();
            _boneJointsFieldMesh.Upload(pointVertices);

            // Render using GLRenderer
            GLRenderer.DrawLinesModern(_boneLinesFieldMesh);
            GLRenderer.DrawPointsModern(_boneJointsFieldMesh, 5.0f);

            // Restore original matrices
            GLRenderer.ProjectionMatrix = savedProjection;
            GLRenderer.ViewMatrix = savedView;
            GLRenderer.ModelMatrix = savedModel;
        }

        /// <summary>
        /// Render battle skeleton bones using modern OpenGL.
        /// </summary>
        public static void RenderBattleSkeletonBones(
            BattleSkeleton bSkeleton,
            BattleFrame bFrame,
            float lineR = 1f, float lineG = 1f, float lineB = 1f,
            float jointR = 1f, float jointG = 0f, float jointB = 0f)
        {
            if (!GLRenderer.IsReady) return;

            ComputeBattleBonePositions(bSkeleton, bFrame,
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
            if (_boneLinesBattleMesh == null)
                _boneLinesBattleMesh = new LineMesh();
            _boneLinesBattleMesh.Upload(lineVertices);

            if (_boneJointsBattleMesh == null)
                _boneJointsBattleMesh = new PointMesh();
            _boneJointsBattleMesh.Upload(pointVertices);

            // Render using GLRenderer
            GLRenderer.DrawLinesModern(_boneLinesBattleMesh);
            GLRenderer.DrawPointsModern(_boneJointsBattleMesh, 5.0f);

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
            _boneLinesFieldMesh?.Dispose();
            _boneLinesFieldMesh = null;

            _boneJointsFieldMesh?.Dispose();
            _boneJointsFieldMesh = null;

            _boneLinesBattleMesh?.Dispose();
            _boneLinesBattleMesh = null;

            _boneJointsBattleMesh?.Dispose();
            _boneJointsBattleMesh = null;
        }
    }
}
