using System;
using System.Drawing;
using System.Windows.Forms;
using OpenTK.Graphics.OpenGL.Compatibility;
using OpenTK.Mathematics;

#nullable enable
namespace KimeraCS
{
    using Rendering;
    using static Rendering.VisualizationHelpers;

    using static FF7Skeleton;
    using static FF7FieldSkeleton;
    using static FF7FieldAnimation;
    using static FF7FieldRSDResource;

    using static FF7BattleSkeleton;
    using static FF7BattleAnimation;
    using static FF7BattleAnimationsPack;

    using static FF7PModel;

    using static Lighting;

    using static Utils;
    using Core;

    public static class ModelDrawing
    {
        public static uint[] tex_ids = new uint[1];
        public const int LETTER_SIZE = 5;
        public const float DEFAULT_NORMAL_SCALE = 1.0f;
        public const int DEFAULT_NORMAL_COLOR = 2;


        //  ---------------------------------------------------------------------------------------------------
        //  ================================== GENERIC FIELD/BATTLE DRAW  =====================================
        //  ---------------------------------------------------------------------------------------------------
        // Cache for normals mesh
        private static LineMesh? _normalsMesh;

        public static void ShowNormals(PGroup Group, PPolygon[] Polys, Vector3[] Verts,
                                       Vector3[] Normals, int[] NormalsIndex,
                                       NormalsDisplayMode showNormals = NormalsDisplayMode.None,
                                       float normalsScale = DEFAULT_NORMAL_SCALE,
                                       int normalsColor = DEFAULT_NORMAL_COLOR)
        {
            if (Group.HiddenQ) return;

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

            // Create or update the normals mesh
            _normalsMesh?.Dispose();
            _normalsMesh = CreateNormalsMesh(Group, Polys, Verts, Normals, NormalsIndex,
                                             showNormals, normalsScale, normalsColor);

            if (_normalsMesh != null)
            {
                GLRenderer.DrawLinesModern(_normalsMesh);
            }

            // Restore original matrices
            GLRenderer.ProjectionMatrix = savedProjection;
            GLRenderer.ViewMatrix = savedView;
            GLRenderer.ModelMatrix = savedModel;
        }

        public static void DrawPModel(ref PModel Model, ref uint[] tex_ids, bool HideHiddenGroupsQ,
                                      RenderingContext? renderContext = null)
        {
            // Get current legacy matrices - these contain the full camera + bone transforms
            double[] projMatrix = new double[16];
            GL.GetDouble(GetPName.ProjectionMatrix,projMatrix);
            var legacyProjection = ToMatrix4(projMatrix);

            double[] mvMatrix = new double[16];
            GL.GetDouble(GetPName.ModelviewMatrix,mvMatrix);
            var legacyModelView = ToMatrix4(mvMatrix);

            // Save original matrices
            var savedProjection = GLRenderer.ProjectionMatrix;
            var savedView = GLRenderer.ViewMatrix;
            var savedModel = GLRenderer.ModelMatrix;

            // Use the legacy matrices directly for full compatibility
            // This ensures bone transforms and camera setup match exactly
            GLRenderer.ProjectionMatrix = legacyProjection;
            GLRenderer.ViewMatrix = Matrix4.Identity;
            GLRenderer.ModelMatrix = legacyModelView;

            GLRenderer.DrawPModelModern(ref Model, tex_ids, HideHiddenGroupsQ);

            // Show normals if enabled (must be done after model drawing)
            if (renderContext != null && renderContext.Options.NormalsDisplayMode != NormalsDisplayMode.None)
            {
                for (int g = 0; g < Model.Header.numGroups; g++)
                {
                    ShowNormals(Model.Groups[g], Model.Polys, Model.Verts,
                                Model.Normals, Model.NormalIndex,
                                renderContext.Options.NormalsDisplayMode,
                                renderContext.Options.NormalsScale,
                                renderContext.Options.NormalsColor);
                }
            }

            // Restore original matrices
            GLRenderer.ProjectionMatrix = savedProjection;
            GLRenderer.ViewMatrix = savedView;
            GLRenderer.ModelMatrix = savedModel;
        }

        public static void DrawPModelWireframe(ref PModel Model, bool HideHiddenGroupsQ)
        {
            // Get current legacy matrices - these contain the full camera + bone transforms
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

            // Use the legacy matrices directly for full compatibility
            GLRenderer.ProjectionMatrix = legacyProjection;
            GLRenderer.ViewMatrix = Matrix4.Identity;
            GLRenderer.ModelMatrix = legacyModelView;

            // Draw wireframe with black color
            GLRenderer.DrawPModelWireframe(ref Model, new Vector3(0, 0, 0), HideHiddenGroupsQ);

            // Restore original matrices
            GLRenderer.ProjectionMatrix = savedProjection;
            GLRenderer.ViewMatrix = savedView;
            GLRenderer.ModelMatrix = savedModel;
        }

        public static void DrawPModelPolygonColors(ref PModel Model, bool HideHiddenGroupsQ)
        {
            // Get current legacy matrices - these contain the full camera + bone transforms
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

            // Use the legacy matrices directly for full compatibility
            GLRenderer.ProjectionMatrix = legacyProjection;
            GLRenderer.ViewMatrix = Matrix4.Identity;
            GLRenderer.ModelMatrix = legacyModelView;

            // Draw with polygon colors
            GLRenderer.DrawPModelPolygonColors(ref Model, HideHiddenGroupsQ);

            // Restore original matrices
            GLRenderer.ProjectionMatrix = savedProjection;
            GLRenderer.ViewMatrix = savedView;
            GLRenderer.ModelMatrix = savedModel;
        }

        public static void DrawPModelDLists(ref PModel Model, ref uint[] tex_ids,
                                            RenderingContext renderContext)
        {
            // Get current legacy matrices - these contain the full camera + bone transforms
            double[] projMatrix = new double[16];
            GL.GetDouble(GetPName.ProjectionMatrix,projMatrix);
            var legacyProjection = ToMatrix4(projMatrix);

            double[] mvMatrix = new double[16];
            GL.GetDouble(GetPName.ModelviewMatrix,mvMatrix);
            var legacyModelView = ToMatrix4(mvMatrix);

            // Save original matrices
            var savedProjection = GLRenderer.ProjectionMatrix;
            var savedView = GLRenderer.ViewMatrix;
            var savedModel = GLRenderer.ModelMatrix;

            // Use the legacy matrices directly for full compatibility
            // This ensures bone transforms and camera setup match exactly
            GLRenderer.ProjectionMatrix = legacyProjection;
            GLRenderer.ViewMatrix = Matrix4.Identity;
            GLRenderer.ModelMatrix = legacyModelView;

            GLRenderer.DrawPModelModern(ref Model, tex_ids, false);

            // Show normals if enabled (must be done after model drawing)
            if (renderContext.Options.NormalsDisplayMode != NormalsDisplayMode.None)
            {
                for (int g = 0; g < Model.Header.numGroups; g++)
                {
                    ShowNormals(Model.Groups[g], Model.Polys, Model.Verts,
                                Model.Normals, Model.NormalIndex,
                                renderContext.Options.NormalsDisplayMode,
                                renderContext.Options.NormalsScale,
                                renderContext.Options.NormalsColor);
                }
            }

            // Restore original matrices
            GLRenderer.ProjectionMatrix = savedProjection;
            GLRenderer.ViewMatrix = savedView;
            GLRenderer.ModelMatrix = savedModel;
        }

        public static void DrawPModelBoundingBox(PModel Model)
        {
            GL.Begin(PrimitiveType.Lines);
            GL.Disable(EnableCap.DepthTest);

            DrawBox(Model.BoundingBox.max_x, Model.BoundingBox.max_y, Model.BoundingBox.max_z,
                    Model.BoundingBox.min_x, Model.BoundingBox.min_y, Model.BoundingBox.min_z,
                    1, 1, 0);

            GL.Enable(EnableCap.DepthTest);
            GL.End();
            // GL.PopMatrix();
        }

 

        //  ---------------------------------------------------------------------------------------------------
        //  ======================================== FIELD DRAW  ==============================================
        //  ---------------------------------------------------------------------------------------------------
        public static int MoveToFieldBone(FieldSkeleton fSkeleton, FieldFrame fFrame, int b_index)
        {
            int iBoneIdx, jsp;
            string[] joint_stack = new string[fSkeleton.bones.Count];

            GL.MatrixMode(MatrixMode.Modelview);

            jsp = 0;
            joint_stack[jsp] = fSkeleton.bones[0].joint_f;

            for (iBoneIdx = 0; iBoneIdx < b_index; iBoneIdx++)
            {
                while (!(fSkeleton.bones[iBoneIdx].joint_f == joint_stack[jsp]) && jsp > 0)
                {
                    GL.PopMatrix();
                    jsp--;
                }
                GL.PushMatrix();

                GL.Rotated(fFrame.rotations[iBoneIdx].beta, 0, 1, 0);
                GL.Rotated(fFrame.rotations[iBoneIdx].alpha, 1, 0, 0);
                GL.Rotated(fFrame.rotations[iBoneIdx].gamma, 0, 0, 1);

                GL.Translated(0, 0, -fSkeleton.bones[iBoneIdx].len);

                jsp++;
                joint_stack[jsp] = fSkeleton.bones[iBoneIdx].joint_i;
            }

            while (!(fSkeleton.bones[b_index].joint_f == joint_stack[jsp]) && jsp > 0)
            {
                GL.PopMatrix();
                jsp--;
            }
            GL.PushMatrix();

            GL.Rotated(fFrame.rotations[b_index].beta, 0, 1, 0);
            GL.Rotated(fFrame.rotations[b_index].alpha, 1, 0, 0);
            GL.Rotated(fFrame.rotations[b_index].gamma, 0, 0, 1);

            return jsp + 1;
        }

        // Cache for bounding box mesh
        private static LineMesh? _boundingBoxMesh;

        public static void DrawBox(float max_x, float max_y, float max_z,
                                   float min_x, float min_y, float min_z,
                                   float red, float green, float blue)
        {
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

            _boundingBoxMesh?.Dispose();
            _boundingBoxMesh = CreateBoundingBoxMesh(max_x, max_y, max_z, min_x, min_y, min_z, red, green, blue);
            GLRenderer.DrawLinesModern(_boundingBoxMesh);

            // Restore original matrices
            GLRenderer.ProjectionMatrix = savedProjection;
            GLRenderer.ViewMatrix = savedView;
            GLRenderer.ModelMatrix = savedModel;
        }

        public static void DrawFieldBonePieceBoundingBox(FieldBone bone, int p_index)
        {
            double[] rot_mat = new double[16];

            GL.Disable(EnableCap.DepthTest);
            GL.MatrixMode(MatrixMode.Modelview);
            GL.Scaled(bone.resizeX, bone.resizeY, bone.resizeZ);

            GL.Translated(bone.fRSDResources[p_index].Model.repositionX,
                          bone.fRSDResources[p_index].Model.repositionY,
                          bone.fRSDResources[p_index].Model.repositionZ);

            BuildMatrixFromQuaternion(bone.fRSDResources[p_index].Model.rotationQuaternion, ref rot_mat);

            GL.MultMatrixd(rot_mat);
            GL.Scaled(bone.fRSDResources[p_index].Model.resizeX,
                              bone.fRSDResources[p_index].Model.resizeY,
                              bone.fRSDResources[p_index].Model.resizeZ);

            DrawBox(bone.fRSDResources[p_index].Model.BoundingBox.max_x,
                    bone.fRSDResources[p_index].Model.BoundingBox.max_y,
                    bone.fRSDResources[p_index].Model.BoundingBox.max_z,
                    bone.fRSDResources[p_index].Model.BoundingBox.min_x,
                    bone.fRSDResources[p_index].Model.BoundingBox.min_y,
                    bone.fRSDResources[p_index].Model.BoundingBox.min_z,
                    0, 1, 0);

            GL.Enable(EnableCap.DepthTest);
        }

        public static void DrawFieldBoneBoundingBox(FieldBone bone)
        {
            int iResourceIdx;

            float max_x, max_y, max_z;
            float min_x, min_y, min_z;

            GL.MatrixMode(MatrixMode.Modelview);
            GL.Scaled(bone.resizeX, bone.resizeY, bone.resizeZ);

            if (bone.nResources == 0)
            {
                GL.Disable(EnableCap.DepthTest);
                
                GL.Color3f(1, 0, 0);
                GL.Begin(PrimitiveType.Lines);
                    GL.Vertex3f(0, 0, 0);
                    GL.Vertex3f(0, 0, -(float)bone.len);
                GL.End();

                GL.Enable(EnableCap.DepthTest);
            }
            else
            {
                max_x = float.NegativeInfinity;
                max_y = float.NegativeInfinity;
                max_z = float.NegativeInfinity;

                min_x = float.PositiveInfinity;
                min_y = float.PositiveInfinity;
                min_z = float.PositiveInfinity;

                for (iResourceIdx = 0; iResourceIdx < bone.nResources; iResourceIdx++)
                {
                    if (max_x < bone.fRSDResources[iResourceIdx].Model.BoundingBox.max_x) 
                            max_x = bone.fRSDResources[iResourceIdx].Model.BoundingBox.max_x;
                    if (max_y < bone.fRSDResources[iResourceIdx].Model.BoundingBox.max_y) 
                            max_y = bone.fRSDResources[iResourceIdx].Model.BoundingBox.max_y;
                    if (max_z < bone.fRSDResources[iResourceIdx].Model.BoundingBox.max_z) 
                            max_z = bone.fRSDResources[iResourceIdx].Model.BoundingBox.max_z;

                    if (min_x > bone.fRSDResources[iResourceIdx].Model.BoundingBox.min_x) 
                            min_x = bone.fRSDResources[iResourceIdx].Model.BoundingBox.min_x;
                    if (min_y > bone.fRSDResources[iResourceIdx].Model.BoundingBox.min_x) 
                            min_y = bone.fRSDResources[iResourceIdx].Model.BoundingBox.min_y;
                    if (min_z > bone.fRSDResources[iResourceIdx].Model.BoundingBox.min_x) 
                            min_z = bone.fRSDResources[iResourceIdx].Model.BoundingBox.min_z;
                }

                GL.Disable(EnableCap.DepthTest);
                DrawBox(max_x, max_y, max_z, min_x, min_y, min_z, 1, 0, 0);
                GL.Enable(EnableCap.DepthTest);
            }
        }

        public static void DrawFieldSkeletonBones(FieldSkeleton fSkeleton, FieldFrame fFrame)
        {
            int iBoneIdx, jsp;
            string[] joint_stack = new string[fSkeleton.bones.Count + 1];
            double[] rot_mat = new double[16];

            jsp = 0;

            GL.MatrixMode(MatrixMode.Modelview);
            GL.PushMatrix();

            GL.Translated(fFrame.rootTranslationX, 0, 0);
            GL.Translated(0, -fFrame.rootTranslationY, 0);
            GL.Translated(0, 0, fFrame.rootTranslationZ);

            BuildRotationMatrixWithQuaternions(fFrame.rootRotationAlpha, fFrame.rootRotationBeta, fFrame.rootRotationGamma, ref rot_mat);

            GL.MultMatrixd(rot_mat);
            GL.PointSize(5f);

            joint_stack[jsp] = fSkeleton.bones[0].joint_f;

            for (iBoneIdx = 0; iBoneIdx < fSkeleton.bones.Count; iBoneIdx++)
            {
                while ((fSkeleton.bones[iBoneIdx].joint_f != joint_stack[jsp]) && jsp > 0)
                {
                    GL.PopMatrix();
                    jsp--;
                }

                GL.PushMatrix();

                // -- Commented in KimeraVB6
                //GL.Rotated(fFrame.rotations[bi].beta, 0, 1, 0);
                //GL.Rotated(fFrame.rotations[bi].alpha, 1, 0, 0);
                //GL.Rotated(fFrame.rotations[bi].gamma, 0, 0, 1);
                BuildRotationMatrixWithQuaternions(fFrame.rotations[iBoneIdx].alpha, 
                                                   fFrame.rotations[iBoneIdx].beta, 
                                                   fFrame.rotations[iBoneIdx].gamma, 
                                                   ref rot_mat);
                GL.MultMatrixd(rot_mat);

                GL.Begin(PrimitiveType.Points);
                    GL.Vertex3f(0, 0, 0);
                    GL.Vertex3f(0, 0, (float)-fSkeleton.bones[iBoneIdx].len);
                GL.End();

                GL.Begin(PrimitiveType.Lines);
                    GL.Vertex3f(0, 0, 0);
                    GL.Vertex3f(0, 0, (float)-fSkeleton.bones[iBoneIdx].len);
                GL.End();

                GL.Translated(0, 0, -fSkeleton.bones[iBoneIdx].len);

                jsp++;
                joint_stack[jsp] = fSkeleton.bones[iBoneIdx].joint_i;
            }

            while (jsp > 0)
            {
                GL.PopMatrix();
                jsp--;
            }
            GL.PopMatrix();
        }

        public static void DrawRSDResource(FieldRSDResource fRSDResource, bool bDListsEnable,
                                           RenderingContext renderContext)
        {
            int iTextureIdx;
            uint[] tex_ids;
            double[] rot_mat = new double[16];

            tex_ids = new uint[fRSDResource.numTextures];

            for (iTextureIdx = 0; iTextureIdx < fRSDResource.numTextures; iTextureIdx++)
            {
                tex_ids[iTextureIdx] = fRSDResource.textures[iTextureIdx].texID;
            }

            GL.MatrixMode(MatrixMode.Modelview);
            GL.PushMatrix();

            GL.Translated(fRSDResource.Model.repositionX, fRSDResource.Model.repositionY, fRSDResource.Model.repositionZ);
            BuildMatrixFromQuaternion(fRSDResource.Model.rotationQuaternion, ref rot_mat);

            GL.MultMatrixd(rot_mat);

            GL.Scaled(fRSDResource.Model.resizeX, fRSDResource.Model.resizeY, fRSDResource.Model.resizeZ);

            SetDefaultOGLRenderState();

            if (!bDListsEnable) DrawPModel(ref fRSDResource.Model, ref tex_ids, false, renderContext);
            else DrawPModelDLists(ref fRSDResource.Model, ref tex_ids, renderContext);

            GL.PopMatrix();
        }

        public static void DrawFieldBone(FieldBone bone, bool bDListsEnable, RenderingContext renderContext)
        {

            int iResourceIdx;

            GL.MatrixMode(MatrixMode.Modelview);
            GL.PushMatrix();

            GL.Scaled(bone.resizeX, bone.resizeY, bone.resizeZ);

            for (iResourceIdx = 0; iResourceIdx < bone.nResources; iResourceIdx++)
                DrawRSDResource(bone.fRSDResources[iResourceIdx], bDListsEnable,
                                renderContext);

            GL.PopMatrix();
        }

        public static void DrawFieldSkeleton(FieldSkeleton fSkeleton, FieldFrame fFrame,
                                             bool bDListsEnable, RenderingContext renderContext)
        {
            int iBoneIdx;
            string[] joint_stack = new string[fSkeleton.bones.Count + 1];
            int jsp;
            double[] rot_mat = new double[16];

            GL.MatrixMode(MatrixMode.Modelview);

            GL.PushMatrix();
            GL.Translated(fFrame.rootTranslationX, 0, 0);
            GL.Translated(0, -fFrame.rootTranslationY, 0);
            GL.Translated(0, 0, fFrame.rootTranslationZ);

            BuildRotationMatrixWithQuaternions(fFrame.rootRotationAlpha, fFrame.rootRotationBeta, fFrame.rootRotationGamma, ref rot_mat);

            GL.MultMatrixd(rot_mat);

            jsp = 0;
            joint_stack[jsp] = fSkeleton.bones[0].joint_f;

            //for (bi = 0; bi < Skeleton.nBones; bi++)
            for (iBoneIdx = 0; iBoneIdx < fSkeleton.bones.Count; iBoneIdx++)
            {
                while (!(fSkeleton.bones[iBoneIdx].joint_f == joint_stack[jsp]) && jsp > 0)
                {
                    GL.PopMatrix();
                    jsp--;
                }

                //if (jsp == 0) SetDefaultOGLRenderState();

                GL.PushMatrix();

                BuildRotationMatrixWithQuaternions(fFrame.rotations[iBoneIdx].alpha, 
                                                   fFrame.rotations[iBoneIdx].beta, 
                                                   fFrame.rotations[iBoneIdx].gamma, 
                                                   ref rot_mat);

                GL.MultMatrixd(rot_mat);

                DrawFieldBone(fSkeleton.bones[iBoneIdx], bDListsEnable, renderContext);

                GL.Translated(0, 0, -fSkeleton.bones[iBoneIdx].len);

                jsp++;
                joint_stack[jsp] = fSkeleton.bones[iBoneIdx].joint_i;
            }

            while (jsp > 0)
            {
                GL.PopMatrix();
                jsp--;
            }
            GL.PopMatrix();
        }


        //  ---------------------------------------------------------------------------------------------------
        //  ======================================== BATTLE DRAW  =============================================
        //  ---------------------------------------------------------------------------------------------------
        public static int MoveToBattleBone(BattleSkeleton bSkeleton, BattleFrame bFrame, int boneIndex)
        {
            int iBoneIdx, jsp, itmpbones;
            int[] joint_stack = new int[bSkeleton.nBones * 4];
            double[] rot_mat = new double[16];

            jsp = 0;
            joint_stack[jsp] = -1;

            if (bSkeleton.nBones > 1) itmpbones = 1;
            else itmpbones = 0;

            GL.MatrixMode(MatrixMode.Modelview);
            GL.PushMatrix();
            GL.Translated(bFrame.startX, bFrame.startY, bFrame.startZ);

            BuildRotationMatrixWithQuaternions(bFrame.bones[0].alpha, bFrame.bones[0].beta, bFrame.bones[0].gamma, ref rot_mat);
            GL.MultMatrixd(rot_mat);

            for (iBoneIdx = 0; iBoneIdx < boneIndex; iBoneIdx++)
            {
                GL.PushName((uint)iBoneIdx);

                while (!(bSkeleton.bones[iBoneIdx].parentBone == joint_stack[jsp]) && jsp > 0)
                {
                    GL.PopMatrix();
                    jsp--;
                }
                GL.PushMatrix();

                // -- Commented in KimeraVB6
                //  GL.Rotated(bFrame.bones[bi + 1].beta, 0, 1, 0);
                //  GL.Rotated(bFrame.bones[bi + 1].alpha, 1, 0, 0);
                //  GL.Rotated(bFrame.bones[bi + 1].gamma, 0, 0, 1);

                BuildRotationMatrixWithQuaternions(bFrame.bones[iBoneIdx + itmpbones].alpha, 
                                                       bFrame.bones[iBoneIdx + itmpbones].beta, 
                                                       bFrame.bones[iBoneIdx + itmpbones].gamma, ref rot_mat);
                GL.MultMatrixd(rot_mat);

                GL.Translated(0, 0, bSkeleton.bones[iBoneIdx].len);

                jsp++;
                joint_stack[jsp] = iBoneIdx;

                GL.PopName();
            }

            while (!(bSkeleton.bones[iBoneIdx].parentBone == joint_stack[jsp]) && jsp > 0)
            {
                GL.PopMatrix();
                jsp--;
            }

            // -- Commented in KimeraVB6
            // GL.PopMatrix();
            //  GL.Rotated(bFrame.bones[boneIndex + itmpbones].beta, 0, 1, 0);
            //  GL.Rotated(bFrame.bones[boneIndex + itmpbones].alpha, 1, 0, 0);
            //  GL.Rotated(bFrame.bones[boneIndex + itmpbones].gamma, 0, 0, 1);

            BuildRotationMatrixWithQuaternions(bFrame.bones[boneIndex + itmpbones].alpha,
                                               bFrame.bones[boneIndex + itmpbones].beta,
                                               bFrame.bones[boneIndex + itmpbones].gamma,
                                               ref rot_mat);
            GL.MultMatrixd(rot_mat);

            return jsp + 1;
        }

        public static int MoveToBattleBoneMiddle(BattleSkeleton bSkeleton, BattleFrame bFrame, int boneIndex)
        {
            int iMoveToBattleBoneMiddleResult;

            iMoveToBattleBoneMiddleResult = MoveToBattleBone(bSkeleton, bFrame, boneIndex);
            GL.Translated(0, 0, bSkeleton.bones[boneIndex].len / 2);

            return iMoveToBattleBoneMiddleResult;
        }

        public static int MoveToBattleBoneEnd(BattleSkeleton bSkeleton, BattleFrame bFrame, int boneIndex)
        {
            int iMoveToBattleBoneEndResult;

            iMoveToBattleBoneEndResult = MoveToBattleBone(bSkeleton, bFrame, boneIndex);
            GL.Translated(0, 0, bSkeleton.bones[boneIndex].len);

            return iMoveToBattleBoneEndResult;
        }

        public static void DrawBattleBoneModelBoundingBox(BattleBone bBone, int partIndex)
        {
            double[] rot_mat = new double[16];

            GL.Disable(EnableCap.DepthTest);
            GL.MatrixMode(MatrixMode.Modelview);
            GL.Scaled(bBone.resizeX, bBone.resizeY, bBone.resizeZ);

            GL.Translated(bBone.Models[partIndex].repositionX, bBone.Models[partIndex].repositionY, bBone.Models[partIndex].repositionZ);

            BuildMatrixFromQuaternion(bBone.Models[partIndex].rotationQuaternion, ref rot_mat);
            GL.MultMatrixd(rot_mat);

            GL.Scaled(bBone.Models[partIndex].resizeX, bBone.Models[partIndex].resizeY, bBone.Models[partIndex].resizeZ);

            DrawBox(bBone.Models[partIndex].BoundingBox.max_x, bBone.Models[partIndex].BoundingBox.max_y, bBone.Models[partIndex].BoundingBox.max_z,
                    bBone.Models[partIndex].BoundingBox.min_x, bBone.Models[partIndex].BoundingBox.min_y, bBone.Models[partIndex].BoundingBox.min_z,
                    0, 1, 0);
            GL.Enable(EnableCap.DepthTest);
        }

        public static void DrawBattleBoneBoundingBox(BattleBone bBone)
        {

            GL.Disable(EnableCap.DepthTest);
            GL.MatrixMode(MatrixMode.Modelview);

            GL.Scaled(bBone.resizeX, bBone.resizeY, bBone.resizeZ);

            if (bBone.hasModel == 1)
            {
                // -- Commented in KimeraVB6
                //GL.Translated(bBone.Models[0].repositionX, bBone.Models[0].repositionY, bBone.Models[0].repositionZ);
                //BuildMatrixFromQuaternion(ref bBone.Models[0].rotationQuaternion, ref rot_mat);
                //GL.MultMatrixd(rot_mat);
                //GL.Scaled(bBone.Models[0].resizeX, bBone.Models[0].resizeY, bBone.Models[0].resizeZ);

                DrawBox(bBone.Models[0].BoundingBox.max_x, bBone.Models[0].BoundingBox.max_y, bBone.Models[0].BoundingBox.max_z,
                        bBone.Models[0].BoundingBox.min_x, bBone.Models[0].BoundingBox.min_y, bBone.Models[0].BoundingBox.min_z,
                        0, 1, 0);
                GL.Enable(EnableCap.DepthTest);
            }
            else
            {
                GL.Color3f(0, 1, 0);
                GL.Begin(PrimitiveType.Lines);
                    GL.Vertex3f(0, 0, 0);
                    GL.Vertex3f(0, 0, bBone.len);
                GL.End();
            }

            GL.Enable(EnableCap.DepthTest);
        }

        public static void DrawBattleWeaponBoundingBox(BattleSkeleton bSkeleton, BattleFrame wpFrame,
                                                       int weaponIndex, int animWeaponIndex)
        {
            double[] rot_mat = new double[16];

            //if (weaponIndex > -1 && bSkeleton.nWeapons > 0)       // -- Commented in KimeraVB6
            if (animWeaponIndex > -1 && bSkeleton.wpModels.Count > 0 && bAnimationsPack.WeaponAnimations.Count > 0)
            {
                GL.PushMatrix();
                GL.Translated(wpFrame.startX, wpFrame.startY, wpFrame.startZ);

                // -- Commented in KimeraVB6
                //GL.Rotated(wpFrame.bones[0].beta, 0, 1, 0);
                //GL.Rotated(wpFrame.bones[0].alpha, 1, 0, 0);
                //GL.Rotated(wpFrame.bones[0].gamma, 0, 0, 1);

                BuildRotationMatrixWithQuaternions(wpFrame.bones[0].alpha, wpFrame.bones[0].beta, wpFrame.bones[0].gamma, ref rot_mat);
                GL.MultMatrixd(rot_mat);

                GL.PushMatrix();

                GL.Translated(bSkeleton.wpModels[weaponIndex].repositionX,
                             bSkeleton.wpModels[weaponIndex].repositionY,
                             bSkeleton.wpModels[weaponIndex].repositionZ);
                
                GL.Rotated(bSkeleton.wpModels[weaponIndex].rotateBeta, 0, 1, 0);
                GL.Rotated(bSkeleton.wpModels[weaponIndex].rotateAlpha, 1, 0, 0);
                GL.Rotated(bSkeleton.wpModels[weaponIndex].rotateGamma, 0, 0, 1);
                
                GL.Scaled(bSkeleton.wpModels[weaponIndex].resizeX, bSkeleton.wpModels[weaponIndex].resizeY, bSkeleton.wpModels[weaponIndex].resizeZ);

                DrawPModelBoundingBox(bSkeleton.wpModels[weaponIndex]);

                GL.PopMatrix();
                GL.PopMatrix();
            }
        }

        public static void SelectBattleBoneAndModel(BattleSkeleton bSkeleton, BattleFrame bFrame, BattleFrame wpFrame,
                                                    int weaponIndex, int boneIndex, int partIndex)
        {
            int i, jsp;

            if (boneIndex > -1 && boneIndex < bSkeleton.nBones)
            {
                jsp = MoveToBattleBone(bSkeleton, bFrame, boneIndex);
                DrawBattleBoneBoundingBox(bSkeleton.bones[boneIndex]);

                if (partIndex > -1)
                    DrawBattleBoneModelBoundingBox(bSkeleton.bones[boneIndex], partIndex);

                for (i = 0; i <= jsp; i++) GL.PopMatrix();
            }
            else
            {
                if (boneIndex == bSkeleton.nBones)
                    DrawBattleWeaponBoundingBox(bSkeleton, wpFrame, weaponIndex, weaponIndex);
            }
        }

        public static void DrawBattleSkeletonBones(BattleSkeleton bSkeleton, BattleFrame bFrame)
        {
            int iBoneIdx, jsp, itmpbones;
            int[] joint_stack;
            double[] rot_mat = new double[16];

            if (bSkeleton.IsBattleLocation) return;

            joint_stack = new int[bSkeleton.nBones + 1];
            jsp = 0;
            joint_stack[jsp] = -1;

            if (bSkeleton.nBones > 1) itmpbones = 1;
            else itmpbones = 0;

            GL.MatrixMode(MatrixMode.Modelview);
            GL.PointSize(5);
            GL.PushMatrix();
            GL.Translated(bFrame.startX, bFrame.startY, bFrame.startZ);

            BuildRotationMatrixWithQuaternions(bFrame.bones[0].alpha, bFrame.bones[0].beta, bFrame.bones[0].gamma, ref rot_mat);
            GL.MultMatrixd(rot_mat);

            for (iBoneIdx = 0; iBoneIdx < bSkeleton.nBones; iBoneIdx++)
            {
                while (!(bSkeleton.bones[iBoneIdx].parentBone == joint_stack[jsp]) && jsp > 0)
                {
                    GL.PopMatrix();
                    jsp--;
                }
                GL.PushMatrix();

                // -- Commented in KimeraVB6
                //GL.Rotated(bFrame.bones[bi + 1].beta, 0, 1, 0);
                //GL.Rotated(bFrame.bones[bi + 1].alpha, 1, 0, 0);
                //GL.Rotated(bFrame.bones[bi + 1].gamma, 0, 0, 1);

                BuildRotationMatrixWithQuaternions(bFrame.bones[iBoneIdx + itmpbones].alpha,
                                                   bFrame.bones[iBoneIdx + itmpbones].beta,
                                                   bFrame.bones[iBoneIdx + itmpbones].gamma,
                                                   ref rot_mat);
                GL.MultMatrixd(rot_mat);

                GL.Begin(PrimitiveType.Points);
                    GL.Vertex3f(0, 0, 0);
                    GL.Vertex3f(0, 0, bSkeleton.bones[iBoneIdx].len);
                GL.End();

                GL.Begin(PrimitiveType.Lines);
                    GL.Vertex3f(0, 0, 0);
                    GL.Vertex3f(0, 0, bSkeleton.bones[iBoneIdx].len);
                GL.End();

                GL.Translated(0, 0, bSkeleton.bones[iBoneIdx].len);

                jsp++;
                joint_stack[jsp] = iBoneIdx;
            }

            if (!bSkeleton.IsBattleLocation)
            {
                while (jsp > 0)
                {
                    GL.PopMatrix();
                    jsp--;
                }
            }
            GL.PopMatrix();
        }

        public static void DrawBattleSkeletonBone(BattleBone bBone, uint[] texIDS, bool bDListsEnable,
                                                  RenderingContext renderContext)
        {
            int iModelIdx;
            PModel tmpbModel = new PModel();

            GL.MatrixMode(MatrixMode.Modelview);
            GL.PushMatrix();
            GL.Scaled(bBone.resizeX, bBone.resizeY, bBone.resizeZ);

            if (bBone.hasModel > 0)
            {

                if (!bDListsEnable)
                {
                    for (iModelIdx = 0; iModelIdx < bBone.nModels; iModelIdx++)
                    {

                        GL.PushMatrix();
                        GL.Translated(bBone.Models[iModelIdx].repositionX, 
                                     bBone.Models[iModelIdx].repositionY, 
                                     bBone.Models[iModelIdx].repositionZ);

                        GL.Rotated(bBone.Models[iModelIdx].rotateAlpha, 1, 0, 0);
                        GL.Rotated(bBone.Models[iModelIdx].rotateBeta, 0, 1, 0);
                        GL.Rotated(bBone.Models[iModelIdx].rotateGamma, 0, 0, 1);

                        GL.Scaled(bBone.Models[iModelIdx].resizeX, 
                                 bBone.Models[iModelIdx].resizeY, 
                                 bBone.Models[iModelIdx].resizeZ);
                        
                        SetDefaultOGLRenderState();

                        tmpbModel = bBone.Models[iModelIdx];
                        DrawPModel(ref tmpbModel, ref texIDS, false, renderContext);
                        bBone.Models[iModelIdx] = tmpbModel;

                        GL.PopMatrix();
                    }
                }
                else
                {
                    for (iModelIdx = 0; iModelIdx < bBone.nModels; iModelIdx++)
                    {
                        GL.PushMatrix();
                        GL.Translated(bBone.Models[iModelIdx].repositionX, 
                                     bBone.Models[iModelIdx].repositionY, 
                                     bBone.Models[iModelIdx].repositionZ);

                        GL.Rotated(bBone.Models[iModelIdx].rotateAlpha, 1, 0, 0);
                        GL.Rotated(bBone.Models[iModelIdx].rotateBeta, 0, 1, 0);
                        GL.Rotated(bBone.Models[iModelIdx].rotateGamma, 0, 0, 1);

                        GL.Scaled(bBone.Models[iModelIdx].resizeX, 
                                 bBone.Models[iModelIdx].resizeY, 
                                 bBone.Models[iModelIdx].resizeZ);

                        tmpbModel = bBone.Models[iModelIdx];
                        DrawPModelDLists(ref tmpbModel, ref texIDS, renderContext);
                        bBone.Models[iModelIdx] = tmpbModel;

                        GL.PopMatrix();
                    }
                }
            }

            GL.PopMatrix();
        }

        public static void DrawBattleSkeleton(BattleSkeleton bSkeleton, BattleFrame bFrame, BattleFrame wpFrame,
                                              int weaponIndex, int animWeaponIndex, bool bDListsEnable,
                                              RenderingContext renderContext)
        {
            int iBoneIdx, jsp, itmpbones;
            int[] joint_stack = new int[bSkeleton.nBones + 1];
            double[] rot_mat = new double[16];

            jsp = 0;
            joint_stack[jsp] = -1;

            if (bSkeleton.nBones > 1) itmpbones = 1;
            else itmpbones = 0;

            GL.MatrixMode(MatrixMode.Modelview);
            GL.PushMatrix();
            GL.Translated(bFrame.startX, bFrame.startY, bFrame.startZ);

            // Debug.Print bFrame.bones[0].alpha; ", "; bFrame.bones[0].Beta; ", "; bFrame.bones[0].Gamma
            BuildRotationMatrixWithQuaternions(bFrame.bones[0].alpha, bFrame.bones[0].beta, bFrame.bones[0].gamma, ref rot_mat);
            GL.MultMatrixd(rot_mat);

            for (iBoneIdx = 0; iBoneIdx < bSkeleton.nBones; iBoneIdx++)
            {
                if (bSkeleton.IsBattleLocation)
                {
                    DrawBattleSkeletonBone(bSkeleton.bones[iBoneIdx], bSkeleton.TexIDS, false, renderContext);
                }
                else
                {
                    while (!(bSkeleton.bones[iBoneIdx].parentBone == joint_stack[jsp]) && jsp > 0)
                    {
                        GL.PopMatrix();
                        jsp--;
                    }

                    GL.PushMatrix();

                    // -- Commented in KimeraVB6
                    //GL.Rotated(bFrame.bones[bi + 1].beta, 0, 1, 0);
                    //GL.Rotated(bFrame.bones[bi + 1].alpha, 1, 0, 0);
                    //GL.Rotated(bFrame.bones[bi + 1].gamma, 0, 0, 1);

                    BuildRotationMatrixWithQuaternions(bFrame.bones[iBoneIdx + itmpbones].alpha,
                                                       bFrame.bones[iBoneIdx + itmpbones].beta,
                                                       bFrame.bones[iBoneIdx + itmpbones].gamma,
                                                       ref rot_mat);
                    GL.MultMatrixd(rot_mat);

                    DrawBattleSkeletonBone(bSkeleton.bones[iBoneIdx], bSkeleton.TexIDS, bDListsEnable,
                                           renderContext);

                    GL.Translated(0, 0, bSkeleton.bones[iBoneIdx].len);

                    jsp++;
                    joint_stack[jsp] = iBoneIdx;
                }
            }

            if (!bSkeleton.IsBattleLocation)
            {
                while (jsp > 0)
                {
                    GL.PopMatrix();
                    jsp--;
                }
            }
            GL.PopMatrix();

            //if (weaponIndex > -1 && bSkeleton.nWeapons > 0)       // -- Commented in KimeraVB6
            if (animWeaponIndex > -1 && bSkeleton.wpModels.Count > 0 && bAnimationsPack.WeaponAnimations.Count > 0)
            {
                GL.PushMatrix();
                GL.Translated(wpFrame.startX, wpFrame.startY, wpFrame.startZ);

                // -- Commented in KimeraVB6
                //GL.Rotated(wpFrame.bones[0].beta, 0, 1, 0);
                //GL.Rotated(wpFrame.bones[0].alpha, 1, 0, 0);
                //GL.Rotated(wpFrame.bones[0].gamma, 0, 0, 1);

                BuildRotationMatrixWithQuaternions(wpFrame.bones[0].alpha, wpFrame.bones[0].beta, wpFrame.bones[0].gamma, ref rot_mat);
                GL.MultMatrixd(rot_mat);

                GL.MatrixMode(MatrixMode.Modelview);
                GL.PushMatrix();

                GL.Translated(bSkeleton.wpModels[weaponIndex].repositionX,
                             bSkeleton.wpModels[weaponIndex].repositionY,
                             bSkeleton.wpModels[weaponIndex].repositionZ);

                GL.Rotated(bSkeleton.wpModels[weaponIndex].rotateAlpha, 1, 0, 0);
                GL.Rotated(bSkeleton.wpModels[weaponIndex].rotateBeta, 0, 1, 0);
                GL.Rotated(bSkeleton.wpModels[weaponIndex].rotateGamma, 0, 0, 1);

                GL.Scaled(bSkeleton.wpModels[weaponIndex].resizeX, bSkeleton.wpModels[weaponIndex].resizeY, bSkeleton.wpModels[weaponIndex].resizeZ);

                SetDefaultOGLRenderState();

                PModel tmpwpModel = new PModel();
                if (bDListsEnable)
                {
                    tmpwpModel = bSkeleton.wpModels[weaponIndex];
                    DrawPModelDLists(ref tmpwpModel, ref bSkeleton.TexIDS, renderContext);
                    bSkeleton.wpModels[weaponIndex] = tmpwpModel;
                }
                else
                {
                    tmpwpModel = bSkeleton.wpModels[weaponIndex];
                    DrawPModel(ref tmpwpModel, ref bSkeleton.TexIDS, false, renderContext);
                    bSkeleton.wpModels[weaponIndex] = tmpwpModel;
                }
                GL.PopMatrix();

                GL.PopMatrix();
            }
        }



        //  ---------------------------------------------------------------------------------------------------
        //  ======================================= SKELETON DRAW  ============================================
        //  ---------------------------------------------------------------------------------------------------

        /// <summary>
        /// Draws the current skeleton/model using the provided context.
        /// If ctx.ModelData is set, uses that model data (fully decoupled).
        /// Otherwise falls back to static model data for backward compatibility.
        /// </summary>
        public static void DrawSkeletonModel(RenderingContext ctx)
        {
            double[] rot_mat = new double[16];
            BattleFrame tmpbFrame;

            Vector3 p_min = new Vector3();
            Vector3 p_max = new Vector3();

            // Use context model data if provided, otherwise fall back to statics
            var modelData = ctx.ModelData;
            PModel pModel = modelData?.PModel ?? fPModel;
            uint[] texIds = modelData?.TextureIds ?? tex_ids;
            FieldSkeleton fieldSkel = modelData?.FieldSkeleton ?? fSkeleton;
            FieldAnimation fieldAnim = modelData?.FieldAnimation ?? fAnimation;
            BattleSkeleton battleSkel = modelData?.BattleSkeleton ?? bSkeleton;
            BattleAnimationsPack battleAnims = modelData?.BattleAnimations ?? bAnimationsPack;

            try
            {
                switch (ctx.ModelType)
                {
                    case K_3DS_MODEL:
                    case K_P_FIELD_MODEL:
                    case K_P_BATTLE_MODEL:
                    case K_P_MAGIC_MODEL:
                        ComputePModelBoundingBox(pModel, ref p_min, ref p_max);

                        SetCameraAroundModel(ref p_min, ref p_max,
                                             ctx.Camera.PanX, ctx.Camera.PanY, ctx.Camera.PanZ + ctx.Camera.Distance,
                                             ctx.Camera.Alpha, ctx.Camera.Beta, ctx.Camera.Gamma, 1, 1, 1);

                        if (ctx.Options.ShowGround)
                        {
                            GL.Disable(EnableCap.Lighting);
                            DrawGround();
                            DrawShadow(ref p_min, ref p_max);
                        }

                        SetLights(ctx.Lighting, (float)(-2 * ComputeSceneRadius(p_min, p_max)));

                        GL.MatrixMode(MatrixMode.Modelview);
                        GL.PushMatrix();

                        GL.Translated(pModel.repositionX,
                                     pModel.repositionY,
                                     pModel.repositionZ);

                        BuildRotationMatrixWithQuaternionsXYZ(pModel.rotateAlpha,
                                                              pModel.rotateBeta,
                                                              pModel.rotateGamma,
                                                              ref rot_mat);

                        GL.MultMatrixd(rot_mat);
                        GL.Scaled(pModel.resizeX,
                                 pModel.resizeY,
                                 pModel.resizeZ);

                        DrawPModel(ref pModel, ref texIds, false, ctx);

                        GL.PopMatrix();

                        break;

                    case K_HRC_SKELETON:
                        ComputeFieldBoundingBox(fieldSkel, fieldAnim.frames[ctx.Animation.CurrentFrame],
                                                ref p_min, ref p_max);

                        SetCameraAroundModel(ref p_min, ref p_max,
                                             ctx.Camera.PanX, ctx.Camera.PanY, ctx.Camera.PanZ + ctx.Camera.Distance,
                                             ctx.Camera.Alpha, ctx.Camera.Beta, ctx.Camera.Gamma, 1, 1, 1);

                        if (ctx.Options.ShowGround)
                        {
                            GL.Disable(EnableCap.Lighting);
                            DrawGround();
                            DrawShadow(ref p_min, ref p_max);
                        }

                        SetLights(ctx.Lighting, (float)(-2 * ComputeSceneRadius(p_min, p_max)));

                        DrawFieldSkeleton(fieldSkel, fieldAnim.frames[ctx.Animation.CurrentFrame],
                                          ctx.Options.EnableDisplayLists, ctx);

                        if (ctx.Options.ShowLastFrameGhost)
                        {

                            GL.ColorMask(true, true, false, true);
                            if (ctx.Animation.CurrentFrame == 0)
                                DrawFieldSkeleton(fieldSkel, fieldAnim.frames[fieldAnim.nFrames - 1],
                                                  ctx.Options.EnableDisplayLists, ctx);
                            else
                                DrawFieldSkeleton(fieldSkel, fieldAnim.frames[ctx.Animation.CurrentFrame - 1],
                                                  ctx.Options.EnableDisplayLists, ctx);

                            GL.ColorMask(true, true, true, true);
                        }

                        GL.Disable(EnableCap.Lighting);

                        if (ctx.Options.ShowBones)
                        {
                            GL.Disable(EnableCap.DepthTest);

                            SkeletonRenderer.RenderFieldSkeletonBones(fieldSkel,
                                fieldAnim.frames[ctx.Animation.CurrentFrame], 0, 1, 0, 1, 0, 0);

                            GL.Enable(EnableCap.DepthTest);
                        }

                        SelectFieldBoneAndPiece(fieldSkel, fieldAnim.frames[ctx.Animation.CurrentFrame],
                                                ctx.Selection.SelectedBone, ctx.Selection.SelectedBonePiece);
                        break;

                    case K_AA_SKELETON:
                    case K_MAGIC_SKELETON:
                        ComputeBattleBoundingBox(battleSkel, battleAnims.SkeletonAnimations[ctx.Animation.AnimationIndex].frames[ctx.Animation.CurrentFrame],
                                                 ref p_min, ref p_max);

                        SetCameraAroundModel(ref p_min, ref p_max,
                                             ctx.Camera.PanX, ctx.Camera.PanY, ctx.Camera.PanZ + ctx.Camera.Distance,
                                             ctx.Camera.Alpha, ctx.Camera.Beta, ctx.Camera.Gamma, 1, 1, 1);

                        if (ctx.Options.ShowGround)
                        {
                            GL.Disable(EnableCap.Lighting);
                            DrawGround();
                            DrawShadow(ref p_min, ref p_max);
                        }

                        SetLights(ctx.Lighting, (float)(-2 * ComputeSceneRadius(p_min, p_max)));

                        tmpbFrame = new BattleFrame();
                        if (battleSkel.wpModels.Count > 0 && battleAnims.WeaponAnimations.Count > 0)
                        {
                            tmpbFrame = battleAnims.WeaponAnimations[ctx.Animation.AnimationIndex].frames[ctx.Animation.CurrentFrame];
                        }

                        DrawBattleSkeleton(battleSkel, battleAnims.SkeletonAnimations[ctx.Animation.AnimationIndex].frames[ctx.Animation.CurrentFrame],
                                           tmpbFrame, 0, ctx.Animation.WeaponAnimationIndex, ctx.Options.EnableDisplayLists,
                                           ctx);

                        if (ctx.Options.ShowLastFrameGhost && !battleSkel.IsBattleLocation)
                        {
                            GL.ColorMask(true, true, false, true);

                            if (ctx.Animation.CurrentFrame == 0)
                            {
                                DrawBattleSkeleton(battleSkel, battleAnims.SkeletonAnimations[ctx.Animation.AnimationIndex].frames[ctx.Animation.CurrentFrame],
                                                   tmpbFrame, 0, ctx.Animation.WeaponAnimationIndex, ctx.Options.EnableDisplayLists,
                                                   ctx);
                            }
                            else
                                DrawBattleSkeleton(battleSkel, battleAnims.SkeletonAnimations[ctx.Animation.AnimationIndex].frames[ctx.Animation.CurrentFrame - 1],
                                                   tmpbFrame, 0, ctx.Animation.WeaponAnimationIndex, ctx.Options.EnableDisplayLists,
                                                   ctx);

                            GL.ColorMask(true, true, true, true);
                        }

                        GL.Disable(EnableCap.Lighting);

                        if (ctx.Options.ShowBones)
                        {
                            GL.Disable(EnableCap.DepthTest);

                            SkeletonRenderer.RenderBattleSkeletonBones(battleSkel,
                                battleAnims.SkeletonAnimations[ctx.Animation.AnimationIndex].frames[ctx.Animation.CurrentFrame],
                                0, 1, 0, 1, 0, 0);

                            GL.Enable(EnableCap.DepthTest);
                        }

                        SelectBattleBoneAndModel(battleSkel, battleAnims.SkeletonAnimations[ctx.Animation.AnimationIndex].frames[ctx.Animation.CurrentFrame],
                            tmpbFrame, ctx.Animation.WeaponAnimationIndex, ctx.Selection.SelectedBone, ctx.Selection.SelectedBonePiece);
                        break;
                }
            }
            catch (Exception ex)
            {
                throw new KimeraException("Error Drawing current model.", ex);
            }
        }



        //  ---------------------------------------------------------------------------------------------------
        //  ============================ DRAWING EXTENSIONS (GROUND/SHADOW)  ==================================
        //  ---------------------------------------------------------------------------------------------------
        public static void DrawGround()
        {
            int gi, lw;
            float r, g, b;

            //  Draw plane
            GL.Color3f(0.9f, 0.9f, 1f);
            GL.Disable(EnableCap.DepthTest);
            
            GL.Begin(PrimitiveType.Quads);
                GL.Vertex4f(300f, 0f, 300f, 0.001f);
                GL.Vertex4f(300f, 0f, -300f, 0.001f);
                GL.Vertex4f(-300f, 0f, -300f, 0.001f);
                GL.Vertex4f(-300f, 0f, 300f, 0.001f);
            GL.End();

            r = 0.9f;
            g = 0.9f;
            b = 1.0f;
            lw = 5;
            //glEnable(GLCapability.GL_LINE_SMOOTH);

            for (gi = 10; gi >= 5; gi--)
            {
                GL.LineWidth(lw);
                GL.Color3f(r, g, b);
                GL.Begin(PrimitiveType.Lines);
                    GL.Vertex4f(0f, 0f, 50f, 0.0001f);
                    GL.Vertex4f(0f, 0f, -50f, 0.0001f);
                    GL.Vertex4f(-50f, 0f, 0f, 0.0001f);
                    GL.Vertex4f(50f, 0f, 0f, 0.0001f);
                GL.End();

                r = 0.9f - 0.9f / 10f * (10 - gi);
                g = 0.9f - 0.9f / 10f * (10 - gi);
                b = 1 - 1f / 10f * (10 - gi);
                lw--;
            }

            GL.LineWidth(1f);
            //glDisable(GLCapability.GL_LINE_SMOOTH);
        }

        public static void DrawShadow(ref Vector3 p_min, ref Vector3 p_max)
        {
            float ground_radius, sub_y, cx, cz;
            int numSegments, si;

            Vector3 p_min_aux, p_max_aux;

            sub_y = p_max.Y;
            p_min_aux = p_min;
            p_max_aux = p_max;
            p_min_aux.Y = 0;
            p_max_aux.Y = 0;

            cx = (p_min.X + p_max.X) / 2;
            cz = (p_min.Z + p_max.Z) / 2;
            ground_radius = CalculateDistance(p_min_aux, p_max_aux) / 2;

            // Draw Shadow
            SetBlendMode(BlendModes.Average);

            numSegments = 20;
            GL.Begin(PrimitiveType.TriangleFan);
                GL.Color4f(0.1f, 0.1f, 0.1f, 0.5f);
                GL.Vertex3f(cx, 0, cz);

                for (si = 0; si <= numSegments; si++)
                {
                    GL.Color4f(0.1f, 0.1f, 0.1f, 0);
                    GL.Vertex3f((float)(ground_radius * Math.Sin(si * 2 * Math.PI / numSegments) + cx), 0,
                               (float)(ground_radius * Math.Cos(si * 2 * Math.PI / numSegments) + cz));
                }
            GL.End();

            GL.Enable(EnableCap.DepthTest);
            GL.Disable(EnableCap.Fog);

            // Draw underlying box (just depth)
            GL.ColorMask(false, false, false, false);
            GL.Color3f(1, 1, 1);
            GL.Begin(PrimitiveType.Quads);
                GL.Vertex3f(p_max.X, 0, p_max.Z);
                GL.Vertex3f(p_max.X, 0, p_min.Z);
                GL.Vertex3f(p_min.X, 0, p_min.Z);
                GL.Vertex3f(p_min.X, 0, p_max.Z);

                GL.Vertex3f(p_max.X, 0, p_max.Z);
                GL.Vertex3f(p_max.X, sub_y, p_max.Z);
                GL.Vertex3f(p_max.X, sub_y, p_min.Z);
                GL.Vertex3f(p_max.X, 0, p_min.Z);

                GL.Vertex3f(p_max.X, 0, p_min.Z);
                GL.Vertex3f(p_max.X, sub_y, p_min.Z);
                GL.Vertex3f(p_min.X, sub_y, p_min.Z);
                GL.Vertex3f(p_min.X, 0, p_min.Z);

                GL.Vertex3f(p_min.X, sub_y, p_max.Z);
                GL.Vertex3f(p_min.X, 0, p_max.Z);
                GL.Vertex3f(p_min.X, 0, p_min.Z);
                GL.Vertex3f(p_min.X, sub_y, p_min.Z);

                GL.Vertex3f(p_max.X, sub_y, p_max.Z);
                GL.Vertex3f(p_max.X, 0, p_max.Z);
                GL.Vertex3f(p_min.X, 0, p_max.Z);
                GL.Vertex3f(p_min.X, sub_y, p_max.Z);
            GL.End();
            GL.ColorMask(true, true, true, true);
        }

        public static void DrawPlane(ref double[] planeTransformation, ref Vector3 planeOriginalPoint1,
                                                                       ref Vector3 planeOriginalPoint2,
                                                                       ref Vector3 planeOriginalPoint3,
                                                                       ref Vector3 planeOriginalPoint4)
        {
            Vector3 p1 = new Vector3();
            Vector3 p2 = new Vector3();
            Vector3 p3 = new Vector3();
            Vector3 p4 = new Vector3();

            GL.Disable(EnableCap.CullFace);

            SetBlendMode(BlendModes.Average);

            MultiplyPoint3DByOGLMatrix(planeTransformation, planeOriginalPoint1, ref p1);
            MultiplyPoint3DByOGLMatrix(planeTransformation, planeOriginalPoint2, ref p2);
            MultiplyPoint3DByOGLMatrix(planeTransformation, planeOriginalPoint3, ref p3);
            MultiplyPoint3DByOGLMatrix(planeTransformation, planeOriginalPoint4, ref p4);

            GL.PolygonMode(TriangleFace.Front, PolygonMode.Fill);
            GL.PolygonMode(TriangleFace.Back, PolygonMode.Line);

            GL.Color4f(1, 0, 0, 0.10f);

            GL.Begin(PrimitiveType.Quads);
                GL.Vertex3f(p1.X, p1.Y, p1.Z);
                GL.Vertex3f(p2.X, p2.Y, p2.Z);
                GL.Vertex3f(p3.X, p3.Y, p3.Z);
                GL.Vertex3f(p4.X, p4.Y, p4.Z);
            GL.End();
        }

        public static void DrawAxes(Control pbIn, int iFrame, int ianimIndex)
        {
            float letterWidth, letterHeight;
            float max_x, max_y, max_z;

            Vector3 pX = new Vector3();
            Vector3 pY = new Vector3();
            Vector3 pZ = new Vector3();

            Vector3 p_max = new Vector3();
            Vector3 p_min = new Vector3();

            GL.Disable(EnableCap.Lighting);

            switch(modelType)
            {
                case K_HRC_SKELETON:
                    ComputeFieldBoundingBox(fSkeleton, fAnimation.frames[iFrame],
                                            ref p_min, ref p_max);
                    break;

                case K_AA_SKELETON:
                case K_MAGIC_SKELETON:
                    ComputeBattleBoundingBox(bSkeleton, bAnimationsPack.SkeletonAnimations[ianimIndex].frames[iFrame],
                                             ref p_min, ref p_max);

                    break;

                case K_3DS_MODEL:
                case K_P_BATTLE_MODEL:
                case K_P_FIELD_MODEL:
                case K_P_MAGIC_MODEL:
                    ComputePModelBoundingBox(fPModel, ref p_min, ref p_max);

                    break;

            }

            //max_x = Math.Abs(p_min.X) > Math.Abs(p_max.X) ? p_min.X : p_max.X;
            //max_y = Math.Abs(p_min.Y) > Math.Abs(p_max.Y) ? p_min.Y : p_max.Y;
            //max_z = Math.Abs(p_min.Z) > Math.Abs(p_max.Z) ? p_min.Z : p_max.Z;

            max_x = Math.Abs(p_max.X - p_min.X);
            max_y = Math.Abs(p_max.Y - p_min.Y);
            max_z = Math.Abs(p_max.Z - p_min.Z);

            if (max_x > max_y && max_x > max_z) max_y = max_z = max_x;
            if (max_y > max_x && max_y > max_z) max_x = max_z = max_y;
            if (max_z > max_x && max_z > max_y) max_x = max_y = max_z;

            GL.Begin(PrimitiveType.Lines);
                GL.Color3f(1, 0, 0);
                GL.Vertex3f(0, 0, 0);
                GL.Vertex3f(max_x, 0, 0);

                if (bSkeleton.IsBattleLocation)
                {
                    GL.Color3f(0, 1, 0);
                    GL.Vertex3f(0, 0, 0);
                    GL.Vertex3f(0, -max_y, 0);

                    GL.Color3f(0, 0, 1);
                    GL.Vertex3f(0, 0, 0);
                    GL.Vertex3f(0, 0, max_z);
                }
                else 
                {
                    GL.Color3f(0, 0, 1);
                    GL.Vertex3f(0, 0, 0);
                    GL.Vertex3f(0, -max_y, 0);

                    GL.Color3f(0, 1, 0);
                    GL.Vertex3f(0, 0, 0);
                    GL.Vertex3f(0, 0, max_z);
                }
            GL.End();

            //  Get projected end of the X axis
            pX.X = max_x;
            pX.Y = 0;
            pX.Z = 0;
            pX = GetProjectedCoords(pX);

            //  Get projected end of the Y axis
            pY.X = 0;
            pY.Y = -max_y;
            pY.Z = 0;
            pY = GetProjectedCoords(pY);

            //  Get projected end of the Z axis
            pZ.X = 0;
            pZ.Y = 0;
            pZ.Z = max_z;
            pZ = GetProjectedCoords(pZ);


            //  Set 2D mode to draw letters
            GL.MatrixMode(MatrixMode.Projection);
            GL.LoadIdentity();
            gluOrtho2D(0, pbIn.ClientRectangle.Width, 0, pbIn.ClientRectangle.Height);
            GL.MatrixMode(MatrixMode.Modelview);
            GL.LoadIdentity();

            letterWidth = LETTER_SIZE;
            letterHeight = (float)(LETTER_SIZE * 1.5);
            GL.Disable(EnableCap.DepthTest);

            GL.Begin(PrimitiveType.Lines);
                //  Draw X
                GL.Color3f(0, 0, 0);
                GL.Vertex2f(pX.X - letterWidth, pX.Y - letterHeight);
                GL.Vertex2f(pX.X + letterWidth, pX.Y + letterHeight);
                GL.Vertex2f(pX.X - letterWidth, pX.Y + letterHeight);
                GL.Vertex2f(pX.X + letterWidth, pX.Y - letterHeight);

                if (bSkeleton.IsBattleLocation)
                {
                    //  Draw Y
                    GL.Color3f(0, 0, 0);
                    GL.Vertex2f(pY.X - letterWidth, pY.Y - letterHeight);
                    GL.Vertex2f(pY.X + letterWidth, pY.Y + letterHeight);
                    GL.Vertex2f(pY.X - letterWidth, pY.Y + letterHeight);
                    GL.Vertex2f(pY.X, pY.Y);

                    //  Draw Z
                    GL.Color3f(0, 0, 0);
                    GL.Vertex2f(pZ.X + letterWidth, pZ.Y + letterHeight);
                    GL.Vertex2f(pZ.X - letterWidth, pZ.Y + letterHeight);

                    GL.Vertex2f(pZ.X + letterWidth, pZ.Y + letterHeight);
                    GL.Vertex2f(pZ.X - letterWidth, pZ.Y - letterHeight);

                    GL.Vertex2f(pZ.X + letterWidth, pZ.Y - letterHeight);
                    GL.Vertex2f(pZ.X - letterWidth, pZ.Y - letterHeight);
                }
                else 
                {
                    //  Draw Y
                    GL.Color3f(0, 0, 0);
                    GL.Vertex2f(pZ.X - letterWidth, pZ.Y - letterHeight);
                    GL.Vertex2f(pZ.X + letterWidth, pZ.Y + letterHeight);
                    GL.Vertex2f(pZ.X - letterWidth, pZ.Y + letterHeight);
                    GL.Vertex2f(pZ.X, pZ.Y);

                    //  Draw Z
                    GL.Color3f(0, 0, 0);
                    GL.Vertex2f(pY.X + letterWidth, pY.Y + letterHeight);
                    GL.Vertex2f(pY.X - letterWidth, pY.Y + letterHeight);

                    GL.Vertex2f(pY.X + letterWidth, pY.Y + letterHeight);
                    GL.Vertex2f(pY.X - letterWidth, pY.Y - letterHeight);

                    GL.Vertex2f(pY.X + letterWidth, pY.Y - letterHeight);
                    GL.Vertex2f(pY.X - letterWidth, pY.Y - letterHeight);
                }
            GL.End();

            GL.Enable(EnableCap.DepthTest);
        }

        public static void KillPrecalculatedLighting(PModel Model, ref PairIB[] translationTableVertex)
        {
            int iVColorIdx;

            for (iVColorIdx = 0; iVColorIdx < Model.Header.numVerts; iVColorIdx++)
            {
                translationTableVertex[iVColorIdx].B = 1;
            }
        }

        /// <summary>
        /// Draws a P-model in the editor using the provided context.
        /// This overload is decoupled from FrmPEditor and can be used in other applications.
        /// </summary>
        public static void DrawPModelEditor(PEditorContext ctx, Control pbIn)
        {
            Vector3 p_min = new Vector3();
            Vector3 p_max = new Vector3();

            GL.Viewport(0, 0, pbIn.ClientRectangle.Width,
                             pbIn.ClientRectangle.Height);
            ClearPanel();

            SetCameraPModel(ctx.EditedModel,
                            ctx.Camera.PanX, ctx.Camera.PanY, ctx.Camera.PanZ + ctx.Camera.Distance,
                            ctx.Camera.Alpha, ctx.Camera.Beta, ctx.Camera.Gamma,
                            1, 1, 1);

            ConcatenateCameraModelView(ctx.Transform.RepositionX, ctx.Transform.RepositionY, ctx.Transform.RepositionZ,
                                       0, 0, 0,
                                       ctx.Transform.ResizeX, ctx.Transform.ResizeY, ctx.Transform.ResizeZ);

            GL.MatrixMode(MatrixMode.Modelview);
            GL.PushMatrix();

            ComputePModelBoundingBox(ctx.EditedModel, ref p_min, ref p_max);
            float modelDiameterNormalized = (-2 * ComputeSceneRadius(p_min, p_max)) / LIGHT_STEPS_PEDITOR;

            if (ctx.EnableLighting)
            {
                GL.Disable(EnableCap.Light0);
                GL.Disable(EnableCap.Light1);
                GL.Disable(EnableCap.Light2);
                GL.Disable(EnableCap.Light3);

                SetLighting(LightName.Light0, modelDiameterNormalized * ctx.LightX,
                                            modelDiameterNormalized * ctx.LightY,
                                            modelDiameterNormalized * ctx.LightZ,
                                            1, 1, 1, false);
            }
            else GL.Disable(EnableCap.Lighting);

            // Sync lighting settings to modern renderer
            GLRenderer.LightingEnabled = ctx.EnableLighting;
            if (ctx.EnableLighting)
            {
                GLRenderer.LightPosition = new Vector3(
                    modelDiameterNormalized * ctx.LightX,
                    modelDiameterNormalized * ctx.LightY,
                    modelDiameterNormalized * ctx.LightZ);
            }

            SetDefaultOGLRenderState();

            // Local model reference for ref parameter
            PModel model = ctx.EditedModel;
            uint[] texIds = ctx.TextureIds ?? tex_ids;

            switch (ctx.DrawMode)
            {
                case DrawMode.K_MESH:
                    DrawPModelWireframe(ref model, true);
                    break;

                case DrawMode.K_PCOLORS:
                    GL.Enable(EnableCap.PolygonOffsetFill);
                    GL.PolygonOffset(1, 1);
                    DrawPModelPolygonColors(ref model, true);
                    GL.Disable(EnableCap.PolygonOffsetFill);

                    DrawPModelWireframe(ref model, true);
                    break;

                case DrawMode.K_VCOLORS:
                    DrawPModel(ref model, ref texIds, true);
                    break;
            }

            GL.PopMatrix();
        }

        public static void DrawAxesPE(Control pbIn, PModel editedPModel)
        {
            float letterWidth, letterHeight;
            float max_x, max_y, max_z;

            Vector3 pX = new Vector3();
            Vector3 pY = new Vector3();
            Vector3 pZ = new Vector3();

            Vector3 p_max = new Vector3();
            Vector3 p_min = new Vector3();

            GL.Disable(EnableCap.Lighting);
            ComputePModelBoundingBox(editedPModel, ref p_min, ref p_max);

            max_x = Math.Abs(p_min.X) > Math.Abs(p_max.X) ? p_min.X : p_max.X;
            max_y = Math.Abs(p_min.Y) > Math.Abs(p_max.Y) ? p_min.Y : p_max.Y;
            max_z = Math.Abs(p_min.Z) > Math.Abs(p_max.Z) ? p_min.Z : p_max.Z;

            GL.Begin(PrimitiveType.Lines);
                GL.Color3f(1, 0, 0);
                GL.Vertex3f(0, 0, 0);
                GL.Vertex3f(2 * max_x, 0, 0);

                GL.Color3f(0, 1, 0);
                GL.Vertex3f(0, 0, 0);
                GL.Vertex3f(0, 2 * max_y, 0);

                GL.Color3f(0, 0, 1);
                GL.Vertex3f(0, 0, 0);
                GL.Vertex3f(0, 0, 2 * max_z);
            GL.End();

            //  Get projected end of the X axis
            pX.X = 2 * max_x;
            pX.Y = 0;
            pX.Z = 0;
            pX = GetProjectedCoords(pX);

            //  Get projected end of the Y axis
            pY.X = 0;
            pY.Y = 2 * max_y;
            pY.Z = 0;
            pY = GetProjectedCoords(pY);

            //  Get projected end of the Z axis
            pZ.X = 0;
            pZ.Y = 0;
            pZ.Z = 2 * max_z;
            pZ = GetProjectedCoords(pZ);


            //  Set 2D mode to draw letters
            GL.MatrixMode(MatrixMode.Projection);
            GL.LoadIdentity();
            gluOrtho2D(0, pbIn.ClientRectangle.Width, 0, pbIn.ClientRectangle.Height);
            GL.MatrixMode(MatrixMode.Modelview);
            GL.LoadIdentity();

            letterWidth = LETTER_SIZE;
            letterHeight = (float)(LETTER_SIZE * 1.5);
            GL.Disable(EnableCap.DepthTest);

            GL.Begin(PrimitiveType.Lines);
                //  Draw X
                GL.Color3f(0, 0, 0);
                GL.Vertex2f(pX.X - letterWidth, pX.Y - letterHeight);
                GL.Vertex2f(pX.X + letterWidth, pX.Y + letterHeight);
                GL.Vertex2f(pX.X - letterWidth, pX.Y + letterHeight);
                GL.Vertex2f(pX.X + letterWidth, pX.Y - letterHeight);

                //  Draw Y
                GL.Color3f(0, 0, 0);
                GL.Vertex2f(pY.X - letterWidth, pY.Y - letterHeight);
                GL.Vertex2f(pY.X + letterWidth, pY.Y + letterHeight);
                GL.Vertex2f(pY.X - letterWidth, pY.Y + letterHeight);
                GL.Vertex2f(pY.X, pY.Y);

                //  Draw Z
                GL.Color3f(0, 0, 0);
                GL.Vertex2f(pZ.X + letterWidth, pZ.Y + letterHeight);
                GL.Vertex2f(pZ.X - letterWidth, pZ.Y + letterHeight);

                GL.Vertex2f(pZ.X + letterWidth, pZ.Y + letterHeight);
                GL.Vertex2f(pZ.X - letterWidth, pZ.Y - letterHeight);

                GL.Vertex2f(pZ.X + letterWidth, pZ.Y - letterHeight);
                GL.Vertex2f(pZ.X - letterWidth, pZ.Y - letterHeight);
            GL.End();

            GL.Enable(EnableCap.DepthTest);
        }

        public static int AddPaintVertex(ref PModel Model, int iGroupIdx, Vector3 vPoint3D, Color vColor)
        {
            int iAddVertexResult = -1;
            //  -------- Warning! Causes the Normals to be inconsistent if lights are disabled.------------------
            //  --------------------------------Must call ComputeNormals ----------------------------------------
            int iVertIdx, iNormalIdx, iTexCoordIdx, baseVerts, baseNormals, baseTexCoords;
            int iNextGroup;

            Model.Header.numVerts++;
            Array.Resize(ref Model.Verts, Model.Header.numVerts);
            Array.Resize(ref Model.Vcolors, Model.Header.numVerts);

            if (Model.Groups[iGroupIdx].texFlag == 1)
            {
                Model.Header.numTexCs++;
                Array.Resize(ref Model.TexCoords, Model.Header.numTexCs);
            }

            Model.Header.numNormals = Model.Header.numVerts;
            Array.Resize(ref Model.Normals, Model.Header.numNormals);
            Model.Header.numNormIdx = Model.Header.numVerts;
            Array.Resize(ref Model.NormalIndex, Model.Header.numNormIdx);
            Model.NormalIndex[Model.Header.numNormIdx - 1] = Model.Header.numNormIdx - 1;

            iNextGroup = GetNextGroup(Model, iGroupIdx);
            //if (iGroupIdx < Model.Header.numGroups - 1)
            if (iNextGroup != -1)
            {
                //baseVerts = Model.Groups[iGroupIdx + 1].offsetVert;
                baseVerts = Model.Groups[iNextGroup].offsetVert;

                for (iVertIdx = Model.Header.numVerts - 1; iVertIdx >= baseVerts; iVertIdx--)
                {
                    Model.Verts[iVertIdx] = Model.Verts[iVertIdx - 1];
                    Model.Vcolors[iVertIdx] = Model.Vcolors[iVertIdx - 1];
                }

                //if (GL.IsEnabled(EnableCap.Lighting))
                //{
                if (Model.Groups[iGroupIdx].texFlag == 1)
                    {
                        //baseNormals = Model.Groups[iGroupIdx + 1].offsetVert;
                        baseNormals = Model.Groups[iNextGroup].offsetVert;

                        for (iNormalIdx = Model.Header.numNormals - 1; iNormalIdx >= baseNormals; iNormalIdx--)
                        {
                            Model.Normals[iNormalIdx] = Model.Normals[iNormalIdx - 1];
                        }
                    }
                //}

                if (Model.Groups[iGroupIdx].texFlag == 1)
                {
                    baseTexCoords = Model.Groups[iGroupIdx].offsetTex + Model.Groups[iGroupIdx].numVert;

                    for (iTexCoordIdx = Model.Header.numTexCs - 1; iTexCoordIdx >= baseTexCoords; iTexCoordIdx--)
                    {
                        Model.TexCoords[iTexCoordIdx] = Model.TexCoords[iTexCoordIdx - 1];
                    }
                }

                while (iNextGroup != -1)
                {
                    Model.Groups[iNextGroup].offsetVert++;

                    if ((Model.Groups[iGroupIdx].texFlag == 1) && (Model.Groups[iNextGroup].texFlag == 1))
                    {
                        Model.Groups[iNextGroup].offsetTex++;
                    }

                    iNextGroup = GetNextGroup(Model, iNextGroup);
                }

            }

            if (iGroupIdx < Model.Header.numGroups)
            {
                Model.Verts[Model.Groups[iGroupIdx].offsetVert + Model.Groups[iGroupIdx].numVert] = vPoint3D;
                Model.Vcolors[Model.Groups[iGroupIdx].offsetVert + Model.Groups[iGroupIdx].numVert] = vColor;
                iAddVertexResult = Model.Groups[iGroupIdx].offsetVert + Model.Groups[iGroupIdx].numVert;
                Model.Groups[iGroupIdx].numVert++;
            }

            return iAddVertexResult;
        }

        public static int PaintVertex(ref PModel Model, int iGroupIdx, int iPolyIdx, int iVertIdx,
                                      byte bR, byte bG, byte bB, bool bTextured)
        {

            int iPaintVertexResult = -1;

            Vector3 tmpVert;
            Color tmpColor;

            if (Model.Vcolors[Model.Polys[iPolyIdx].Verts[iVertIdx] +
                              Model.Groups[iGroupIdx].offsetVert].R == bR &&
                Model.Vcolors[Model.Polys[iPolyIdx].Verts[iVertIdx] +
                              Model.Groups[iGroupIdx].offsetVert].G == bG &&
                Model.Vcolors[Model.Polys[iPolyIdx].Verts[iVertIdx] +
                              Model.Groups[iGroupIdx].offsetVert].B == bB)
                        iPaintVertexResult = Model.Polys[iPolyIdx].Verts[iVertIdx];

            if (iPaintVertexResult == -1)
            {

                tmpVert.X = Model.Verts[Model.Polys[iPolyIdx].Verts[iVertIdx] +
                                        Model.Groups[iGroupIdx].offsetVert].X;
                tmpVert.Y = Model.Verts[Model.Polys[iPolyIdx].Verts[iVertIdx] +
                                        Model.Groups[iGroupIdx].offsetVert].Y;
                tmpVert.Z = Model.Verts[Model.Polys[iPolyIdx].Verts[iVertIdx] +
                                        Model.Groups[iGroupIdx].offsetVert].Z;

                tmpColor = Color.FromArgb(255, bR, bG, bB);

                iPaintVertexResult = AddPaintVertex(ref Model, iGroupIdx, tmpVert, tmpColor) - Model.Groups[iGroupIdx].offsetVert;

                Model.Normals[Model.NormalIndex[iPaintVertexResult]] =
                    new Vector3(Model.Normals[Model.NormalIndex[Model.Polys[iPolyIdx].Verts[iVertIdx]]+
                                                                Model.Groups[iGroupIdx].offsetVert].X,
                                Model.Normals[Model.NormalIndex[Model.Polys[iPolyIdx].Verts[iVertIdx]] +
                                                                Model.Groups[iGroupIdx].offsetVert].Y,
                                Model.Normals[Model.NormalIndex[Model.Polys[iPolyIdx].Verts[iVertIdx]] +
                                                                Model.Groups[iGroupIdx].offsetVert].Z);

                if (bTextured)
                {
                    Model.TexCoords[Model.Groups[iGroupIdx].offsetTex + iPaintVertexResult] =
                        new Vector2(Model.TexCoords[Model.Polys[iPolyIdx].Verts[iVertIdx] + Model.Groups[iGroupIdx].offsetTex].X,
                                    Model.TexCoords[Model.Polys[iPolyIdx].Verts[iVertIdx] + Model.Groups[iGroupIdx].offsetTex].Y);
                }

                //    // -- Commented in KimeraVB6
                //    // Model.Normals[iPaintVertexResult + Model.Groups[iGroupIdx].offsetVert] = Model.Normals[iVertIdx];
                //}
                //else
                //{
                //    // Debug.Print "Substituido por: " + Str$(PaintVertex);
                //}
            }
            
            return iPaintVertexResult;
        }

        //  ------------------------------WARNINGS!----------------------------------
        //  -------*Can causes the Normals to be inconsistent (call ComputeNormals).--
        //  -------*Can causes inconsistent edges (call ComputeEdges).----------------
        //  -------*Can cause unused vertices (call KillUnusedVertices).--------------
        public static void PaintPolygon(ref PModel Model, int iPolyIdx, byte bR, byte bG, byte bB)
        {
            int iGroupIdx, iVertIdx;

            iGroupIdx = GetPolygonGroup(Model, iPolyIdx);

            for (iVertIdx = 0; iVertIdx < 3; iVertIdx++)
            {
                Model.Polys[iPolyIdx].Verts[iVertIdx] = 
                    (ushort)PaintVertex(ref Model, iGroupIdx, iPolyIdx, iVertIdx, bR, bG, bB,
                                        Model.Groups[iGroupIdx].texFlag == 1);
                //  'Debug.Print "Vert(:", .Verts(vi), ",", Group, ")", obj.Verts(.Verts(vi) + obj.Groups(Group).offVert).X, obj.Verts(.Verts(vi) + obj.Groups(Group).offVert).Y, obj.Verts(.Verts(vi) + obj.Groups(Group).offVert).Z
            }

            Model.Pcolors[iPolyIdx] = Color.FromArgb(Model.Pcolors[iPolyIdx].A, bR, bG, bB);
        }

        public static Color ComputePolyColor(PModel Model, int iPolyIdx)
        {
            int iGroupIdx, iVertIdx;
            int iTmpA = 0, iTmpR = 0, iTmpG = 0, iTmpB = 0;

            iGroupIdx = GetPolygonGroup(Model, iPolyIdx);

            for (iVertIdx = 0; iVertIdx < 3; iVertIdx++)
            {

                iTmpA += Model.Vcolors[Model.Polys[iPolyIdx].Verts[iVertIdx] +
                                       Model.Groups[iGroupIdx].offsetVert].A;
                iTmpR += Model.Vcolors[Model.Polys[iPolyIdx].Verts[iVertIdx] +
                                       Model.Groups[iGroupIdx].offsetVert].R;
                iTmpG += Model.Vcolors[Model.Polys[iPolyIdx].Verts[iVertIdx] +
                                       Model.Groups[iGroupIdx].offsetVert].G;
                iTmpB += Model.Vcolors[Model.Polys[iPolyIdx].Verts[iVertIdx] +
                                       Model.Groups[iGroupIdx].offsetVert].B;
            }

            return Color.FromArgb(iTmpA / 3, iTmpR / 3, iTmpG / 3, iTmpB / 3);
        }




    }
}

