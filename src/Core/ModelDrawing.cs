using System;
using System.Collections.Generic;
using System.Drawing;
using System.Windows.Forms;
using OpenTK.Graphics.OpenGL.Compatibility;
using OpenTK.Mathematics;

namespace KimeraCS
{
    using Rendering;
    using static Rendering.VisualizationHelpers;

    using static FrmSkeletonEditor;
    using static FrmPEditor;

    using static FF7Skeleton;
    using static FF7FieldSkeleton;
    using static FF7FieldAnimation;
    using static FF7FieldRSDResource;

    using static FF7BattleSkeleton;
    using static FF7BattleAnimation;

    using static FF7PModel;

    using static Lighting;

    using static Utils;
    using KimeraCS.Core;

    public static class ModelDrawing
    {
        public static uint[] tex_ids = new uint[1];


        //  ---------------------------------------------------------------------------------------------------
        //  ================================== GENERIC FIELD/BATTLE DRAW  =====================================
        //  ---------------------------------------------------------------------------------------------------
        // Cache for normals mesh
        private static LineMesh _normalsMesh;

        public static void ShowNormals(PGroup Group, PPolygon[] Polys, Point3D[] Verts,
                                       Point3D[] Normals, int[] NormalsIndex)
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
                                             bShowVertexNormals, bShowFaceNormals,
                                             fNormalsScale, iNormalsColor);

            if (_normalsMesh != null)
            {
                GLRenderer.DrawLinesModern(_normalsMesh);
            }

            // Restore original matrices
            GLRenderer.ProjectionMatrix = savedProjection;
            GLRenderer.ViewMatrix = savedView;
            GLRenderer.ModelMatrix = savedModel;
        }

        public static void DrawGroup(PGroup Group, PPolygon[] Polys, Point3D[] Verts,
                                     Color[] Vcolors, Point3D[] Normals, int[] NormalsIndex,
                                     Point2D[] TexCoords, PHundret Hundret, bool HideHiddenQ)
        {

            if (Group.HiddenQ && HideHiddenQ) return;

            int iPolyIdx = 0, iVertIdx = 0;
            float x, y, z;

            try
            {
                GL.Begin(PrimitiveType.Triangles);
                GL.ColorMaterial(TriangleFace.FrontAndBack, ColorMaterialParameter.AmbientAndDiffuse);

                for (iPolyIdx = Group.offsetPoly; iPolyIdx < Group.offsetPoly + Group.numPoly; iPolyIdx++)
                {
                    for (iVertIdx = 0; iVertIdx < 3; iVertIdx++)
                    {
                        if (Hundret.blend_mode == 0 && !(Hundret.shademode == 1) && !bSkeleton.IsBattleLocation)
                            GL.Color4f(Vcolors[Polys[iPolyIdx].Verts[iVertIdx] + Group.offsetVert].R / 255.0f,
                                      Vcolors[Polys[iPolyIdx].Verts[iVertIdx] + Group.offsetVert].G / 255.0f,
                                      Vcolors[Polys[iPolyIdx].Verts[iVertIdx] + Group.offsetVert].B / 255.0f,
                                      0.5f);
                        else
                            GL.Color4f(Vcolors[Polys[iPolyIdx].Verts[iVertIdx] + Group.offsetVert].R / 255.0f,
                                      Vcolors[Polys[iPolyIdx].Verts[iVertIdx] + Group.offsetVert].G / 255.0f,
                                      Vcolors[Polys[iPolyIdx].Verts[iVertIdx] + Group.offsetVert].B / 255.0f,
                                      1.0f);

                        if (Normals.Length > 0)
                            GL.Normal3f(Normals[NormalsIndex[Polys[iPolyIdx].Verts[iVertIdx] + Group.offsetVert]].x,
                                       Normals[NormalsIndex[Polys[iPolyIdx].Verts[iVertIdx] + Group.offsetVert]].y,
                                       Normals[NormalsIndex[Polys[iPolyIdx].Verts[iVertIdx] + Group.offsetVert]].z);

                        if (TexCoords != null)
                            if (Group.texFlag == 1 && TexCoords.Length > 0)
                            {
                                x = TexCoords[Group.offsetTex + Polys[iPolyIdx].Verts[iVertIdx]].x;
                                y = TexCoords[Group.offsetTex + Polys[iPolyIdx].Verts[iVertIdx]].y;
                                GL.TexCoord2f(x, y);
                            }

                        x = Verts[Polys[iPolyIdx].Verts[iVertIdx] + Group.offsetVert].x;
                        y = Verts[Polys[iPolyIdx].Verts[iVertIdx] + Group.offsetVert].y;
                        z = Verts[Polys[iPolyIdx].Verts[iVertIdx] + Group.offsetVert].z;
                        GL.Vertex3f(x, y, z);
                    }
                }

                GL.End();

                // Let's try here to render the normals
                if (bShowVertexNormals || bShowFaceNormals)
                    ShowNormals(Group, Polys, Verts, Normals, NormalsIndex);

            }

            catch (Exception ex)
            {
                throw new KimeraException(
                    "Exception in DrawGroup procedure. Group: " + Group.realGID +
                    ", Polygon (iPolyIdx): " + iPolyIdx +
                    ", Vertex (iVertIdx): " + iVertIdx +
                    ", offsetVertex: " + Group.offsetVert +
                    ", offsetPolygon: " + Group.offsetPoly, ex);
            }
        }

        public static bool IsColorKey(PModel Model, int iGroupIdx)
        {
            bool bIsColorKey = false;

            switch (modelType)
            {
                case K_HRC_SKELETON:
                    if (fSkeleton.textures_pool != null)
                        if (fSkeleton.textures_pool.Count > 0 && Model.Groups[iGroupIdx].texID >= 0)
                            if (fSkeleton.textures_pool[Model.Groups[iGroupIdx].texID].ColorKeyFlag == 1)
                                bIsColorKey = true;

                    break;

                case K_AA_SKELETON:
                case K_MAGIC_SKELETON:
                    if (bSkeleton.textures.Count > 0 && Model.Groups[iGroupIdx].texID >= 0)
                        if (bSkeleton.textures[Model.Groups[iGroupIdx].texID].ColorKeyFlag == 1)
                            bIsColorKey = true;

                    break;
            }

            return bIsColorKey;
        }

        public static void DrawPModel(ref PModel Model, ref uint[] tex_ids, bool HideHiddenGroupsQ)
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
            if (bShowVertexNormals || bShowFaceNormals)
            {
                for (int g = 0; g < Model.Header.numGroups; g++)
                {
                    ShowNormals(Model.Groups[g], Model.Polys, Model.Verts,
                                Model.Normals, Model.NormalIndex);
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
            GLRenderer.DrawPModelWireframe(ref Model, new OpenTK.Mathematics.Vector3(0, 0, 0), HideHiddenGroupsQ);

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

        public static void DrawGroupDList(ref PGroup Group)
        {
            GL.CallList(Group.DListNum);
        }

        public static void DrawPModelDLists(ref PModel Model, ref uint[] tex_ids)
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
            GLRenderer.ViewMatrix = OpenTK.Mathematics.Matrix4.Identity;
            GLRenderer.ModelMatrix = legacyModelView;

            GLRenderer.DrawPModelModern(ref Model, tex_ids, false);

            // Show normals if enabled (must be done after model drawing)
            if (bShowVertexNormals || bShowFaceNormals)
            {
                for (int g = 0; g < Model.Header.numGroups; g++)
                {
                    ShowNormals(Model.Groups[g], Model.Polys, Model.Verts,
                                Model.Normals, Model.NormalIndex);
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
        private static LineMesh _boundingBoxMesh;

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
                max_x = -(float)INFINITY_SINGLE;
                max_y = -(float)INFINITY_SINGLE;
                max_z = -(float)INFINITY_SINGLE;

                min_x = (float)INFINITY_SINGLE;
                min_y = (float)INFINITY_SINGLE;
                min_z = (float)INFINITY_SINGLE;

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

        public static void DrawRSDResource(FieldRSDResource fRSDResource, bool bDListsEnable)
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

            if (!bDListsEnable) DrawPModel(ref fRSDResource.Model, ref tex_ids, false);
            else DrawPModelDLists(ref fRSDResource.Model, ref tex_ids);

            GL.PopMatrix();
        }

        public static void DrawFieldBone(FieldBone bone, bool bDListsEnable)
        {

            int iResourceIdx;

            GL.MatrixMode(MatrixMode.Modelview);
            GL.PushMatrix();

            GL.Scaled(bone.resizeX, bone.resizeY, bone.resizeZ);

            for (iResourceIdx = 0; iResourceIdx < bone.nResources; iResourceIdx++)
                DrawRSDResource(bone.fRSDResources[iResourceIdx], bDListsEnable);

            GL.PopMatrix();
        }

        public static void DrawFieldSkeleton(FieldSkeleton fSkeleton, FieldFrame fFrame, bool bDListsEnable)
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

                DrawFieldBone(fSkeleton.bones[iBoneIdx], bDListsEnable);

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

        public static void DrawBattleWeaponBoundingBox(BattleSkeleton bSkeleton, BattleFrame wpFrame, int weaponIndex)
        {
            double[] rot_mat = new double[16];

            //if (weaponIndex > -1 && bSkeleton.nWeapons > 0)       // -- Commented in KimeraVB6
            if (ianimWeaponIndex > -1 && bSkeleton.wpModels.Count > 0 && bAnimationsPack.WeaponAnimations.Count > 0)
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

        public static void DrawBattleSkeletonBone(BattleBone bBone, uint[] texIDS, bool bDListsEnable)
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
                        DrawPModel(ref tmpbModel, ref texIDS, false);
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
                        DrawPModelDLists(ref tmpbModel, ref texIDS);
                        bBone.Models[iModelIdx] = tmpbModel;

                        GL.PopMatrix();
                    }
                }
            }

            GL.PopMatrix();
        }

        public static void DrawBattleSkeleton(BattleSkeleton bSkeleton, BattleFrame bFrame, BattleFrame wpFrame,
                                              int weaponIndex, bool bDListsEnable)
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
                    DrawBattleSkeletonBone(bSkeleton.bones[iBoneIdx], bSkeleton.TexIDS, false);
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

                    DrawBattleSkeletonBone(bSkeleton.bones[iBoneIdx], bSkeleton.TexIDS, bDListsEnable);

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
            if (ianimWeaponIndex > -1 && bSkeleton.wpModels.Count > 0 && bAnimationsPack.WeaponAnimations.Count > 0)
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
                    DrawPModelDLists(ref tmpwpModel, ref bSkeleton.TexIDS);
                    bSkeleton.wpModels[weaponIndex] = tmpwpModel;
                }
                else
                {
                    tmpwpModel = bSkeleton.wpModels[weaponIndex];
                    DrawPModel(ref tmpwpModel, ref bSkeleton.TexIDS, false);
                    bSkeleton.wpModels[weaponIndex] = tmpwpModel;
                }
                GL.PopMatrix();

                GL.PopMatrix();
            }
        }



        //  ---------------------------------------------------------------------------------------------------
        //  ======================================= SKELETON DRAW  ============================================
        //  ---------------------------------------------------------------------------------------------------
        public static void DrawSkeletonModel(bool bDListsEnable)
        {

            double[] rot_mat = new double[16];
            BattleFrame tmpbFrame;

            Point3D p_min = new Point3D();
            Point3D p_max = new Point3D();

            try
            {
                switch (modelType)
                {
                    case K_3DS_MODEL:
                    case K_P_FIELD_MODEL:
                    case K_P_BATTLE_MODEL:
                    case K_P_MAGIC_MODEL:
                        ComputePModelBoundingBox(fPModel, ref p_min, ref p_max);

                        SetCameraAroundModel(ref p_min, ref p_max,
                                             panX, panY, (float)(panZ + DIST),
                                             (float)alpha, (float)beta, (float)gamma, 1, 1, 1);

                        if (bShowGround)
                        {
                            GL.Disable(EnableCap.Lighting);
                            DrawGround();
                            DrawShadow(ref p_min, ref p_max);
                        }

                        SetLights();

                        GL.MatrixMode(MatrixMode.Modelview);
                        GL.PushMatrix();

                        GL.Translated(fPModel.repositionX,
                                     fPModel.repositionY,
                                     fPModel.repositionZ);

                        BuildRotationMatrixWithQuaternionsXYZ(fPModel.rotateAlpha,
                                                              fPModel.rotateBeta,
                                                              fPModel.rotateGamma,
                                                              ref rot_mat);

                        GL.MultMatrixd(rot_mat);
                        GL.Scaled(fPModel.resizeX,
                                 fPModel.resizeY,
                                 fPModel.resizeZ);

                        DrawPModel(ref fPModel, ref tex_ids, false);

                        GL.PopMatrix();

                        break;

                    case K_HRC_SKELETON:
                        ComputeFieldBoundingBox(fSkeleton, fAnimation.frames[iCurrentFrameScroll], 
                                                ref p_min, ref p_max);

                        SetCameraAroundModel(ref p_min, ref p_max, 
                                             panX, panY, (float)(panZ + DIST),
                                             (float)alpha, (float)beta, (float)gamma, 1, 1, 1);

                        if (bShowGround)
                        {
                            GL.Disable(EnableCap.Lighting);
                            DrawGround();
                            DrawShadow(ref p_min, ref p_max);
                        }

                        SetLights();

                        DrawFieldSkeleton(fSkeleton, fAnimation.frames[iCurrentFrameScroll], bDListsEnable);

                        if (bShowLastFrameGhost)
                        {

                            GL.ColorMask(true, true, false, true);
                            if (iCurrentFrameScroll == 0)
                                DrawFieldSkeleton(fSkeleton, fAnimation.frames[fAnimation.nFrames - 1],
                                                  bDListsEnable);
                            else
                                DrawFieldSkeleton(fSkeleton, fAnimation.frames[iCurrentFrameScroll - 1],
                                                  bDListsEnable);

                            GL.ColorMask(true, true, true, true);
                        }

                        GL.Disable(EnableCap.Lighting);

                        if (bShowBones)
                        {
                            GL.Disable(EnableCap.DepthTest);

                            SkeletonRenderer.RenderFieldSkeletonBones(fSkeleton,
                                fAnimation.frames[iCurrentFrameScroll], 0, 1, 0, 1, 0, 0);

                            GL.Enable(EnableCap.DepthTest);
                        }

                        SelectFieldBoneAndPiece(fSkeleton, fAnimation.frames[iCurrentFrameScroll],
                                                SelectedBone, SelectedBonePiece);
                        break;

                    case K_AA_SKELETON:
                    case K_MAGIC_SKELETON:
                        //if (!bSkeleton.IsBattleLocation)
                        //int     animIndex = Int32.Parse(cbBattleAnimation.Items[cbBattleAnimation.SelectedIndex].ToString());

                        ComputeBattleBoundingBox(bSkeleton, bAnimationsPack.SkeletonAnimations[ianimIndex].frames[iCurrentFrameScroll],
                                                 ref p_min, ref p_max);

                        SetCameraAroundModel(ref p_min, ref p_max,
                                             panX, panY, (float)(panZ + DIST), 
                                             (float)alpha, (float)beta, (float)gamma, 1, 1, 1);

                        if (bShowGround)
                        {
                            GL.Disable(EnableCap.Lighting);
                            DrawGround();
                            DrawShadow(ref p_min, ref p_max);
                        }

                        SetLights();

                        tmpbFrame = new BattleFrame();
                        if (bSkeleton.wpModels.Count > 0 && bAnimationsPack.WeaponAnimations.Count > 0)
                        {
                            tmpbFrame = bAnimationsPack.WeaponAnimations[ianimIndex].frames[iCurrentFrameScroll];
                        }

                        DrawBattleSkeleton(bSkeleton, bAnimationsPack.SkeletonAnimations[ianimIndex].frames[iCurrentFrameScroll],
                                           tmpbFrame, ianimWeaponIndex, bDListsEnable);

                        if (bShowLastFrameGhost && !bSkeleton.IsBattleLocation)
                        {
                            GL.ColorMask(true, true, false, true);

                            if (iCurrentFrameScroll == 0)
                            {
                                DrawBattleSkeleton(bSkeleton, bAnimationsPack.SkeletonAnimations[ianimIndex].frames[iCurrentFrameScroll],
                                                   tmpbFrame, ianimWeaponIndex, bDListsEnable);
                            }
                            else
                                DrawBattleSkeleton(bSkeleton, bAnimationsPack.SkeletonAnimations[ianimIndex].frames[iCurrentFrameScroll - 1],
                                                   tmpbFrame, ianimWeaponIndex, bDListsEnable);

                            GL.ColorMask(true, true, true, true);
                        }

                        GL.Disable(EnableCap.Lighting);

                        if (bShowBones)
                        {
                            GL.Disable(EnableCap.DepthTest);

                            SkeletonRenderer.RenderBattleSkeletonBones(bSkeleton,
                                bAnimationsPack.SkeletonAnimations[ianimIndex].frames[iCurrentFrameScroll],
                                0, 1, 0, 1, 0, 0);

                            GL.Enable(EnableCap.DepthTest);
                        }

                        SelectBattleBoneAndModel(bSkeleton, bAnimationsPack.SkeletonAnimations[ianimIndex].frames[iCurrentFrameScroll],
                            tmpbFrame, ianimWeaponIndex, SelectedBone, SelectedBonePiece);
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

        public static void DrawShadow(ref Point3D p_min, ref Point3D p_max)
        {
            float ground_radius, sub_y, cx, cz;
            int numSegments, si;

            Point3D p_min_aux, p_max_aux;

            sub_y = p_max.y;
            p_min_aux = p_min;
            p_max_aux = p_max;
            p_min_aux.y = 0;
            p_max_aux.y = 0;

            cx = (p_min.x + p_max.x) / 2;
            cz = (p_min.z + p_max.z) / 2;
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
                    GL.Vertex3f((float)(ground_radius * Math.Sin(si * 2 * PI / numSegments) + cx), 0,
                               (float)(ground_radius * Math.Cos(si * 2 * PI / numSegments) + cz));
                }
            GL.End();

            GL.Enable(EnableCap.DepthTest);
            GL.Disable(EnableCap.Fog);

            // Draw underlying box (just depth)
            GL.ColorMask(false, false, false, false);
            GL.Color3f(1, 1, 1);
            GL.Begin(PrimitiveType.Quads);
                GL.Vertex3f(p_max.x, 0, p_max.z);
                GL.Vertex3f(p_max.x, 0, p_min.z);
                GL.Vertex3f(p_min.x, 0, p_min.z);
                GL.Vertex3f(p_min.x, 0, p_max.z);

                GL.Vertex3f(p_max.x, 0, p_max.z);
                GL.Vertex3f(p_max.x, sub_y, p_max.z);
                GL.Vertex3f(p_max.x, sub_y, p_min.z);
                GL.Vertex3f(p_max.x, 0, p_min.z);

                GL.Vertex3f(p_max.x, 0, p_min.z);
                GL.Vertex3f(p_max.x, sub_y, p_min.z);
                GL.Vertex3f(p_min.x, sub_y, p_min.z);
                GL.Vertex3f(p_min.x, 0, p_min.z);

                GL.Vertex3f(p_min.x, sub_y, p_max.z);
                GL.Vertex3f(p_min.x, 0, p_max.z);
                GL.Vertex3f(p_min.x, 0, p_min.z);
                GL.Vertex3f(p_min.x, sub_y, p_min.z);

                GL.Vertex3f(p_max.x, sub_y, p_max.z);
                GL.Vertex3f(p_max.x, 0, p_max.z);
                GL.Vertex3f(p_min.x, 0, p_max.z);
                GL.Vertex3f(p_min.x, sub_y, p_max.z);
            GL.End();
            GL.ColorMask(true, true, true, true);
        }

        public static void DrawPlane(ref double[] planeTransformation, ref Point3D planeOriginalPoint1,
                                                                       ref Point3D planeOriginalPoint2,
                                                                       ref Point3D planeOriginalPoint3,
                                                                       ref Point3D planeOriginalPoint4)
        {
            Point3D p1 = new Point3D();
            Point3D p2 = new Point3D();
            Point3D p3 = new Point3D();
            Point3D p4 = new Point3D();

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
                GL.Vertex3f(p1.x, p1.y, p1.z);
                GL.Vertex3f(p2.x, p2.y, p2.z);
                GL.Vertex3f(p3.x, p3.y, p3.z);
                GL.Vertex3f(p4.x, p4.y, p4.z);
            GL.End();
        }

        public static void DrawAxes(Control pbIn, int iFrame)
        {
            float letterWidth, letterHeight;
            float max_x, max_y, max_z;

            Point3D pX = new Point3D();
            Point3D pY = new Point3D();
            Point3D pZ = new Point3D();

            Point3D p_max = new Point3D();
            Point3D p_min = new Point3D();

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

            //max_x = Math.Abs(p_min.x) > Math.Abs(p_max.x) ? p_min.x : p_max.x;
            //max_y = Math.Abs(p_min.y) > Math.Abs(p_max.y) ? p_min.y : p_max.y;
            //max_z = Math.Abs(p_min.z) > Math.Abs(p_max.z) ? p_min.z : p_max.z;

            max_x = Math.Abs(p_max.x - p_min.x);
            max_y = Math.Abs(p_max.y - p_min.y);
            max_z = Math.Abs(p_max.z - p_min.z);

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
            pX.x = max_x;
            pX.y = 0;
            pX.z = 0;
            pX = GetProjectedCoords(pX);

            //  Get projected end of the Y axis
            pY.x = 0;
            pY.y = -max_y;
            pY.z = 0;
            pY = GetProjectedCoords(pY);

            //  Get projected end of the Z axis
            pZ.x = 0;
            pZ.y = 0;
            pZ.z = max_z;
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
                GL.Vertex2f(pX.x - letterWidth, pX.y - letterHeight);
                GL.Vertex2f(pX.x + letterWidth, pX.y + letterHeight);
                GL.Vertex2f(pX.x - letterWidth, pX.y + letterHeight);
                GL.Vertex2f(pX.x + letterWidth, pX.y - letterHeight);

                if (bSkeleton.IsBattleLocation)
                {
                    //  Draw Y
                    GL.Color3f(0, 0, 0);
                    GL.Vertex2f(pY.x - letterWidth, pY.y - letterHeight);
                    GL.Vertex2f(pY.x + letterWidth, pY.y + letterHeight);
                    GL.Vertex2f(pY.x - letterWidth, pY.y + letterHeight);
                    GL.Vertex2f(pY.x, pY.y);

                    //  Draw Z
                    GL.Color3f(0, 0, 0);
                    GL.Vertex2f(pZ.x + letterWidth, pZ.y + letterHeight);
                    GL.Vertex2f(pZ.x - letterWidth, pZ.y + letterHeight);

                    GL.Vertex2f(pZ.x + letterWidth, pZ.y + letterHeight);
                    GL.Vertex2f(pZ.x - letterWidth, pZ.y - letterHeight);

                    GL.Vertex2f(pZ.x + letterWidth, pZ.y - letterHeight);
                    GL.Vertex2f(pZ.x - letterWidth, pZ.y - letterHeight);
                }
                else 
                {
                    //  Draw Y
                    GL.Color3f(0, 0, 0);
                    GL.Vertex2f(pZ.x - letterWidth, pZ.y - letterHeight);
                    GL.Vertex2f(pZ.x + letterWidth, pZ.y + letterHeight);
                    GL.Vertex2f(pZ.x - letterWidth, pZ.y + letterHeight);
                    GL.Vertex2f(pZ.x, pZ.y);

                    //  Draw Z
                    GL.Color3f(0, 0, 0);
                    GL.Vertex2f(pY.x + letterWidth, pY.y + letterHeight);
                    GL.Vertex2f(pY.x - letterWidth, pY.y + letterHeight);

                    GL.Vertex2f(pY.x + letterWidth, pY.y + letterHeight);
                    GL.Vertex2f(pY.x - letterWidth, pY.y - letterHeight);

                    GL.Vertex2f(pY.x + letterWidth, pY.y - letterHeight);
                    GL.Vertex2f(pY.x - letterWidth, pY.y - letterHeight);
                }
            GL.End();

            GL.Enable(EnableCap.DepthTest);
        }




        //  ---------------------------------------------------------------------------------------------------
        //  ======================================= PEDITOR DRAW  =============================================
        //  ---------------------------------------------------------------------------------------------------
        public static void DrawPModelPolys(PModel Model)
        {
            int iGroupIdx, iPolyIdx, iVertIdx;

            GL.ShadeModel(ShadingModel.Flat);

            GL.PolygonMode(TriangleFace.Front, PolygonMode.Line);
            GL.PolygonMode(TriangleFace.Back, PolygonMode.Fill);
            GL.Enable(EnableCap.ColorMaterial);

            for (iGroupIdx = 0; iGroupIdx < Model.Header.numGroups; iGroupIdx++)
            {
                if (!Model.Groups[iGroupIdx].HiddenQ)
                {

                    // We will apply Group update values for Reposition/Resize/Rotate.
                    double[] rot_mat = new double[16];

                    GL.MatrixMode(MatrixMode.Modelview);
                    GL.PushMatrix();

                    GL.Translated(Model.Groups[iGroupIdx].repGroupX,
                                 Model.Groups[iGroupIdx].repGroupY,
                                 Model.Groups[iGroupIdx].repGroupZ);

                    BuildRotationMatrixWithQuaternionsXYZ(Model.Groups[iGroupIdx].rotGroupAlpha,
                                                          Model.Groups[iGroupIdx].rotGroupBeta,
                                                          Model.Groups[iGroupIdx].rotGroupGamma,
                                                          ref rot_mat);

                    GL.MultMatrixd(rot_mat);
                    GL.Scaled(Model.Groups[iGroupIdx].rszGroupX,
                             Model.Groups[iGroupIdx].rszGroupY,
                             Model.Groups[iGroupIdx].rszGroupZ);

                    for (iPolyIdx = Model.Groups[iGroupIdx].offsetPoly;
                         iPolyIdx < Model.Groups[iGroupIdx].offsetPoly + Model.Groups[iGroupIdx].numPoly;
                         iPolyIdx++)
                    {
                        GL.Color4f(Model.Pcolors[iPolyIdx].R / 255.0f, 
                                  Model.Pcolors[iPolyIdx].G / 255.0f, 
                                  Model.Pcolors[iPolyIdx].B / 255.0f, 
                                  Model.Pcolors[iPolyIdx].A / 255.0f);

                        GL.ColorMaterial(TriangleFace.FrontAndBack, ColorMaterialParameter.AmbientAndDiffuse);

                        GL.Begin(PrimitiveType.Triangles);
                        for (iVertIdx = 0; iVertIdx < 3; iVertIdx++)
                        {
                            //GL.Normal3f(Model.Normals[Model.Polys[iPolyIdx].Normals[iVertIdx]].x,
                            //           Model.Normals[Model.Polys[iPolyIdx].Normals[iVertIdx]].y,
                            //           Model.Normals[Model.Polys[iPolyIdx].Normals[iVertIdx]].z);

                            GL.Normal3f(Model.Normals[Model.NormalIndex[Model.Polys[iPolyIdx].Verts[iVertIdx] +
                                                                       Model.Groups[iGroupIdx].offsetVert]].x,
                                       Model.Normals[Model.NormalIndex[Model.Polys[iPolyIdx].Verts[iVertIdx] +
                                                                       Model.Groups[iGroupIdx].offsetVert]].y,
                                       Model.Normals[Model.NormalIndex[Model.Polys[iPolyIdx].Verts[iVertIdx] +
                                                                       Model.Groups[iGroupIdx].offsetVert]].z);

                            GL.Vertex3f(Model.Verts[Model.Polys[iPolyIdx].Verts[iVertIdx] + 
                                                   Model.Groups[iGroupIdx].offsetVert].x,
                                       Model.Verts[Model.Polys[iPolyIdx].Verts[iVertIdx] + 
                                                   Model.Groups[iGroupIdx].offsetVert].y,
                                       Model.Verts[Model.Polys[iPolyIdx].Verts[iVertIdx] + 
                                                   Model.Groups[iGroupIdx].offsetVert].z);
                        }
                        GL.End();
                    }

                    GL.PopMatrix();
                }

                // Let's try here to render the normals
                if (bShowVertexNormals || bShowFaceNormals)
                    ShowNormals(Model.Groups[iGroupIdx], Model.Polys, Model.Verts, 
                                Model.Normals, Model.NormalIndex);
            }

            //  GL.PopMatrix();  -- Commented in KimeraVB6
        }

        public static void DrawPModelMesh(PModel Model)
        {
            int iGroupIdx, iPolyIdx, iVertIdx;

            // -- Commented in KimeraVB6
            //  glMatrixMode GL_MODELVIEW
            //  glPushMatrix
            //  With obj
            //      glScalef .ResizeX, .ResizeY, .ResizeZ
            //      glRotatef .RotateAlpha, 1, 0, 0
            //      glRotatef .RotateBeta, 0, 1, 0
            //      glRotatef .RotateGamma, 0, 0, 1
            //      glTranslatef .RepositionX, .RepositionY, .RepositionZ
            //  End With

            GL.PolygonMode(TriangleFace.Front, PolygonMode.Line);
            GL.PolygonMode(TriangleFace.Back, PolygonMode.Line);
            GL.Color3f(0, 0, 0);

            for (iGroupIdx = 0; iGroupIdx < Model.Header.numGroups; iGroupIdx++)
            {

                if (!Model.Groups[iGroupIdx].HiddenQ)
                {

                    // We will apply Group update values for Reposition/Resize/Rotate.
                    double[] rot_mat = new double[16];

                    GL.MatrixMode(MatrixMode.Modelview);
                    GL.PushMatrix();

                    GL.Translated(Model.Groups[iGroupIdx].repGroupX,
                                 Model.Groups[iGroupIdx].repGroupY,
                                 Model.Groups[iGroupIdx].repGroupZ);

                    BuildRotationMatrixWithQuaternionsXYZ(Model.Groups[iGroupIdx].rotGroupAlpha,
                                                          Model.Groups[iGroupIdx].rotGroupBeta,
                                                          Model.Groups[iGroupIdx].rotGroupGamma,
                                                          ref rot_mat);

                    GL.MultMatrixd(rot_mat);
                    GL.Scaled(Model.Groups[iGroupIdx].rszGroupX,
                             Model.Groups[iGroupIdx].rszGroupY,
                             Model.Groups[iGroupIdx].rszGroupZ);

                    for (iPolyIdx = Model.Groups[iGroupIdx].offsetPoly; 
                         iPolyIdx < Model.Groups[iGroupIdx].offsetPoly + Model.Groups[iGroupIdx].numPoly;
                         iPolyIdx++)
                    {
                        GL.Begin(PrimitiveType.Triangles);
                        for (iVertIdx = 0; iVertIdx < 3; iVertIdx++)
                        {
                            GL.Vertex3f(Model.Verts[Model.Polys[iPolyIdx].Verts[iVertIdx] + 
                                                   Model.Groups[iGroupIdx].offsetVert].x,
                                       Model.Verts[Model.Polys[iPolyIdx].Verts[iVertIdx] + 
                                                   Model.Groups[iGroupIdx].offsetVert].y,
                                       Model.Verts[Model.Polys[iPolyIdx].Verts[iVertIdx] + 
                                                   Model.Groups[iGroupIdx].offsetVert].z);
                        }
                        GL.End();
                    }

                    GL.PopMatrix();
                }
            }

            //  GL.PopMatrix();      -- Commented in KimeraVB6
        }

        public static void KillPrecalculatedLighting(PModel Model, ref PairIB[] translationTableVertex)
        {
            int iVColorIdx;

            for (iVColorIdx = 0; iVColorIdx < Model.Header.numVerts; iVColorIdx++)
            {
                translationTableVertex[iVColorIdx].B = 1;
            }
        }

        public static void DrawPModelEditor(bool bEnableLighting, Control pbIn)
        {

            Point3D p_min = new Point3D();
            Point3D p_max = new Point3D();

            GL.Viewport(0, 0, pbIn.ClientRectangle.Width,
                             pbIn.ClientRectangle.Height);
            ClearPanel();

            SetCameraPModel(EditedPModel,
                            panXPE, panYPE, panZPE + DISTPE,
                            alphaPE, betaPE, gammaPE,
                            1, 1, 1);

            ConcatenateCameraModelView(repXPE, repYPE, repZPE,
                                       0, 0, 0,
                                       rszXPE, rszYPE, rszZPE);

            GL.MatrixMode(MatrixMode.Modelview);
            GL.PushMatrix();

            ComputePModelBoundingBox(EditedPModel, ref p_min, ref p_max);
            float modelDiameterNormalized = (-2 * ComputeSceneRadius(p_min, p_max)) / FrmPEditor.LIGHT_STEPS;

            if (bEnableLighting)
            {
                GL.Disable(EnableCap.Light0);
                GL.Disable(EnableCap.Light1);
                GL.Disable(EnableCap.Light2);
                GL.Disable(EnableCap.Light3);

                SetLighting(LightName.Light0, modelDiameterNormalized * FrmPEditor.iLightX,
                                            modelDiameterNormalized * FrmPEditor.iLightY,
                                            modelDiameterNormalized * FrmPEditor.iLightZ,
                                            1, 1, 1, false);
            }
            else GL.Disable(EnableCap.Lighting);

            // Sync lighting settings to modern renderer
            GLRenderer.LightingEnabled = bEnableLighting;
            if (bEnableLighting)
            {
                GLRenderer.LightPosition = new OpenTK.Mathematics.Vector3(
                    modelDiameterNormalized * FrmPEditor.iLightX,
                    modelDiameterNormalized * FrmPEditor.iLightY,
                    modelDiameterNormalized * FrmPEditor.iLightZ);
            }

            SetDefaultOGLRenderState();

            switch (drawMode)
            {
                case K_MESH:
                    DrawPModelWireframe(ref EditedPModel, true);
                    break;

                case K_PCOLORS:
                    GL.Enable(EnableCap.PolygonOffsetFill);
                    GL.PolygonOffset(1, 1);
                    DrawPModelPolygonColors(ref EditedPModel, true);
                    GL.Disable(EnableCap.PolygonOffsetFill);

                    DrawPModelWireframe(ref EditedPModel, true);
                    break;

                case K_VCOLORS:
                    DrawPModel(ref EditedPModel, ref tex_ids, true);
                    break;
            }

            GL.PopMatrix();
        }

        public static void DrawAxesPE(Control pbIn)
        {
            float letterWidth, letterHeight;
            float max_x, max_y, max_z;

            Point3D pX = new Point3D();
            Point3D pY = new Point3D();
            Point3D pZ = new Point3D();

            Point3D p_max = new Point3D();
            Point3D p_min = new Point3D();

            GL.Disable(EnableCap.Lighting);
            ComputePModelBoundingBox(EditedPModel, ref p_min, ref p_max);

            max_x = Math.Abs(p_min.x) > Math.Abs(p_max.x) ? p_min.x : p_max.x;
            max_y = Math.Abs(p_min.y) > Math.Abs(p_max.y) ? p_min.y : p_max.y;
            max_z = Math.Abs(p_min.z) > Math.Abs(p_max.z) ? p_min.z : p_max.z;

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
            pX.x = 2 * max_x;
            pX.y = 0;
            pX.z = 0;
            pX = GetProjectedCoords(pX);

            //  Get projected end of the Y axis
            pY.x = 0;
            pY.y = 2 * max_y;
            pY.z = 0;
            pY = GetProjectedCoords(pY);

            //  Get projected end of the Z axis
            pZ.x = 0;
            pZ.y = 0;
            pZ.z = 2 * max_z;
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
                GL.Vertex2f(pX.x - letterWidth, pX.y - letterHeight);
                GL.Vertex2f(pX.x + letterWidth, pX.y + letterHeight);
                GL.Vertex2f(pX.x - letterWidth, pX.y + letterHeight);
                GL.Vertex2f(pX.x + letterWidth, pX.y - letterHeight);

                //  Draw Y
                GL.Color3f(0, 0, 0);
                GL.Vertex2f(pY.x - letterWidth, pY.y - letterHeight);
                GL.Vertex2f(pY.x + letterWidth, pY.y + letterHeight);
                GL.Vertex2f(pY.x - letterWidth, pY.y + letterHeight);
                GL.Vertex2f(pY.x, pY.y);

                //  Draw Z
                GL.Color3f(0, 0, 0);
                GL.Vertex2f(pZ.x + letterWidth, pZ.y + letterHeight);
                GL.Vertex2f(pZ.x - letterWidth, pZ.y + letterHeight);

                GL.Vertex2f(pZ.x + letterWidth, pZ.y + letterHeight);
                GL.Vertex2f(pZ.x - letterWidth, pZ.y - letterHeight);

                GL.Vertex2f(pZ.x + letterWidth, pZ.y - letterHeight);
                GL.Vertex2f(pZ.x - letterWidth, pZ.y - letterHeight);
            GL.End();

            GL.Enable(EnableCap.DepthTest);
        }

        public static int GetEqualGroupVertices(PModel Model, int iActualVertIdx, ref List<int> lstVerts)
        {
            int iGroupIdx, iVertIdx;
            Point3D tmpUP3DVert;

            tmpUP3DVert = Model.Verts[iActualVertIdx];
            iGroupIdx = GetVertexGroup(Model, iActualVertIdx);

            for (iVertIdx = Model.Groups[iGroupIdx].offsetVert; 
                 iVertIdx < Model.Groups[iGroupIdx].offsetVert + Model.Groups[iGroupIdx].numVert;
                 iVertIdx++)
            {
                if (ComparePoints3D(Model.Verts[iVertIdx], tmpUP3DVert))
                {
                    lstVerts.Add(iVertIdx);
                    //  Debug.Print "Intended("; n_verts; ")"; Str$(vi)
                }
            }

            return lstVerts.Count;
        }

        public static int AddPaintVertex(ref PModel Model, int iGroupIdx, Point3D vPoint3D, Color vColor)
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

            Point3D tmpVert;
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

                tmpVert.x = Model.Verts[Model.Polys[iPolyIdx].Verts[iVertIdx] +
                                        Model.Groups[iGroupIdx].offsetVert].x;
                tmpVert.y = Model.Verts[Model.Polys[iPolyIdx].Verts[iVertIdx] +
                                        Model.Groups[iGroupIdx].offsetVert].y;
                tmpVert.z = Model.Verts[Model.Polys[iPolyIdx].Verts[iVertIdx] +
                                        Model.Groups[iGroupIdx].offsetVert].z;

                tmpColor = Color.FromArgb(255, bR, bG, bB);

                iPaintVertexResult = AddPaintVertex(ref Model, iGroupIdx, tmpVert, tmpColor) - Model.Groups[iGroupIdx].offsetVert;

                Model.Normals[Model.NormalIndex[iPaintVertexResult]] =
                    new Point3D(Model.Normals[Model.NormalIndex[Model.Polys[iPolyIdx].Verts[iVertIdx]]+
                                                                Model.Groups[iGroupIdx].offsetVert].x,
                                Model.Normals[Model.NormalIndex[Model.Polys[iPolyIdx].Verts[iVertIdx]] +
                                                                Model.Groups[iGroupIdx].offsetVert].y,
                                Model.Normals[Model.NormalIndex[Model.Polys[iPolyIdx].Verts[iVertIdx]] +
                                                                Model.Groups[iGroupIdx].offsetVert].z);

                if (bTextured)
                {
                    Model.TexCoords[Model.Groups[iGroupIdx].offsetTex + iPaintVertexResult] =
                        new Point2D(Model.TexCoords[Model.Polys[iPolyIdx].Verts[iVertIdx] + Model.Groups[iGroupIdx].offsetTex].x,
                                    Model.TexCoords[Model.Polys[iPolyIdx].Verts[iVertIdx] + Model.Groups[iGroupIdx].offsetTex].y);
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
                //  'Debug.Print "Vert(:", .Verts(vi), ",", Group, ")", obj.Verts(.Verts(vi) + obj.Groups(Group).offVert).x, obj.Verts(.Verts(vi) + obj.Groups(Group).offVert).y, obj.Verts(.Verts(vi) + obj.Groups(Group).offVert).z
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

