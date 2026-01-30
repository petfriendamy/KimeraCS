using OpenTK.Graphics.OpenGL.Compatibility;
using OpenTK.Mathematics;
using KimeraCS.Core;


namespace KimeraCS.Rendering
{
    using static VisualizationHelpers;

    using static FF7Skeleton;
    using static FF7PModel;
    using static Lighting;
    using static Utils;

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

        private static void ShowNormals(PGroup Group, PPolygon[] Polys, Vector3[] Verts,
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

        private static void DrawPModel(ref PModel Model, ref uint[] tex_ids, bool HideHiddenGroupsQ,
                                      RenderingOptions renderOptions)
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
            if (renderOptions.NormalsDisplayMode != NormalsDisplayMode.None)
            {
                for (int g = 0; g < Model.Header.numGroups; g++)
                {
                    ShowNormals(Model.Groups[g], Model.Polys, Model.Verts,
                                Model.Normals, Model.NormalIndex,
                                renderOptions.NormalsDisplayMode,
                                renderOptions.NormalsScale,
                                renderOptions.NormalsColor);
                }
            }

            // Restore original matrices
            GLRenderer.ProjectionMatrix = savedProjection;
            GLRenderer.ViewMatrix = savedView;
            GLRenderer.ModelMatrix = savedModel;
        }

        private static void DrawPModelWireframe(ref PModel Model, bool HideHiddenGroupsQ)
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

        private static void DrawPModelPolygonColors(ref PModel Model, bool HideHiddenGroupsQ)
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

        private static void DrawPModelDLists(ref PModel Model, ref uint[] tex_ids,
                                            RenderingOptions renderOptions)
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
            if (renderOptions.NormalsDisplayMode != NormalsDisplayMode.None)
            {
                for (int g = 0; g < Model.Header.numGroups; g++)
                {
                    ShowNormals(Model.Groups[g], Model.Polys, Model.Verts,
                                Model.Normals, Model.NormalIndex,
                                renderOptions.NormalsDisplayMode,
                                renderOptions.NormalsScale,
                                renderOptions.NormalsColor);
                }
            }

            // Restore original matrices
            GLRenderer.ProjectionMatrix = savedProjection;
            GLRenderer.ViewMatrix = savedView;
            GLRenderer.ModelMatrix = savedModel;
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


        //  ---------------------------------------------------------------------------------------------------
        //  ====================================== SKELETON DRAW  =============================================
        //  ---------------------------------------------------------------------------------------------------

        /// <summary>
        /// Draws a single model attached to a unified bone.
        /// Handles both field (quaternion rotation) and battle (Euler rotation) models.
        /// </summary>
        private static void DrawUnifiedBoneModel(UnifiedBoneModel boneModel, uint[] globalTexIds,
                                                  bool useQuaternion, bool bDListsEnable,
                                                  RenderingOptions options)
        {
            double[] rot_mat = new double[16];
            var model = boneModel.Model;

            // Build texture IDs from the model's own textures if available (field models),
            // otherwise use global texture IDs (battle models)
            uint[] texIds;
            if (boneModel.Textures != null && boneModel.Textures.Count > 0)
            {
                texIds = new uint[boneModel.Textures.Count];
                for (int i = 0; i < boneModel.Textures.Count; i++)
                {
                    texIds[i] = boneModel.Textures[i].texID;
                }
            }
            else
            {
                texIds = globalTexIds;
            }

            GL.MatrixMode(MatrixMode.Modelview);
            GL.PushMatrix();

            // Apply model position offset
            GL.Translated(model.repositionX, model.repositionY, model.repositionZ);

            // Apply rotation - quaternion for field, Euler for battle
            if (useQuaternion)
            {
                BuildMatrixFromQuaternion(model.rotationQuaternion, ref rot_mat);
                GL.MultMatrixd(rot_mat);
            }
            else
            {
                GL.Rotated(model.rotateAlpha, 1, 0, 0);
                GL.Rotated(model.rotateBeta, 0, 1, 0);
                GL.Rotated(model.rotateGamma, 0, 0, 1);
            }

            // Apply model scale
            GL.Scaled(model.resizeX, model.resizeY, model.resizeZ);

            SetDefaultOGLRenderState();

            // Draw the model
            var tmpModel = model;
            if (bDListsEnable)
                DrawPModelDLists(ref tmpModel, ref texIds, options);
            else
                DrawPModel(ref tmpModel, ref texIds, false, options);

            GL.PopMatrix();
        }

        /// <summary>
        /// Draws all models attached to a unified bone.
        /// </summary>
        private static void DrawUnifiedBone(UnifiedBone bone, uint[] texIds,
                                            bool useQuaternion, bool bDListsEnable,
                                            RenderingOptions options)
        {
            if (!bone.HasModel) return;

            GL.MatrixMode(MatrixMode.Modelview);
            GL.PushMatrix();

            // Apply bone-level scale
            GL.Scaled(bone.Scale.X, bone.Scale.Y, bone.Scale.Z);

            // Draw each model attached to this bone
            foreach (var boneModel in bone.Models)
            {
                DrawUnifiedBoneModel(boneModel, texIds, useQuaternion, bDListsEnable, options);
            }

            GL.PopMatrix();
        }

        /// <summary>
        /// Draws a unified skeleton with its attached models.
        /// Works for both field and battle skeletons through the unified format.
        /// </summary>
        private static void DrawUnifiedSkeleton(UnifiedSkeleton skeleton, UnifiedFrame frame,
                                               uint[] texIds, bool bDListsEnable,
                                               RenderingOptions options)
        {
            if (skeleton == null || frame == null || skeleton.Bones == null || skeleton.Bones.Count == 0)
                return;

            double[] rot_mat = new double[16];
            bool useQuaternion = skeleton.SourceType == SkeletonSourceType.Field;
            float boneDir = skeleton.BoneDirection == BoneDirection.NegativeZ ? -1f : 1f;

            // Stack for hierarchical transforms
            var parentStack = new System.Collections.Generic.Stack<int>();

            GL.MatrixMode(MatrixMode.Modelview);
            GL.PushMatrix();

            // Apply root transform
            GL.Translated(frame.RootTranslation.X, frame.RootTranslation.Y, frame.RootTranslation.Z);

            BuildRotationMatrixWithQuaternions(
                frame.RootRotation.Alpha,
                frame.RootRotation.Beta,
                frame.RootRotation.Gamma,
                ref rot_mat);
            GL.MultMatrixd(rot_mat);

            parentStack.Push(-1);

            for (int i = 0; i < skeleton.Bones.Count; i++)
            {
                var bone = skeleton.Bones[i];

                // Pop matrices until we find the parent
                while (parentStack.Count > 0 && parentStack.Peek() != bone.ParentIndex)
                {
                    GL.PopMatrix();
                    parentStack.Pop();
                }

                GL.PushMatrix();
                parentStack.Push(i);

                // Apply bone rotation
                if (i < frame.BoneRotations.Count)
                {
                    var rot = frame.BoneRotations[i];
                    BuildRotationMatrixWithQuaternions(rot.Alpha, rot.Beta, rot.Gamma, ref rot_mat);
                    GL.MultMatrixd(rot_mat);
                }

                // Draw this bone's models
                DrawUnifiedBone(bone, texIds, useQuaternion, bDListsEnable, options);

                // Translate along bone for children
                GL.Translated(0, 0, boneDir * bone.Length);
            }

            // Pop remaining matrices
            while (parentStack.Count > 0)
            {
                GL.PopMatrix();
                parentStack.Pop();
            }

            GL.PopMatrix();
        }

        /// <summary>
        /// Draws weapons for a unified skeleton (battle only).
        /// </summary>
        private static void DrawUnifiedWeapon(UnifiedSkeleton skeleton, UnifiedFrame weaponFrame,
                                             int weaponIndex, uint[] texIds, bool bDListsEnable,
                                             RenderingOptions options)
        {
            if (skeleton == null || skeleton.Weapons == null || weaponIndex < 0 ||
                weaponIndex >= skeleton.Weapons.Count || weaponFrame == null)
                return;

            double[] rot_mat = new double[16];
            var weapon = skeleton.Weapons[weaponIndex];

            GL.MatrixMode(MatrixMode.Modelview);
            GL.PushMatrix();

            // Apply weapon frame transform
            GL.Translated(weaponFrame.RootTranslation.X, weaponFrame.RootTranslation.Y, weaponFrame.RootTranslation.Z);

            BuildRotationMatrixWithQuaternions(
                weaponFrame.RootRotation.Alpha,
                weaponFrame.RootRotation.Beta,
                weaponFrame.RootRotation.Gamma,
                ref rot_mat);
            GL.MultMatrixd(rot_mat);

            GL.PushMatrix();

            // Apply weapon model transform
            GL.Translated(weapon.repositionX, weapon.repositionY, weapon.repositionZ);

            GL.Rotated(weapon.rotateAlpha, 1, 0, 0);
            GL.Rotated(weapon.rotateBeta, 0, 1, 0);
            GL.Rotated(weapon.rotateGamma, 0, 0, 1);

            GL.Scaled(weapon.resizeX, weapon.resizeY, weapon.resizeZ);

            SetDefaultOGLRenderState();

            // Draw the weapon model
            var tmpModel = weapon;
            if (bDListsEnable)
                DrawPModelDLists(ref tmpModel, ref texIds, options);
            else
                DrawPModel(ref tmpModel, ref texIds, false, options);

            GL.PopMatrix();
            GL.PopMatrix();
        }


        /// <summary>
        /// Draws the current skeleton/model using the provided context.
        /// </summary>
        public static void DrawSkeletonModel(RenderingContext ctx)
        {
            double[] rot_mat = new double[16];

            Vector3 p_min = new Vector3();
            Vector3 p_max = new Vector3();

            var modelData = ctx.ModelData;
            if (modelData != null)
            {
                PModel pModel = modelData.PModel;
                uint[] texIds = modelData.TextureIds;
                UnifiedSkeleton? skeleton = modelData.Skeleton;
                UnifiedAnimation? fieldAnim = modelData.Animation;
                UnifiedAnimationPack? battleAnims = modelData.AnimationPack;

                try
                {
                    switch (ctx.ModelType)
                    {
                        case ModelType.ImportedModel:
                        case ModelType.PFieldModel:
                        case ModelType.PBattleModel:
                        case ModelType.PMagicModel:
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

                            DrawPModel(ref pModel, ref texIds, false, ctx.Options);

                            GL.PopMatrix();

                            break;

                        case ModelType.HRCSkeleton:
                            if (skeleton != null && fieldAnim != null)
                            {
                                skeleton.ComputeBoundingBox(fieldAnim.Frames[ctx.Animation.CurrentFrame],
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

                                var currentFrame = GetCurrentFrame(0, ctx.Animation.CurrentFrame);
                                if (currentFrame != null)
                                {
                                    DrawUnifiedSkeleton(skeleton, currentFrame, texIds, ctx.Options.EnableDisplayLists, ctx.Options);

                                    if (ctx.Options.ShowLastFrameGhost)
                                    {
                                        GL.ColorMask(true, true, false, true);
                                        int ghostFrameIndex = ctx.Animation.CurrentFrame == 0
                                            ? fieldAnim.FrameCount - 1
                                            : ctx.Animation.CurrentFrame - 1;
                                        var ghostFrame = GetCurrentFrame(0, ghostFrameIndex);
                                        if (ghostFrame != null)
                                        {
                                            DrawUnifiedSkeleton(skeleton, ghostFrame, texIds, ctx.Options.EnableDisplayLists, ctx.Options);
                                            GL.ColorMask(true, true, true, true);
                                        }
                                    }

                                    GL.Disable(EnableCap.Lighting);

                                    if (ctx.Options.ShowBones)
                                    {
                                        GL.Disable(EnableCap.DepthTest);
                                        SkeletonRenderer.RenderSkeletonBones(skeleton, currentFrame,
                                            0, 1, 0, 1, 0, 0);
                                        GL.Enable(EnableCap.DepthTest);
                                    }

                                    skeleton.SelectBoneAndModel(currentFrame,
                                                                ctx.Selection.SelectedBone, ctx.Selection.SelectedBonePiece);
                                }
                            }
                            break;

                        case ModelType.AASkeleton:
                        case ModelType.MagicSkeleton:
                            if (skeleton != null && battleAnims != null)
                            {
                                skeleton.ComputeBoundingBox(battleAnims.SkeletonAnimations[ctx.Animation.AnimationIndex].Frames[ctx.Animation.CurrentFrame],
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

                                // Use unified skeleton format for drawing
                                var currentFrame = GetCurrentFrame(ctx.Animation.AnimationIndex, ctx.Animation.CurrentFrame);
                                if (currentFrame != null)
                                {
                                    // Use skeleton's texture IDs for battle models
                                    var battleTexIds = skeleton.TextureIDs;
                                    DrawUnifiedSkeleton(skeleton, currentFrame, battleTexIds, ctx.Options.EnableDisplayLists, ctx.Options);

                                    // Draw weapon using unified format
                                    if (skeleton.Weapons != null && skeleton.Weapons.Count > 0)
                                    {
                                        var weaponFrame = GetCurrentWeaponFrame(ctx.Animation.AnimationIndex, ctx.Animation.CurrentFrame);
                                        if (weaponFrame != null)
                                        {
                                            DrawUnifiedWeapon(skeleton, weaponFrame, ctx.Animation.WeaponAnimationIndex,
                                                battleTexIds, ctx.Options.EnableDisplayLists, ctx.Options);
                                        }
                                    }

                                    if (ctx.Options.ShowLastFrameGhost && !skeleton.IsBattleLocation && modelData.AnimationPack != null)
                                    {
                                        GL.ColorMask(true, true, false, true);
                                        int ghostFrameIndex = ctx.Animation.CurrentFrame == 0
                                            ? modelData.AnimationPack.SkeletonAnimations[ctx.Animation.AnimationIndex].FrameCount - 1
                                            : ctx.Animation.CurrentFrame - 1;
                                        var ghostFrame = GetCurrentFrame(ctx.Animation.AnimationIndex, ghostFrameIndex);
                                        if (ghostFrame != null)
                                        {
                                            DrawUnifiedSkeleton(skeleton, ghostFrame, battleTexIds, ctx.Options.EnableDisplayLists, ctx.Options);
                                            GL.ColorMask(true, true, true, true);
                                        }
                                    }

                                    GL.Disable(EnableCap.Lighting);

                                    if (ctx.Options.ShowBones)
                                    {
                                        GL.Disable(EnableCap.DepthTest);
                                        SkeletonRenderer.RenderSkeletonBones(skeleton, currentFrame,
                                            0, 1, 0, 1, 0, 0);
                                        GL.Enable(EnableCap.DepthTest);
                                    }

                                    {
                                        var selectWeaponFrame = GetCurrentWeaponFrame(ctx.Animation.AnimationIndex, ctx.Animation.CurrentFrame);
                                        skeleton.SelectBoneAndModel(currentFrame,
                                            ctx.Selection.SelectedBone, ctx.Selection.SelectedBonePiece,
                                            selectWeaponFrame, ctx.Animation.WeaponAnimationIndex);
                                    }
                                }
                            }
                            break;
                    }
                }
                catch (Exception ex)
                {
                    throw new KimeraException("Error Drawing current model.", ex);
                }
            }
            else
            {
                throw new KimeraException("No model data provided!");
            }
        }



        //  ---------------------------------------------------------------------------------------------------
        //  ============================ DRAWING EXTENSIONS (GROUND/SHADOW)  ==================================
        //  ---------------------------------------------------------------------------------------------------
        private static void DrawGround()
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

        private static void DrawShadow(ref Vector3 p_min, ref Vector3 p_max)
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
            SetBlendMode(BlendMode.Average);

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

            SetBlendMode(BlendMode.Average);

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
                case ModelType.HRCSkeleton:
                    if (skeleton != null && animation != null)
                        skeleton.ComputeBoundingBox(animation.Frames[iFrame],
                                                    ref p_min, ref p_max);
                    break;

                case ModelType.AASkeleton:
                case ModelType.MagicSkeleton:
                    if (skeleton != null && animationPack != null)
                        skeleton.ComputeBoundingBox(animationPack.SkeletonAnimations[ianimIndex].Frames[iFrame],
                                                    ref p_min, ref p_max);

                    break;

                case ModelType.ImportedModel:
                case ModelType.PBattleModel:
                case ModelType.PFieldModel:
                case ModelType.PMagicModel:
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

                if (skeleton != null && skeleton.IsBattleLocation)
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

                if (skeleton != null && skeleton.IsBattleLocation)
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
                case Core.DrawMode.Mesh:
                    DrawPModelWireframe(ref model, true);
                    break;

                case Core.DrawMode.PolygonColors:
                    GL.Enable(EnableCap.PolygonOffsetFill);
                    GL.PolygonOffset(1, 1);
                    DrawPModelPolygonColors(ref model, true);
                    GL.Disable(EnableCap.PolygonOffsetFill);

                    DrawPModelWireframe(ref model, true);
                    break;

                case Core.DrawMode.VertexColors:
                    DrawPModel(ref model, ref texIds, true, RenderingOptions.Default);
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

        private static int AddPaintVertex(ref PModel Model, int iGroupIdx, Vector3 vPoint3D, Color vColor)
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

        private static int PaintVertex(ref PModel Model, int iGroupIdx, int iPolyIdx, int iVertIdx,
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

