namespace KimeraCS.Core
{
    using static FF7FieldSkeleton;
    using static FF7FieldAnimation;

    using static FF7BattleSkeleton;
    using static FF7BattleAnimationsPack;

    using static FF7FieldRSDResource;
    using static FF7PModel;
    using static FF7TMDModel;

    using static FF7TEXTexture;

    using static Utils;

    public static class FF7Skeleton
    {
        // The currently loaded model type
        public static ModelType modelType = ModelType.None;

        public static UnifiedSkeleton? skeleton;
        public static UnifiedAnimation? animation;           // Field animations
        public static UnifiedAnimationPack? animationPack;   // Battle animations

        public static PModel fPModel;
        public static TMDModel mTMDModel;

        public static bool IsTMDModel;
        public static bool IsRSDResource;

        public static bool bLoaded;

        /// <summary>
        /// Populate the unified skeleton and animation from the loaded FF7 format data.
        /// Call this after loading a skeleton to enable unified rendering.
        /// </summary>
        /// <param name="fSkeleton">The field skeleton to unify.</param>
        /// <param name="fAnimation">The field animation to unify.</param>
        private static void PopulateUnifiedFormat(FieldSkeleton fSkeleton, FieldAnimation fAnimation)
        {
            skeleton = new UnifiedSkeleton(fSkeleton);
            animation = new UnifiedAnimation(fAnimation);
            animationPack = null;
        }

        /// <summary>
        /// Populate the unified skeleton and animation from the loaded FF7 format data.
        /// Call this after loading a skeleton to enable unified rendering.
        /// </summary>
        /// <param name="bSkeleton">The batt;e skeleton to unify.</param>
        /// <param name="fAnimation">The battle animations pack to unify.</param>
        private static void PopulateUnifiedFormat(BattleSkeleton bSkeleton, BattleAnimationsPack bAnimationsPack, ModelType type)
        {
            skeleton = new UnifiedSkeleton(bSkeleton);
            animation = null;
            animationPack = new UnifiedAnimationPack(bAnimationsPack, bSkeleton.nBones);
        }

        /// <summary>
        /// Get the current unified frame for rendering.
        /// </summary>
        public static UnifiedFrame? GetCurrentFrame(int animationIndex, int frameIndex)
        {
            if (animation != null)
            {
                return animation.GetFrame(frameIndex);
            }

            if (animationPack != null)
            {
                var anim = animationPack.GetSkeletonAnimation(animationIndex);
                return anim?.GetFrame(frameIndex);
            }

            return null;
        }

        /// <summary>
        /// Get the current weapon frame for rendering (battle only).
        /// </summary>
        public static UnifiedFrame? GetCurrentWeaponFrame(int animationIndex, int frameIndex)
        {
            if (animationPack == null)
                return null;

            var anim = animationPack.GetWeaponAnimation(animationIndex);
            return anim?.GetFrame(frameIndex);
        }

        //
        // Global Skeleton/Model functions/procedures
        //
        public static int LoadSkeleton(string strFileName, bool loadGeometryQ,
                                       bool isLimitBreak, bool ignoreMissingPFiles,
                                       bool repairPolys, bool removeTextureCoords)
        {
            int iloadSkeletonResult = 1;

            try
            {

                // First we destroy the previous loaded Field Skeleton.
                if (bLoaded && ((int)modelType >= 3 && (int)modelType <= 5))
                {
                    if (DestroySkeleton() != 1) iloadSkeletonResult = -2;
                }

                modelType = GetSkeletonType(strFileName);

                if ((int)modelType >= 3 && (int)modelType <= 5)
                {
                    // LOAD Skeleton
                    // We load the Field Skeleton into memory.
                    switch (modelType)
                    {
                        case ModelType.HRCSkeleton:
                            // Field Skeleton (.hrc)
                            var fSkeleton = new FieldSkeleton(strFileName, loadGeometryQ,
                                                          ignoreMissingPFiles, repairPolys,
                                                          removeTextureCoords);
                            var fAnimation = new FieldAnimation(fSkeleton, "DUMMY.A", false);
                            PopulateUnifiedFormat(fSkeleton, fAnimation);
                            break;

                        case ModelType.AASkeleton:
                            // Battle Skeleton (aa)
                            var bSkeleton = new BattleSkeleton(strFileName, isLimitBreak, true, repairPolys,
                                                           removeTextureCoords);

                            // Normally we will have the ??DA file with the Animation Pack.
                            // Location Battle Models has NOT ??DA file.
                            // But editing models, it is possible we work without it. So, we will make something
                            // similiar as we did with Field Models, but we will check if ??DA file for the model exists.
                            var bAnimationsPack = new BattleAnimationsPack(bSkeleton, modelType, strFileName);
                            PopulateUnifiedFormat(bSkeleton, bAnimationsPack, modelType);
                            break;

                        case ModelType.MagicSkeleton:
                            // Magic Skeleton (.d)
                            bSkeleton = new BattleSkeleton(strFileName, true, repairPolys, removeTextureCoords);

                            // Normally we will have the *.A00 file with the Animation Pack.
                            // But editing models, it is possible we work without it. So, we will make something
                            // similiar as we did with Field Models, but we will check if *.A00 file for the model exists.
                            bAnimationsPack = new BattleAnimationsPack(bSkeleton, modelType, strFileName);

                            break;
                    }
                    bLoaded = true;
                }
                else
                {
                    iloadSkeletonResult = 0;  // No known skeleton
                }
            }
            catch (Exception ex)
            {
                throw new KimeraException("Error loading skeleton.", ex);
            }

            return iloadSkeletonResult;
        }

        public static int LoadSkeletonFromDB(string strFileName, 
                                                  string strAnimFileName, 
                                                  bool loadGeometryQ,
                                                  bool isLimitBreak,
                                                  bool ignoreMissingPFiles,
                                                  bool repairPolys,
                                                  bool removeTextureCoords)
        {
            int iloadSkeletonResult = 1;

            try
            {

                // First we destroy the previous loaded Field Skeleton.
                if (bLoaded && ((int)modelType >= 3 && (int)modelType <= 5))
                {
                    if (DestroySkeleton() != 1) iloadSkeletonResult = -2;
                }

                modelType = GetSkeletonType(strFileName);

                if ((int)modelType >= 3 && (int)modelType <= 5)
                {
                    // LOAD Skeleton
                    switch (modelType)
                    {
                        case ModelType.HRCSkeleton:
                            // We load the Field Skeleton into memory.

                            // Field Skeleton (.hrc)
                            var fSkeleton = new FieldSkeleton(strFileName, loadGeometryQ, ignoreMissingPFiles,
                                                          repairPolys, removeTextureCoords);
                            var result = LoadAnimationFromDB(fSkeleton, strAnimFileName);
                            if (result.Item1 != null)
                            {
                                PopulateUnifiedFormat(fSkeleton, (FieldAnimation)result.Item1);
                                iloadSkeletonResult = result.Item2;
                            }
                            break;

                        case ModelType.AASkeleton:
                            // Battle Skeleton (aa)
                            var bSkeleton = new BattleSkeleton(strFileName, isLimitBreak, true, repairPolys,
                                                           removeTextureCoords);

                            // Normally we will have the ??DA file with the Animation Pack.
                            // Location Battle Models has NOT ??DA file.
                            // But editing models, it is possible we work without it. So, we will make something
                            // similiar as we did with Field Models, but we will check if ??DA file for the model exists.
                            var bAnimationsPack = new BattleAnimationsPack(bSkeleton, modelType, strFileName);
                            PopulateUnifiedFormat(bSkeleton, bAnimationsPack, modelType);
                            break;

                        case ModelType.MagicSkeleton:
                            // Magic Skeleton (.d)
                            bSkeleton = new BattleSkeleton(strFileName, true, repairPolys, removeTextureCoords);

                            // Normally we will have the ??DA file with the Animation Pack.
                            // Location Battle Models has NOT ??DA file.
                            // But editing models, it is possible we work without it. So, we will make something
                            // similiar as we did with Field Models, but we will check if ??DA file for the model exists.
                            bAnimationsPack = new BattleAnimationsPack(bSkeleton, modelType, strFileName);
                            PopulateUnifiedFormat(bSkeleton, bAnimationsPack, modelType);
                            break;
                    }

                    bLoaded = true;
                }
                else
                {
                    iloadSkeletonResult = 0;  // No known skeleton
                }
            }
            catch (Exception ex)
            {
                strGlobalExceptionMessage = ex.Message;

                iloadSkeletonResult = -1;  // Error loading skeleton
            }

            return iloadSkeletonResult;
        }

        public static int DestroySkeleton()
        {
            int iDestroySkeletonResult = 1;

            skeleton = null;
            animation = null;
            animationPack = null;

            bLoaded = false;
            return iDestroySkeletonResult;
        }

        public static ModelType GetSkeletonType(string strFileName)
        {
            ModelType iSkeletonType = ModelType.None;
            string tmpFileName;

            if (strFileName.Length > 0)
            {
                switch (Path.GetExtension(strFileName).ToUpper())
                {
                    case ".HRC":
                        iSkeletonType = ModelType.HRCSkeleton;
                        break;

                    case "":
                        tmpFileName = Path.GetFileNameWithoutExtension(strFileName).ToUpper();
                        if (tmpFileName.Length > 2)
                        {
                            if (tmpFileName[tmpFileName.Length - 1] == 'A' && tmpFileName[tmpFileName.Length - 2] == 'A')
                                iSkeletonType = ModelType.AASkeleton;
                        }
                        break;

                    case ".D":
                        iSkeletonType = ModelType.MagicSkeleton;
                        break;
                }
            }

            return iSkeletonType;
        }

        public static int LoadRSDResourceModel(string strRSDFolder, string strRSDName,
                                               bool ignoreMissingPFiles, bool repairPolys,
                                               bool removeTextureCoords)
        {
            int iLoadRSDResourceModelResult = 0;
            string strfAnimation = "";

            FieldBone tmpfBone;
            FieldRSDResource tmpfRSDResource;
            List<TEX> textures_pool = [];

            try
            {
                // Create fSkeleton with 1 bone
                var fSkeleton = new FieldSkeleton()
                {
                    nBones = 1,
                    fileName = strRSDName,
                    name = strRSDName,

                    bones = new List<FieldBone>(),
                };

                // Load RSD Resource
                tmpfBone = new FieldBone()
                {
                    len = 1,
                    joint_f = "null",
                    joint_i = "root",
                    nResources = 1,
                    fRSDResources = new List<FieldRSDResource>(),

                    resizeX = 1,
                    resizeY = 1,
                    resizeZ = 1,
                };

                tmpfRSDResource = new FieldRSDResource(strRSDName, ref textures_pool, strRSDFolder,
                                                       ignoreMissingPFiles, repairPolys, removeTextureCoords);
                tmpfBone.fRSDResources.Add(tmpfRSDResource);

                fSkeleton.bones.Add(tmpfBone);

                // Create the Dummy animation (normally individual RSD Resource has not own animation)
                var fAnimation = new FieldAnimation(fSkeleton,
                                                strRSDFolder + "\\" + strfAnimation,
                                                false);

                modelType = ModelType.HRCSkeleton;
                IsRSDResource = true;
                bLoaded = true;

                // Populate unified format for rendering
                PopulateUnifiedFormat(fSkeleton, fAnimation);
            }
            catch (Exception ex)
            {
                throw new FileLoadException("There has been some error loading RSD Resource: " + strRSDName + ".",
                                            strRSDName, ex);
            }

            return iLoadRSDResourceModelResult;
        }



        ////////////////////////////////////////////////////////////////////////////////////////////////
        // Global Animation functions/procedures
        ////////////////////////////////////////////////////////////////////////////////////////////////
        public static int WriteAnimation(string strFileName)
        {
            int isaveAnimationResult = 0;

            try
            {
                switch (modelType)
                {
                    case ModelType.HRCSkeleton:
                        var fAnimation = animation?.ToFieldAnimation();
                        if (fAnimation != null)
                        {
                            WriteFieldAnimation((FieldAnimation)fAnimation, strFileName);
                            isaveAnimationResult = 1;
                        }
                        break;

                    case ModelType.AASkeleton:
                    case ModelType.MagicSkeleton:
                        if (animationPack != null)
                        {
                            var bAnimationsPack = animationPack.ToBattleAnimationsPack();
                            isaveAnimationResult = WriteBattleAnimationsPack(ref bAnimationsPack, modelType, strFileName);
                        }
                        break;
                }
            }
            catch (Exception ex)
            {
                strGlobalExceptionMessage = ex.Message;

                isaveAnimationResult = -1;
            }

            return isaveAnimationResult;
        }

        public static bool NumAnimFramesIsOne(int index)
        {
            bool iNumAnimFramesIsOneResult = false;
            switch (modelType)
            {
                case ModelType.HRCSkeleton:
                    if (animation?.FrameCount == 1) iNumAnimFramesIsOneResult = true;
                    break;

                case ModelType.AASkeleton:
                case ModelType.MagicSkeleton:
                    if (animationPack?.SkeletonAnimations[index].FrameCount == 1) iNumAnimFramesIsOneResult = true;
                    break;
            }

            return iNumAnimFramesIsOneResult;
        }

        public static int ReadFrameData(string strFileName, bool bMerge)
        {
            int iinputFrameData = 0;
            if (skeleton != null)
            {
                try
                {
                    switch (modelType)
                    {
                        case ModelType.HRCSkeleton:
                            if (animation != null)
                            {
                                var fSkeleton = skeleton.ToFieldSkeleton();
                                var fAnimation = animation.ToFieldAnimation();
                                ReadFieldFrameData(fSkeleton, ref fAnimation, strFileName, bMerge);
                                animation = new UnifiedAnimation(fAnimation);
                                iinputFrameData = 1;
                            }
                            break;

                        case ModelType.AASkeleton:
                        case ModelType.MagicSkeleton:
                            //iinputFrameData = WriteBattleFrameDataPack(ref bAnimationsPack, strFileName);
                            break;
                    }
                }
                catch (Exception ex)
                {
                    strGlobalExceptionMessage = ex.Message;

                    iinputFrameData = -1;
                }
            }
            return iinputFrameData;
        }

        public static int WriteFrameData(string strFileName)
        {
            int ioutputFrameData = 0;
            if (skeleton != null)
            {
                try
                {
                    switch (modelType)
                    {
                        case ModelType.HRCSkeleton:
                            if (animation != null)
                            {
                                var fSkeleton = skeleton.ToFieldSkeleton();
                                var fAnimation = animation.ToFieldAnimation();
                                WriteFieldFrameData(fSkeleton, fAnimation, strFileName);
                                animation = new UnifiedAnimation(fAnimation);
                                ioutputFrameData = 1;
                            }
                            break;

                        case ModelType.AASkeleton:
                        case ModelType.MagicSkeleton:
                            //ioutputFrameData = WriteBattleAnimationsPack(ref bAnimationsPack, strFileName);
                            break;
                    }
                }
                catch (Exception ex)
                {
                    strGlobalExceptionMessage = ex.Message;

                    ioutputFrameData = -1;
                }
            }
            return ioutputFrameData;
        }

        public static int ReadFrameDataSelective(string strFileName)
        {
            int iinputFrameData = 0;
            if (skeleton != null)
            {
                try
                {
                    switch (modelType)
                    {
                        case ModelType.HRCSkeleton:
                            if (animation != null)
                            {
                                var fSkeleton = skeleton.ToFieldSkeleton();
                                var fAnimation = animation.ToFieldAnimation();
                                ReadFieldFrameDataSelective(fSkeleton, ref fAnimation, strFileName);
                                animation = new UnifiedAnimation(fAnimation);
                                iinputFrameData = 1;
                            }
                            break;

                        case ModelType.AASkeleton:
                        case ModelType.MagicSkeleton:
                            //iinputFrameData = WriteBattleFrameDataPack(ref bAnimationsPack, strFileName);
                            break;
                    }
                }
                catch (Exception ex)
                {
                    strGlobalExceptionMessage = ex.Message;

                    iinputFrameData = -1;
                }
            }
            return iinputFrameData;
        }
    }
}
