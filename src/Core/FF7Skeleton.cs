using System;
using System.IO;
using System.Collections.Generic;

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
        public static ModelType modelType = ModelType.K_NONE;

        // Animation constants for skeleton
        public const int K_FRAME_BONE_ROTATION = 0;
        public const int K_FRAME_ROOT_ROTATION = 1;
        public const int K_FRAME_ROOT_TRANSLATION = 2;

        // Global vars
        public static FieldSkeleton fSkeleton;
        public static FieldAnimation fAnimation;

        public static BattleSkeleton bSkeleton;
        public static BattleAnimationsPack bAnimationsPack;

        public static PModel fPModel;
        public static TMDModel mTMDModel;

        public static bool IsTMDModel;
        public static bool IsRSDResource;

        public static bool bLoaded;

        //
        // Global Skeleton/Model functions/procedures
        //
        public static int LoadSkeleton(string strFileName, bool loadGeometryQ, bool isLimitBreak,
                                       bool ignoreMissingPFiles, bool repairPolys,
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
                    // We load the Field Skeleton into memory.
                    switch (modelType)
                    {
                        case ModelType.K_HRC_SKELETON:
                            string strAnimationName = "";

                            // Field Skeleton (.hrc)
                            fSkeleton = new FieldSkeleton(strFileName, loadGeometryQ, ignoreMissingPFiles,
                                                          repairPolys, removeTextureCoords);

                            // We try to find some compatible Field Animation for the Field Skeleton.
                            // If there is no compatible field animation we have this var:   strGlobalFieldAnimationName = ""
                            iloadSkeletonResult = SearchFirstCompatibleFieldAnimationFileName(fSkeleton, 
                                                                                              Path.GetDirectoryName(strFileName), 
                                                                                              ref strAnimationName);

                            if (iloadSkeletonResult == 1)
                                fAnimation = new FieldAnimation(fSkeleton,
                                                                Path.GetDirectoryName(strFileName) + "\\" + strAnimationName, 
                                                                strAnimationName != "DUMMY.A");

                            break;

                        case ModelType.K_AA_SKELETON:
                            // Battle Skeleton (aa)
                            bSkeleton = new BattleSkeleton(strFileName, isLimitBreak, true, repairPolys,
                                                           removeTextureCoords);

                            // Normally we will have the ??DA file with the Animation Pack.
                            // Location Battle Models has NOT ??DA file.
                            // But editing models, it is possible we work without it. So, we will make something
                            // similiar as we did with Field Models, but we will check if ??DA file for the model exists.
                            bAnimationsPack = new BattleAnimationsPack(bSkeleton, strFileName);

                            break;

                        case ModelType.K_MAGIC_SKELETON:
                            // Magic Skeleton (.d)
                            bSkeleton = new BattleSkeleton(strFileName, true, repairPolys, removeTextureCoords);

                            // Normally we will have the *.A00 file with the Animation Pack.
                            // But editing models, it is possible we work without it. So, we will make something
                            // similiar as we did with Field Models, but we will check if *.A00 file for the model exists.
                            bAnimationsPack = new BattleAnimationsPack(bSkeleton, strFileName);

                            break;
                    }
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

        public static int LoadFieldSkeletonFromDB(string strFileName, 
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
                        case ModelType.K_HRC_SKELETON:
                            // We load the Field Skeleton into memory.

                            // Field Skeleton (.hrc)
                            fSkeleton = new FieldSkeleton(strFileName, loadGeometryQ, ignoreMissingPFiles,
                                                          repairPolys, removeTextureCoords);

                            iloadSkeletonResult = LoadAnimationFromDB(strAnimFileName);
                            break;

                        case ModelType.K_AA_SKELETON:
                            // Battle Skeleton (aa)
                            bSkeleton = new BattleSkeleton(strFileName, isLimitBreak, true, repairPolys,
                                                           removeTextureCoords);

                            // Normally we will have the ??DA file with the Animation Pack.
                            // Location Battle Models has NOT ??DA file.
                            // But editing models, it is possible we work without it. So, we will make something
                            // similiar as we did with Field Models, but we will check if ??DA file for the model exists.
                            bAnimationsPack = new BattleAnimationsPack(bSkeleton, strFileName);
                            break;

                        case ModelType.K_MAGIC_SKELETON:
                            // Magic Skeleton (.d)
                            bSkeleton = new BattleSkeleton(strFileName, true, repairPolys, removeTextureCoords);

                            // Normally we will have the ??DA file with the Animation Pack.
                            // Location Battle Models has NOT ??DA file.
                            // But editing models, it is possible we work without it. So, we will make something
                            // similiar as we did with Field Models, but we will check if ??DA file for the model exists.
                            bAnimationsPack = new BattleAnimationsPack(bSkeleton, strFileName);
                            break;
                    }
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

            try
            {
                switch (modelType)
                {
                    case ModelType.K_HRC_SKELETON:
                        DestroyFieldSkeleton(fSkeleton);
                        break;

                    case ModelType.K_AA_SKELETON:
                    case ModelType.K_MAGIC_SKELETON:
                        DestroyBattleSkeleton(bSkeleton);
                        break;
                }
            }
            catch
            {
                iDestroySkeletonResult = -1;
            }

            return iDestroySkeletonResult;
        }

        public static ModelType GetSkeletonType(string strFileName)
        {
            ModelType iSkeletonType = ModelType.K_NONE;
            string tmpFileName;

            if (strFileName.Length > 0)
            {
                switch (Path.GetExtension(strFileName).ToUpper())
                {
                    case ".HRC":
                        iSkeletonType = ModelType.K_HRC_SKELETON;
                        break;

                    case "":
                        tmpFileName = Path.GetFileNameWithoutExtension(strFileName).ToUpper();
                        if (tmpFileName.Length > 2)
                        {
                            if (tmpFileName[tmpFileName.Length - 1] == 'A' && tmpFileName[tmpFileName.Length - 2] == 'A')
                                iSkeletonType = ModelType.K_AA_SKELETON;
                        }
                        break;

                    case ".D":
                        iSkeletonType = ModelType.K_MAGIC_SKELETON;
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
            List<TEX> textures_pool = new List<TEX>();

            try
            {
                // Create fSkeleton with 1 bone
                fSkeleton = new FieldSkeleton()
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
                fAnimation = new FieldAnimation(fSkeleton,
                                                strRSDFolder + "\\" + strfAnimation,
                                                false);

                modelType = ModelType.K_HRC_SKELETON;
                IsRSDResource = true;
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
                    case ModelType.K_HRC_SKELETON:
                        WriteFieldAnimation(fAnimation, strFileName);

                        isaveAnimationResult = 1;
                        break;

                    case ModelType.K_AA_SKELETON:
                    case ModelType.K_MAGIC_SKELETON:
                        isaveAnimationResult = WriteBattleAnimationsPack(ref bAnimationsPack, strFileName);
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
                case ModelType.K_HRC_SKELETON:
                    if (fAnimation.nFrames == 1) iNumAnimFramesIsOneResult = true;
                    break;

                case ModelType.K_AA_SKELETON:
                case ModelType.K_MAGIC_SKELETON:
                    if (bAnimationsPack.SkeletonAnimations[index].numFramesShort == 1) iNumAnimFramesIsOneResult = true;
                    break;
            }

            return iNumAnimFramesIsOneResult;
        }

        public static int ReadFrameData(string strFileName, bool bMerge)
        {
            int iinputFrameData = 0;

            try
            {
                switch (modelType)
                {
                    case ModelType.K_HRC_SKELETON:
                        ReadFieldFrameData(fSkeleton, ref fAnimation, strFileName, bMerge);

                        iinputFrameData = 1;
                        break;

                    case ModelType.K_AA_SKELETON:
                    case ModelType.K_MAGIC_SKELETON:
                        //iinputFrameData = WriteBattleFrameDataPack(ref bAnimationsPack, strFileName);
                        break;
                }
            }
            catch (Exception ex)
            {
                strGlobalExceptionMessage = ex.Message;

                iinputFrameData = -1;
            }

            return iinputFrameData;
        }

        public static int WriteFrameData(string strFileName)
        {
            int ioutputFrameData = 0;

            try
            {
                switch (modelType)
                {
                    case ModelType.K_HRC_SKELETON:
                        WriteFieldFrameData(fSkeleton, fAnimation, strFileName);

                        ioutputFrameData = 1;
                        break;

                    case ModelType.K_AA_SKELETON:
                    case ModelType.K_MAGIC_SKELETON:
                        //ioutputFrameData = WriteBattleAnimationsPack(ref bAnimationsPack, strFileName);
                        break;
                }
            }
            catch (Exception ex)
            {
                strGlobalExceptionMessage = ex.Message;

                ioutputFrameData = -1;
            }

            return ioutputFrameData;
        }

        public static int ReadFrameDataSelective(string strFileName)
        {
            int iinputFrameData = 0;

            try
            {
                switch (modelType)
                {
                    case ModelType.K_HRC_SKELETON:
                        ReadFieldFrameDataSelective(fSkeleton, ref fAnimation, strFileName);

                        iinputFrameData = 1;
                        break;

                    case ModelType.K_AA_SKELETON:
                    case ModelType.K_MAGIC_SKELETON:
                        //iinputFrameData = WriteBattleFrameDataPack(ref bAnimationsPack, strFileName);
                        break;
                }
            }
            catch (Exception ex)
            {
                strGlobalExceptionMessage = ex.Message;

                iinputFrameData = -1;
            }

            return iinputFrameData;
        }



    }
}
