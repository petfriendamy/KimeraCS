using System;
using System.Collections.Generic;
using System.IO;

using OpenTK.Graphics.OpenGL.Compatibility;
using OpenTK.Mathematics;

namespace KimeraCS
{
    using static FF7Skeleton;
    using static FF7BattleAnimation;

    using static FF7TEXTexture;
    using static FF7PModel;

    using static Utils;

    public static class FF7BattleSkeleton
    {

        //
        // Battle Skeleton Structure
        //
        public struct BattleSkeleton
        {
            public string fileName;
            public int skeletonType;                       //  0 - Enemy Model, 1 - Battle Location, 2 - PC Battle Model?
            public int unk1;                               //  Always 1?
            public int unk2;                               //  Always 1?
            public int nBones;
            public int unk3;                               //  Always 0?
            public int nJoints;
            public int nTextures;
            public int nsSkeletonAnims;
            public int unk4;                               //  Num Skeleton Anims + 2?
            public int nWeapons;
            public int nsWeaponsAnims;
            public int unk5;                               //  Always 0?
            public int unk6;                             //  Global len?

            public List<BattleBone> bones;
            public List<TEX> textures;
            public List<PModel> wpModels;
            public uint[] TexIDS;
            public bool IsBattleLocation;
            public bool CanHaveLimitBreak;

            //  Constructor for the Battle Skeleton (battle.lgp files with ??AA filename format)
            public BattleSkeleton(string strFullFileName, bool isLimitBreak, bool loadGeometryQ,
                                  bool repairPolys, bool removeTextureCoords)
            {
                int pSuffix1, pSuffix2, pSuffix2End;
                string baseBattleSkeletonName;
                string weaponFileName;
                string strFileDirectoryName;
                int bi, ti, ji;

                byte[] fileBuffer;

                PModel tmpWPModel;
                TEX tmpTEX;

                BattleBone tmpbBone;

                fileName = Path.GetFileName(strFullFileName).ToUpper();
                strFileDirectoryName = Path.GetDirectoryName(strFullFileName);

                // Let's read Main Battle Skeleton part into memory.
                fileBuffer = File.ReadAllBytes(strFullFileName);

                textures = new List<TEX>();
                bones = new List<BattleBone>();
                wpModels = new List<PModel>();

                // Read memory fileBuffer
                using (var fileMemory = new MemoryStream(fileBuffer))
                {
                    using (var memReader = new BinaryReader(fileMemory))
                    {
                        skeletonType = memReader.ReadInt32();
                        unk1 = memReader.ReadInt32();
                        unk2 = memReader.ReadInt32();
                        nBones = memReader.ReadInt32();

                        unk3 = memReader.ReadInt32();
                        nJoints = memReader.ReadInt32();
                        nTextures = memReader.ReadInt32();
                        nsSkeletonAnims = memReader.ReadInt32();

                        unk4 = memReader.ReadInt32();
                        nWeapons = memReader.ReadInt32();
                        nsWeaponsAnims = memReader.ReadInt32();
                        unk5 = memReader.ReadInt32();
                        unk6 = memReader.ReadInt32();

                        CanHaveLimitBreak = isLimitBreak;
                        baseBattleSkeletonName = fileName.Substring(0, 2);
                        pSuffix1 = 'A';

                        if (nBones == 0)
                        {
                            // It's Battle Location
                            IsBattleLocation = true;

                            pSuffix1 = 'A';
                            pSuffix2 = 'M';

                            for (ji = 0; ji < nJoints; ji++)
                            {
                                if (pSuffix2 > 'Z')
                                {
                                    pSuffix2 = 'A';
                                    pSuffix1++;
                                }

                                if (loadGeometryQ)
                                {
                                    tmpbBone = new BattleBone() 
                                    {
                                        Models = new List<PModel>(),
                                    };

                                    LoadBattleLocationPiece(ref tmpbBone, nBones, 
                                                            strFileDirectoryName, 
                                                            baseBattleSkeletonName + Convert.ToChar(pSuffix1) +
                                                                                     Convert.ToChar(pSuffix2),
                                                            repairPolys, removeTextureCoords);
                                    nBones++;
                                    bones.Add(tmpbBone);

                                    pSuffix2++;
                                }
                            }
                        }
                        else
                        {
                            //  It's a character battle model
                            IsBattleLocation = false;

                            // Read Battle Bones files
                            pSuffix2 = 'M';

                            for (bi = 0; bi < nBones; bi++)
                            {
                                if (pSuffix2 > 'Z')
                                {
                                    pSuffix1++;
                                    pSuffix2 = 'A';
                                }

                                bones.Add(new BattleBone(memReader, strFileDirectoryName, 
                                                         baseBattleSkeletonName + Convert.ToChar(pSuffix1) +
                                                                                  Convert.ToChar(pSuffix2), 
                                                         loadGeometryQ, repairPolys, removeTextureCoords));

                                pSuffix2++;                                
                            }

                            //  Read Battle Weapon files
                            pSuffix2End = 'K' + nWeapons;
;
                            if (nWeapons > 0)
                            {
                                for (pSuffix2 = 'K'; pSuffix2 < pSuffix2End; pSuffix2++)
                                {
                                    weaponFileName = baseBattleSkeletonName + 'C' + Convert.ToChar(pSuffix2);

                                    if (File.Exists(strFileDirectoryName + "\\" + weaponFileName))
                                    {
                                        if (loadGeometryQ)
                                        {
                                            tmpWPModel = new PModel();

                                            LoadPModel(ref tmpWPModel, strFileDirectoryName, weaponFileName, true,
                                                       repairPolys, removeTextureCoords);
                                            wpModels.Add(tmpWPModel);
                                        }

                                        //  Debug.Print "Loaded weapon model " + weaponFileName
                                    }
                                    else
                                    {
                                        tmpWPModel = new PModel();
                                        wpModels.Add(tmpWPModel);
                                    }
                                }
                            }
                        }

                        //  Read Battle Textures files
                        TexIDS = new uint[nTextures];

                        if (loadGeometryQ)
                        {
                            textures = new List<TEX>();

                            ti = 0;

                            pSuffix2End = 'C' + nTextures;

                            for (pSuffix2 = 'C'; pSuffix2 < pSuffix2End; pSuffix2++)
                            {
                                tmpTEX = new TEX() 
                                {
                                    TEXfileName = baseBattleSkeletonName.ToUpper() + "A" + Convert.ToChar(pSuffix2),
                                };
                                
                                if (ReadTEXTexture(ref tmpTEX, strFileDirectoryName + "\\" + tmpTEX.TEXfileName) == 0)
                                {
                                    LoadTEXTexture(ref tmpTEX);
                                    LoadBitmapFromTEXTexture(ref tmpTEX);
                                }

                                TexIDS[ti] = tmpTEX.texID;

                                textures.Add(tmpTEX);

                                ti++;
                            }
                        }
                    }
                }
            }


            //  Constructor for the Magic Skeleton (magic.lgp files with .D extension)
            public BattleSkeleton(string strFullFileName, bool loadGeometryQ, bool repairPolys,
                                  bool removeTextureCoords)
            {
                int bi, ti;
                string baseMagicSkeletonName, strFileDirectoryName;
                string pSuffix, tSuffix;
                byte[] fileBuffer;
                TEX tmpTEX;

                fileName = Path.GetFileName(strFullFileName).ToUpper();
                strFileDirectoryName = Path.GetDirectoryName(strFullFileName);

                // Let's read Main Battle Skeleton part into memory.
                fileBuffer = File.ReadAllBytes(strFullFileName);

                IsBattleLocation = false;
                CanHaveLimitBreak = false;

                bones = new List<BattleBone>();
                textures = new List<TEX>();
                wpModels = new List<PModel>();      // This is not used but we need it to initialize in the struct constructor

                // Read memory fileBuffer
                baseMagicSkeletonName = fileName.Substring(0, fileName.IndexOf('.'));

                using (var fileMemory = new MemoryStream(fileBuffer))
                {
                    using (var memReader = new BinaryReader(fileMemory))
                    {
                        skeletonType = memReader.ReadInt32();
                        unk1 = memReader.ReadInt32();
                        unk2 = memReader.ReadInt32();
                        nBones = memReader.ReadInt32();

                        unk3 = memReader.ReadInt32();
                        nJoints = memReader.ReadInt32();
                        nTextures = memReader.ReadInt32();
                        nsSkeletonAnims = memReader.ReadInt32();

                        unk4 = memReader.ReadInt32();
                        nWeapons = memReader.ReadInt32();
                        nsWeaponsAnims = memReader.ReadInt32();
                        unk5 = memReader.ReadInt32();
                        unk6 = memReader.ReadInt32();

                        //  Read Magic Bones files (P?? model files)
                        for (bi = 0; bi < nBones; bi++)
                        {
                            pSuffix = ".P" + bi.ToString("00");

                            // Have in mind that not all the models exists.
                            bones.Add(new BattleBone(memReader, strFileDirectoryName, baseMagicSkeletonName + pSuffix,
                                                     loadGeometryQ, repairPolys, removeTextureCoords));
                        }

                        //  Read Magic Texture files (T?? tex files)
                        TexIDS = new uint[nTextures];

                        if (loadGeometryQ)
                        {
                            for (ti = 0; ti < nTextures; ti++)
                            {
                                tSuffix = ".T" + ti.ToString("00");

                                tmpTEX = new TEX()
                                {
                                    TEXfileName = baseMagicSkeletonName.ToUpper() + tSuffix,
                                };

                                if (ReadTEXTexture(ref tmpTEX, strFileDirectoryName + "\\" + tmpTEX.TEXfileName) == 0)
                                {
                                    LoadTEXTexture(ref tmpTEX);
                                    LoadBitmapFromTEXTexture(ref tmpTEX);
                                }

                                TexIDS[ti] = tmpTEX.texID;

                                textures.Add(tmpTEX);
                            }
                        }
                    }
                }
            }

            //checks if the internal polys need to be repaired
            public PModel? CheckPolys()
            {
                foreach (var bone in bones)
                {
                    foreach (var p in bone.Models)
                    {
                        var curr = p;
                        int result = FF7PModel.CheckPolys(ref curr);
                        if (result >= 0) return curr;
                    }
                }
                return null;
            }

            //attempts to repair internal polys
            public void RepairPolys()
            {
                foreach (var bone in bones)
                {
                    foreach (var p in bone.Models)
                    {
                        var curr = p;
                        FF7PModel.RepairPolys(ref curr);
                    }
                }
            }
        }



        //
        // Battle Skeleton Bone Structure
        //
        public struct BattleBone
        {
            public int parentBone;
            public float len;
            public int hasModel;
            public List<PModel> Models;
            //  -------------Extra Atributes----------------
            public int nModels;
            public float resizeX;
            public float resizeY;
            public float resizeZ;

            public BattleBone(BinaryReader memReader, string strDirectoryName, string modelName,
                              bool loadGeometryQ, bool repairPolys, bool removeTextureCoords)
            {
                PModel tmpbPModel;

                parentBone = memReader.ReadInt32();
                len = memReader.ReadSingle();
                hasModel = memReader.ReadInt32();
                nModels = 0;

                Models = new List<PModel>();

                if (hasModel != 0)
                {
                    if (loadGeometryQ)
                    {
                        nModels = 1;

                        tmpbPModel = new PModel();
                        LoadPModel(ref tmpbPModel, strDirectoryName, modelName, true, repairPolys, removeTextureCoords);
                        Models.Add(tmpbPModel);
                    }
                }

                resizeX = 1;
                resizeY = 1;
                resizeZ = 1;
            }
        }

        public static void LoadBattleLocationPiece(ref BattleBone bBone, int boneIndex, string strDirectoryName,
                                                   string modelName, bool repairPolys, bool removeTextureCoords)
        {
            PModel bLocBone;
            
            bBone.parentBone = boneIndex;
            bBone.hasModel = 1;
            bBone.nModels = 1;

            bLocBone = new PModel();
            LoadPModel(ref bLocBone, strDirectoryName, modelName, true, repairPolys, removeTextureCoords);
            bBone.Models.Add(bLocBone);

            bBone.len = ComputeDiameter(bLocBone.BoundingBox) / 2;
            bBone.resizeX = 1;
            bBone.resizeY = 1;
            bBone.resizeZ = 1;
        }

        //
        // Battle Skeleton functions
        //
        public static void ComputeBattleBoneBoundingBox(BattleBone bBone, ref Point3D p_min, ref Point3D p_max)
        {
            int mi;
            double[] MV_matrix = new double[16];

            Point3D p_min_aux;
            Point3D p_max_aux;
            Point3D p_min_aux_trans = new Point3D();
            Point3D p_max_aux_trans = new Point3D();

            // Build base transform using pure math
            Matrix4 baseMatrix = Matrix4.CreateScale(bBone.resizeX, bBone.resizeY, bBone.resizeZ);

            if (bBone.hasModel == 1)
            {
                p_max.x = (float)-INFINITY_SINGLE;
                p_max.y = (float)-INFINITY_SINGLE;
                p_max.z = (float)-INFINITY_SINGLE;

                p_min.x = (float)INFINITY_SINGLE;
                p_min.y = (float)INFINITY_SINGLE;
                p_min.z = (float)INFINITY_SINGLE;

                for (mi = 0; mi < bBone.nModels; mi++)
                {
                    // Build model transform
                    Matrix4 modelMatrix = baseMatrix;
                    modelMatrix *= Matrix4.CreateTranslation(bBone.Models[mi].repositionX, bBone.Models[mi].repositionY, bBone.Models[mi].repositionZ);
                    modelMatrix *= Matrix4.CreateRotationY(MathHelper.DegreesToRadians((float)bBone.Models[mi].rotateBeta));
                    modelMatrix *= Matrix4.CreateRotationX(MathHelper.DegreesToRadians((float)bBone.Models[mi].rotateAlpha));
                    modelMatrix *= Matrix4.CreateRotationZ(MathHelper.DegreesToRadians((float)bBone.Models[mi].rotateGamma));
                    modelMatrix *= Matrix4.CreateScale(bBone.resizeX, bBone.resizeY, bBone.resizeZ);

                    p_min_aux.x = bBone.Models[mi].BoundingBox.min_x;
                    p_min_aux.y = bBone.Models[mi].BoundingBox.min_y;
                    p_min_aux.z = bBone.Models[mi].BoundingBox.min_z;

                    p_max_aux.x = bBone.Models[mi].BoundingBox.max_x;
                    p_max_aux.y = bBone.Models[mi].BoundingBox.max_y;
                    p_max_aux.z = bBone.Models[mi].BoundingBox.max_z;

                    MV_matrix = Matrix4ToDoubleArray(modelMatrix);

                    ComputeTransformedBoxBoundingBox(MV_matrix, ref p_min_aux, ref p_max_aux, ref p_min_aux_trans, ref p_max_aux_trans);

                    if (p_max.x < p_max_aux_trans.x) p_max.x = p_max_aux_trans.x;
                    if (p_max.y < p_max_aux_trans.y) p_max.y = p_max_aux_trans.y;
                    if (p_max.z < p_max_aux_trans.z) p_max.z = p_max_aux_trans.z;

                    if (p_min.x > p_min_aux_trans.x) p_min.x = p_min_aux_trans.x;
                    if (p_min.y > p_min_aux_trans.y) p_min.y = p_min_aux_trans.y;
                    if (p_min.z > p_min_aux_trans.z) p_min.z = p_min_aux_trans.z;
                }
            }
            else
            {
                p_max.x = 0;
                p_max.y = 0;
                p_max.z = 0;

                p_min.x = 0;
                p_min.y = 0;
                p_min.z = 0;
            }
        }

        public static void ComputeBattleBoundingBox(BattleSkeleton bSkeleton, BattleFrame bFrame, ref Point3D p_min, ref Point3D p_max)
        {
            double[] rot_mat = new double[16];
            double[] MV_matrix = new double[16];

            Point3D p_max_bone = new Point3D();
            Point3D p_min_bone = new Point3D();
            Point3D p_max_bone_trans = new Point3D();
            Point3D p_min_bone_trans = new Point3D();

            int[] joint_stack = new int[bSkeleton.nBones * 4];
            Matrix4[] matrixStack = new Matrix4[bSkeleton.nBones + 2];
            int matrixStackPtr = 0;
            int jsp, bi, iframeCnt;

            jsp = 0;
            joint_stack[jsp] = -1;

            p_max.x = -(float)INFINITY_SINGLE;
            p_max.y = -(float)INFINITY_SINGLE;
            p_max.z = -(float)INFINITY_SINGLE;
            p_min.x = (float)INFINITY_SINGLE;
            p_min.y = (float)INFINITY_SINGLE;
            p_min.z = (float)INFINITY_SINGLE;

            // Build initial transform using pure math
            Matrix4 currentMatrix = Matrix4.Identity;
            currentMatrix *= Matrix4.CreateTranslation((float)bFrame.startX, (float)bFrame.startY, (float)bFrame.startZ);

            BuildRotationMatrixWithQuaternions(bFrame.bones[0].alpha, bFrame.bones[0].beta, bFrame.bones[0].gamma, ref rot_mat);
            Matrix4 rotMatrix = DoubleArrayToMatrix4(rot_mat);
            currentMatrix *= rotMatrix;

            matrixStack[matrixStackPtr++] = currentMatrix;

            for (bi = 0; bi < bSkeleton.nBones; bi++)
            {
                while (!(bSkeleton.bones[bi].parentBone == joint_stack[jsp]) && jsp > 0)
                {
                    matrixStackPtr--;
                    currentMatrix = matrixStack[matrixStackPtr];
                    jsp--;
                }
                matrixStack[matrixStackPtr++] = currentMatrix;

                if (bSkeleton.nBones > 1) iframeCnt = 1;
                else iframeCnt = 0;
                BuildRotationMatrixWithQuaternions(bFrame.bones[bi + iframeCnt].alpha,
                                                   bFrame.bones[bi + iframeCnt].beta,
                                                   bFrame.bones[bi + iframeCnt].gamma, ref rot_mat);
                rotMatrix = DoubleArrayToMatrix4(rot_mat);
                currentMatrix *= rotMatrix;

                ComputeBattleBoneBoundingBox(bSkeleton.bones[bi], ref p_min_bone, ref p_max_bone);

                MV_matrix = Matrix4ToDoubleArray(currentMatrix);

                ComputeTransformedBoxBoundingBox(MV_matrix, ref p_min_bone, ref p_max_bone, ref p_min_bone_trans, ref p_max_bone_trans);

                if (p_max.x < p_max_bone_trans.x) p_max.x = p_max_bone_trans.x;
                if (p_max.y < p_max_bone_trans.y) p_max.y = p_max_bone_trans.y;
                if (p_max.z < p_max_bone_trans.z) p_max.z = p_max_bone_trans.z;

                if (p_min.x > p_min_bone_trans.x) p_min.x = p_min_bone_trans.x;
                if (p_min.y > p_min_bone_trans.y) p_min.y = p_min_bone_trans.y;
                if (p_min.z > p_min_bone_trans.z) p_min.z = p_min_bone_trans.z;

                currentMatrix *= Matrix4.CreateTranslation(0, 0, bSkeleton.bones[bi].len);
                jsp++;
                joint_stack[jsp] = bi;
            }
        }

        public static float ComputeBattleDiameter(BattleSkeleton bSkeleton)
        {
            int bi;

            float maxPath, currentPath;
            int jsp;
            int[] joint_stack = new int[bSkeleton.nBones + 1];

            maxPath = 0;
            currentPath = 0;
            jsp = 0;

            if (bSkeleton.IsBattleLocation)
            {
                for (bi = 0; bi < bSkeleton.nBones; bi++)
                    if (bSkeleton.bones[bi].len > maxPath) maxPath = bSkeleton.bones[bi].len;
            }
            else
            {
                if (bSkeleton.nBones == 1 && bSkeleton.bones[0].len <= 0)
                    maxPath = bSkeleton.bones[0].Models[0].diameter;
                else
                {
                    joint_stack[jsp] = -1;

                    for (bi = 0; bi < bSkeleton.nBones; bi++)
                    {
                        while (!(bSkeleton.bones[bi].parentBone == joint_stack[jsp]) && jsp > 0)
                        {
                            currentPath += bSkeleton.bones[joint_stack[jsp]].len;
                            jsp--;
                        }

                        currentPath -= bSkeleton.bones[bi].len;
                        if (currentPath > maxPath) maxPath = currentPath;
                        jsp++;
                        joint_stack[jsp] = bi;
                    }
                }
            }

            return maxPath;
        }

        public static string GetBattleModelTextureFilename(BattleSkeleton bSkeleton, int nTex)
        {
            return bSkeleton.fileName.Substring(0, 2) + 'A' + Convert.ToChar('C' + nTex);
        }

        public static void AddBattleBoneModel(ref BattleBone bBone, ref PModel Model)
        {           
            if (bBone.nModels > 0)
            {
                if (modelType == K_AA_SKELETON)
                {
                    Model.fileName = bBone.Models[0].fileName + (bBone.nModels - 1).ToString();
                }
                else
                {
                    Model.fileName = bBone.Models[0].fileName.Substring(0, bBone.Models[0].fileName.IndexOf('.')) + ".P" + (bBone.nModels - 1).ToString();
                }
            }

            bBone.Models.Add(Model);
            bBone.nModels++;
            bBone.hasModel = 1;
        }

        public static void RemoveBattleBoneModel(ref BattleBone bBone, ref int b_index)
        {
            if (b_index < bBone.nModels)
            {
                bBone.nModels--;
                bBone.Models.RemoveAt(b_index);

                if (bBone.nModels == 0) bBone.hasModel = 0;
            }
        }



        //  ---------------------------------------------------------------------------------------------------
        //  ============================================= SAVING ==============================================
        //  ---------------------------------------------------------------------------------------------------
        public static void MergeBattleBoneModels(ref BattleBone bBone)
        {
            int mi;
            PModel tmpModel;

            for (mi = 1; mi < bBone.nModels; mi++)

            {
                tmpModel = bBone.Models[0];
                MergePModels(ref tmpModel, bBone.Models[mi]);
                bBone.Models[0] = tmpModel;
            }
        }

        public static void ApplyBattleBoneChanges(ref BattleBone bBone) {
            int mi;
            PModel tmpModel;

            for (mi = 0; mi < bBone.nModels; mi++)
            {
                if (bBone.hasModel == 1)
                {
                    if (GL.IsEnabled(EnableCap.Lighting)) 
                    {
                        tmpModel = bBone.Models[mi];
                        ApplyCurrentVColors(ref tmpModel);
                        bBone.Models[mi] = tmpModel;
                    }

                    GL.MatrixMode(MatrixMode.Modelview);
                    GL.PushMatrix();

                    SetCameraModelViewQuat(bBone.Models[mi].repositionX, bBone.Models[mi].repositionY, bBone.Models[mi].repositionZ,
                                           bBone.Models[mi].rotationQuaternion,
                                           bBone.Models[mi].resizeX, bBone.Models[mi].resizeY, bBone.Models[mi].resizeZ);

                    GL.Scaled(bBone.resizeX, bBone.resizeY, bBone.resizeZ);

                    tmpModel = bBone.Models[mi];
                    ApplyPChanges(ref tmpModel, false);
                    bBone.Models[mi] = tmpModel;

                    GL.MatrixMode(MatrixMode.Modelview);
                    GL.PopMatrix();
                }
            }

            MergeBattleBoneModels(ref bBone);

            if (bBone.nModels > 1) while (bBone.nModels > 1) bBone.Models.RemoveAt(bBone.Models.Count - 1);
            bBone.nModels = 1;
        }

        public static void ApplyBattleWeaponChanges(ref PModel wpModel)
        {
            if (GL.IsEnabled(EnableCap.Lighting)) ApplyCurrentVColors(ref wpModel);

            GL.MatrixMode(MatrixMode.Modelview);
            GL.PushMatrix();

            SetCameraModelView(wpModel.repositionX, wpModel.repositionY, wpModel.repositionZ,
                               wpModel.rotateAlpha, wpModel.rotateBeta, wpModel.rotateGamma,
                               wpModel.resizeX, wpModel.resizeY, wpModel.resizeZ);

            GL.Scaled(wpModel.resizeX, wpModel.resizeY, wpModel.resizeZ);

            ApplyPChanges(ref wpModel, true);

            GL.MatrixMode(MatrixMode.Modelview);
            GL.PopMatrix();
        }

        public static void ApplyBattleChanges(ref BattleSkeleton bSkeleton, BattleFrame bFrame, BattleFrame bwpFrame)
        {
            int bi, wi, jsp;
            int[] joint_stack = new int[bSkeleton.nBones + 1];
            double[] rot_mat = new double[16];
            BattleBone tmpbBone;

            jsp = 0;
            joint_stack[0] = -1;

            GL.MatrixMode(MatrixMode.Modelview);
            GL.PushMatrix();
            // glLoadIdentity();
            GL.Translated(bFrame.startX, bFrame.startY, bFrame.startZ);

            BuildRotationMatrixWithQuaternions(bFrame.bones[0].alpha, bFrame.bones[0].beta, bFrame.bones[0].gamma, ref rot_mat);
            GL.MultMatrixd(rot_mat);

            for (bi = 0; bi < bSkeleton.nBones; bi++)
            {
                while (!(bSkeleton.bones[bi].parentBone == joint_stack[jsp]) && jsp > 0)
                {
                    GL.PopMatrix();
                    jsp--;
                }
                GL.PushMatrix();

                //GL.Rotated(bFrame.bones[bi + 1].beta, 0, 1, 0);
                //GL.Rotated(bFrame.bones[bi + 1].alpha, 1, 0, 0);
                //GL.Rotated(bFrame.bones[bi + 1].gamma, 0, 0, 1);

                BuildRotationMatrixWithQuaternions(bFrame.bones[bi + (bSkeleton.nBones > 1 ? 1 : 0)].alpha,
                                                   bFrame.bones[bi + (bSkeleton.nBones > 1 ? 1 : 0)].beta,
                                                   bFrame.bones[bi + (bSkeleton.nBones > 1 ? 1 : 0)].gamma,
                                                   ref rot_mat);
                GL.MultMatrixd(rot_mat);

                if (bSkeleton.bones[bi].hasModel == 1)
                {
                    tmpbBone = bSkeleton.bones[bi];
                    ApplyBattleBoneChanges(ref tmpbBone);
                    bSkeleton.bones[bi] = tmpbBone;
                }

                GL.Translated(0, 0, bSkeleton.bones[bi].len);

                jsp++;
                joint_stack[jsp] = bi;
            }

            while (jsp > 0)
            {
                GL.PopMatrix();
                jsp--;
            }
            GL.PopMatrix();

            if (bSkeleton.wpModels.Count > 0)
            {
                PModel tmpwpModel;

                GL.MatrixMode(MatrixMode.Modelview);
                GL.PushMatrix();
                //glLoadIdentity();
                GL.Translated(bwpFrame.startX, bwpFrame.startY, bwpFrame.startZ);
                GL.MultMatrixd(rot_mat);

                BuildRotationMatrixWithQuaternions(bwpFrame.bones[0].alpha, bwpFrame.bones[0].beta, bwpFrame.bones[0].gamma, ref rot_mat);

                for (wi = 0; wi < bSkeleton.nWeapons; wi++)
                {
                    if (bSkeleton.wpModels[wi].Polys != null)
                    {
                        tmpwpModel = bSkeleton.wpModels[wi];
                        ApplyBattleWeaponChanges(ref tmpwpModel);
                        bSkeleton.wpModels[wi] = tmpwpModel;
                    }
                }
                GL.PopMatrix();
            }
        }

        public static void CreateDListsFromBattleSkeletonBone(ref BattleBone bBone)
        {
            int mi;
            PModel tmpModel;

            for (mi = 0; mi < bBone.nModels; mi++)
            {
                tmpModel = bBone.Models[mi];
                CreateDListsFromPModel(ref tmpModel);
                bBone.Models[mi] = tmpModel;
            }
        }

        public static void CreateDListsFromBattleSkeleton(ref BattleSkeleton bSkeleton)
        {
            int bi;
            BattleBone tmpbBone;

            for (bi = 0; bi < bSkeleton.bones.Count; bi++)
            {
                tmpbBone = bSkeleton.bones[bi];
                CreateDListsFromBattleSkeletonBone(ref tmpbBone);
                bSkeleton.bones[bi] = tmpbBone;
            }
        }

        public static void WriteBattleBone(ref BattleBone bBone, string strModelFileName)
        {
            PModel tmpModel;

            if (bBone.hasModel == 1)
            {
                tmpModel = bBone.Models[0];
                WriteGlobalPModel(ref tmpModel, strModelFileName);
                bBone.Models[0] = tmpModel;
            }

            bBone.resizeX = 1;
            bBone.resizeY = 1;
            bBone.resizeZ = 1;
        }

        public static void WriteBattleSkeleton(ref BattleSkeleton bSkeleton, string strFileName)
        {
            int iBoneIdx, iWeaponIdx, iTextureIdx;
            int pSuffix1, pSuffix2;
            string strBaseFileName, strFullDirectoryName;
            byte[] fileBuffer = new byte[(13 * 4) + (bSkeleton.nBones * 12)];
            BattleBone tmpbBone;
            PModel tmpModel;

            strBaseFileName = Path.GetFileNameWithoutExtension(strFileName).Substring(0, 2);
            strFullDirectoryName = Path.GetDirectoryName(strFileName);

            // Writer Main Battle file (AA)
            using (MemoryStream fileMemory = new MemoryStream(fileBuffer))
            {
                using (var memWriter = new BinaryWriter(fileMemory))
                {
                    memWriter.Write(bSkeleton.skeletonType);
                    memWriter.Write(bSkeleton.unk1);
                    memWriter.Write(bSkeleton.unk2);

                    if (bSkeleton.IsBattleLocation)
                        memWriter.Write((int)0);
                    else
                        memWriter.Write(bSkeleton.nBones);

                    memWriter.Write(bSkeleton.unk3);
                    memWriter.Write(bSkeleton.nJoints);
                    memWriter.Write(bSkeleton.nTextures);
                    memWriter.Write(bSkeleton.nsSkeletonAnims);

                    memWriter.Write(bSkeleton.unk4);
                    memWriter.Write(bSkeleton.nWeapons);
                    memWriter.Write(bSkeleton.nsWeaponsAnims);
                    memWriter.Write(bSkeleton.unk5);
                    memWriter.Write(bSkeleton.unk6);

                    for (iBoneIdx = 0; iBoneIdx < bSkeleton.nBones; iBoneIdx++)
                    {
                        memWriter.Write(bSkeleton.bones[iBoneIdx].parentBone);
                        memWriter.Write(bSkeleton.bones[iBoneIdx].len);
                        memWriter.Write(bSkeleton.bones[iBoneIdx].hasModel);
                    }
                }
            }
            File.WriteAllBytes(strFullDirectoryName + "\\" + strBaseFileName + "AA", fileBuffer);


            // Write Battle Bones files (AM->CJ)
            pSuffix1 = 'A';
            pSuffix2 = 'M';

            for (iBoneIdx = 0; iBoneIdx < bSkeleton.nBones; iBoneIdx++)
            {
                if (pSuffix2 > 'Z')
                {
                    pSuffix1++;
                    pSuffix2 = 'A';
                }

                tmpbBone = bSkeleton.bones[iBoneIdx];
                WriteBattleBone(ref tmpbBone, strFullDirectoryName + "\\" + 
                                              strBaseFileName + Convert.ToChar(pSuffix1) + 
                                              Convert.ToChar(pSuffix2));
                bSkeleton.bones[iBoneIdx] = tmpbBone;

                pSuffix2++;
            }


            // Write Battle Weapon files (CK->CZ)
            pSuffix1 = 'C';
            pSuffix2 = 'K';

            for (iWeaponIdx = 0; iWeaponIdx < bSkeleton.nWeapons; iWeaponIdx++)
            {
                if (bSkeleton.wpModels[iWeaponIdx].Polys != null)
                {
                    tmpModel = bSkeleton.wpModels[iWeaponIdx];
                    WriteGlobalPModel(ref tmpModel, strFullDirectoryName + "\\" + 
                                                    strBaseFileName + Convert.ToChar(pSuffix1) +
                                                    Convert.ToChar(pSuffix2 + iWeaponIdx));
                    bSkeleton.wpModels[iWeaponIdx] = tmpModel;
                }
            }


            // Write Battle Texture files (AC->AL)
            pSuffix1 = 'A';
            pSuffix2 = 'C';

            for (iTextureIdx = 0; iTextureIdx < bSkeleton.nTextures; iTextureIdx++)
            {
                WriteTEXTexture(bSkeleton.textures[iTextureIdx], strFullDirectoryName + "\\" + 
                                                                 strBaseFileName + 
                                                                 Convert.ToChar(pSuffix1) + 
                                                                 Convert.ToChar(pSuffix2 + iTextureIdx));
            }
        }

        public static void WriteMagicSkeleton(ref BattleSkeleton bSkeleton, string strFileName)
        {
            int iBoneIdx, iTextureIdx;
            string pSuffix, tSuffix, strBaseFileName, strFullDirectoryName;
            byte[] fileBuffer = new byte[(13 * 4) + (bSkeleton.nBones * 12)];
            BattleBone tmpbBone;

            strBaseFileName = Path.GetFileNameWithoutExtension(strFileName);
            strFullDirectoryName = Path.GetDirectoryName(strFileName);

            // Writer Main Magic file (.D)
            using (MemoryStream fileMemory = new MemoryStream(fileBuffer))
            {
                using (var memWriter = new BinaryWriter(fileMemory))
                {
                    memWriter.Write(bSkeleton.skeletonType);
                    memWriter.Write(bSkeleton.unk1);
                    memWriter.Write(bSkeleton.unk2);
                    memWriter.Write(bSkeleton.nBones);

                    memWriter.Write(bSkeleton.unk3);
                    memWriter.Write(bSkeleton.nJoints);
                    memWriter.Write(bSkeleton.nTextures);
                    memWriter.Write(bSkeleton.nsSkeletonAnims);

                    memWriter.Write(bSkeleton.unk4);
                    memWriter.Write(bSkeleton.nWeapons);
                    memWriter.Write(bSkeleton.nsWeaponsAnims);
                    memWriter.Write(bSkeleton.unk5);
                    memWriter.Write(bSkeleton.unk6);

                    for (iBoneIdx = 0; iBoneIdx < bSkeleton.nBones; iBoneIdx++)
                    {
                        memWriter.Write(bSkeleton.bones[iBoneIdx].parentBone);
                        memWriter.Write(bSkeleton.bones[iBoneIdx].len);
                        memWriter.Write(bSkeleton.bones[iBoneIdx].hasModel);
                    }
                }
            }
            File.WriteAllBytes(strFullDirectoryName + "\\" + strBaseFileName + ".D", fileBuffer);


            // Write Battle Bones files (.P??)
            for (iBoneIdx = 0; iBoneIdx < bSkeleton.nBones; iBoneIdx++)
            {
                pSuffix = ".P" + iBoneIdx.ToString("00");

                tmpbBone = bSkeleton.bones[iBoneIdx];
                WriteBattleBone(ref tmpbBone, strFullDirectoryName + "\\" + strBaseFileName + pSuffix);
                bSkeleton.bones[iBoneIdx] = tmpbBone;
            }


            // Write Battle Texture files (.T??)
            for (iTextureIdx = 0; iTextureIdx < bSkeleton.nTextures; iTextureIdx++)
            {
                tSuffix = ".T" + iTextureIdx.ToString("00");
                WriteTEXTexture(bSkeleton.textures[iTextureIdx], strFullDirectoryName + "\\" + 
                                                                 strBaseFileName + tSuffix);
            }
        }



        //  ---------------------------------------------------------------------------------------------------
        //  ============================================= DESTROY =============================================
        //  ---------------------------------------------------------------------------------------------------
        public static void DestroyBattleBoneResources(ref BattleBone bBone)
        {
            int iResourceIdx;
            PModel tmpModel;

            for (iResourceIdx = 0; iResourceIdx < bBone.nModels; iResourceIdx++)
            {
                tmpModel = bBone.Models[iResourceIdx];
                DestroyPModelResources(ref tmpModel);
                bBone.Models[iResourceIdx] = tmpModel;
            }

            if (bBone.Models != null) bBone.Models.Clear();
        }

        public static void DestroyBattleSkeleton(BattleSkeleton bSkeleton)
        {
            int iBoneIdx, iWeaponIdx, iTextureIdxbi;
            BattleBone tmpbBone;
            PModel tmpModel;

            if (bSkeleton.nBones > 0)
            {
                // Free skeleton models
                for (iBoneIdx = 0; iBoneIdx < bSkeleton.nBones; iBoneIdx++)
                {
                    tmpbBone = bSkeleton.bones[iBoneIdx];
                    DestroyBattleBoneResources(ref tmpbBone);
                    bSkeleton.bones[iBoneIdx] = tmpbBone;
                }

                if (bSkeleton.bones != null) bSkeleton.bones.Clear();

                // Free weapon models
                for (iWeaponIdx = 0; iWeaponIdx < bSkeleton.wpModels.Count; iWeaponIdx++)
                {
                    if (bSkeleton.wpModels[iWeaponIdx].Polys != null)
                    {
                        tmpModel = bSkeleton.wpModels[iWeaponIdx];
                        DestroyPModelResources(ref tmpModel);
                        bSkeleton.wpModels[iWeaponIdx] = tmpModel;
                    }
                }

                if (bSkeleton.wpModels != null) bSkeleton.wpModels.Clear();

                // Free textures
                int[] lsttexID = new int[1];
                for (iTextureIdxbi = 0; iTextureIdxbi < bSkeleton.textures.Count; iTextureIdxbi++)
                {
                    lsttexID[0] = (int)bSkeleton.textures[iTextureIdxbi].texID;
                    GL.DeleteTextures(1, lsttexID);
                    bSkeleton.textures[iTextureIdxbi].bitmap?.Dispose();
                }

                if (bSkeleton.textures != null) bSkeleton.textures.Clear();
            }
        }



        //  ---------------------------------------------------------------------------------------------------
        //  ========================================== COPY SKELETON ==========================================
        //  ---------------------------------------------------------------------------------------------------
        public static BattleBone CopybBone(BattleBone bBoneIn)
        {
            BattleBone bBoneOut;

            bBoneOut = new BattleBone()
            {
                hasModel = bBoneIn.hasModel,
                len = bBoneIn.len,
                nModels = bBoneIn.nModels,
                parentBone = bBoneIn.parentBone,
                resizeX = bBoneIn.resizeX,
                resizeY = bBoneIn.resizeY,
                resizeZ = bBoneIn.resizeZ,

                Models = new List<PModel>(),
            };

            foreach (PModel bModel in bBoneIn.Models) bBoneOut.Models.Add(CopyPModel(bModel));

            return bBoneOut;
        }

        public static BattleSkeleton CopybSkeleton(BattleSkeleton bSkeletonIn)
        {
            BattleSkeleton bSkeletonOut;

            bSkeletonOut = new BattleSkeleton()
            {
                fileName = bSkeletonIn.fileName,
                IsBattleLocation = bSkeletonIn.IsBattleLocation,
                CanHaveLimitBreak = bSkeletonIn.CanHaveLimitBreak,
                nBones = bSkeletonIn.nBones,
                nJoints = bSkeletonIn.nJoints,
                nsSkeletonAnims = bSkeletonIn.nsSkeletonAnims,
                nsWeaponsAnims = bSkeletonIn.nsWeaponsAnims,
                nTextures = bSkeletonIn.nTextures,
                nWeapons = bSkeletonIn.nWeapons,
                skeletonType = bSkeletonIn.skeletonType,
                TexIDS = bSkeletonIn.TexIDS,
                unk1 = bSkeletonIn.unk1,
                unk2 = bSkeletonIn.unk2,
                unk3 = bSkeletonIn.unk3,
                unk4 = bSkeletonIn.unk4,
                unk5 = bSkeletonIn.unk5,
                unk6 = bSkeletonIn.unk6,

                bones = new List<BattleBone>(),
                wpModels = new List<PModel>(),
                textures = new List<TEX>(),
            };

            foreach (BattleBone itmbBone in bSkeletonIn.bones) bSkeletonOut.bones.Add(CopybBone(itmbBone));          
            foreach (PModel itmbwpModel in bSkeletonIn.wpModels) bSkeletonOut.wpModels.Add(CopyPModel(itmbwpModel));            
            foreach (TEX itmbTex in bSkeletonIn.textures) bSkeletonOut.textures.Add(itmbTex);

            return bSkeletonOut;
        }





    }
}
