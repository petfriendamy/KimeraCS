using OpenTK.Graphics.OpenGL.Compatibility;
using OpenTK.Mathematics;

namespace KimeraCS.Core
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
            public SkeletonType skeletonType;              //  0 - Enemy Model, 1 - Battle Location, 2 - PC Battle Model?
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
                strFileDirectoryName = (Path.GetDirectoryName(strFullFileName) ?? string.Empty);

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
                        skeletonType = (SkeletonType)memReader.ReadInt32();
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
                strFileDirectoryName = (Path.GetDirectoryName(strFullFileName) ?? string.Empty);

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
                        skeletonType = (SkeletonType)memReader.ReadInt32();
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
            strFullDirectoryName = (Path.GetDirectoryName(strFileName) ?? string.Empty);

            // Writer Main Battle file (AA)
            using (MemoryStream fileMemory = new MemoryStream(fileBuffer))
            {
                using (var memWriter = new BinaryWriter(fileMemory))
                {
                    memWriter.Write((int)bSkeleton.skeletonType);
                    memWriter.Write(bSkeleton.unk1);
                    memWriter.Write(bSkeleton.unk2);

                    if (bSkeleton.IsBattleLocation)
                        memWriter.Write(0);
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
            strFullDirectoryName = (Path.GetDirectoryName(strFileName) ?? string.Empty);

            // Writer Main Magic file (.D)
            using (MemoryStream fileMemory = new MemoryStream(fileBuffer))
            {
                using (var memWriter = new BinaryWriter(fileMemory))
                {
                    memWriter.Write((int)bSkeleton.skeletonType);
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
    }
}
