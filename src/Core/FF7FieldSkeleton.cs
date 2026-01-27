using System.Globalization;
using System.Text;

using OpenTK.Graphics.OpenGL.Compatibility;
using OpenTK.Mathematics;
using KimeraCS.Rendering;

namespace KimeraCS.Core
{
    using static FF7FieldAnimation;
    using static FF7FieldRSDResource;

    using static FF7TEXTexture;
    using static FF7PModel;

    using static Utils;

    public static class FF7FieldSkeleton
    {

        //
        // Field Skeleton Structure
        //
        public struct FieldSkeleton
        {
            public string fileName;
            public string name;
            public int nBones;

            public List<FieldBone> bones;
            public List<TEX> textures_pool;

            public FieldSkeleton(string strfileName,
                                 bool loadGeometryQ = true,
                                 bool ignoreMissingPFiles = true,
                                 bool repairPolys = false,
                                 bool removeTextureCoords = false)
            {
                string strFileDirectoryName = (Path.GetDirectoryName(strfileName) ?? string.Empty);

                textures_pool = new List<TEX>();

                fileName = Path.GetFileNameWithoutExtension(strfileName).ToUpper();

                // Let's read HRC file into memory.
                string[] hrcString = File.ReadAllLines(strfileName);

                // Let's process the lines.
                // Skeleton name
                name = hrcString[1].Substring(10);

                // Skeleton number of bones
                nBones = Int32.Parse(hrcString[2].Substring(7));
                // Check model without skeleton

                bones = new List<FieldBone>();

                // Populate bones.
                // There is always "root" even there are no bones in some models.

                //  We will get the rest of bones rows checking if we are reading a correct line from
                //  .HRC file. There are some Field Models that can have '#' lines not useful among other things maybe.
                int i = 4;
                int numLine = 0;
                string rowOne = "", rowTwo = "", rowThree = "", rowFour = "";

                while (i < hrcString.Length)
                {
                    hrcString[i] = hrcString[i].Trim();

                    if (hrcString[i] != "" && hrcString[i][0] != '#' && hrcString[i][0] != ' ')
                    {
                        switch (numLine)
                        {
                            case 0:
                                rowOne = hrcString[i];
                                numLine++;
                                break;
                            case 1:
                                rowTwo = hrcString[i];
                                numLine++;
                                break;
                            case 2:
                                rowThree = hrcString[i];
                                numLine++;
                                break;
                            case 3:
                                rowFour = hrcString[i];
                                numLine = 0;
                                break;
                        }
                        
                        if (numLine == 0)
                        {
                            bones.Add(new FieldBone(rowOne,
                                                    rowTwo,
                                                    rowThree,
                                                    rowFour,
                                                    ref textures_pool,
                                                    loadGeometryQ,
                                                    strFileDirectoryName,
                                                    ignoreMissingPFiles,
                                                    repairPolys,
                                                    removeTextureCoords));
                        }
                    }

                    i++;
                }
            }

            //checks if the internal polys need to be repaired
            public PModel? CheckPolys()
            {
                foreach (var bone in bones)
                {
                    foreach (var rsd in bone.fRSDResources)
                    {
                        var curr = rsd.Model;
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
                    foreach (var rsd in bone.fRSDResources)
                    {
                        var curr = rsd.Model;
                        FF7PModel.RepairPolys(ref curr);
                    }
                }
            }
        }

        //
        // Field Skeleton Bone Structure
        //
        public struct FieldBone
        {
            public int nResources;
            public List<FieldRSDResource> fRSDResources;
            public double len;
            public string joint_i;
            public string joint_f;
            // added attributes
            public float resizeX;
            public float resizeY;
            public float resizeZ;

            public FieldBone(string inJointI, string inJointF, string inLen, string inRSDLine,
                             ref List<TEX> texturesPool, bool loadgeometryQ, string strFolderName,
                             bool ignoreMissingPFiles, bool repairPolys, bool removeTextureCoords)
            {
                string[] rsdRes = inRSDLine.Split();
                int i;

                len = double.Parse(inLen, CultureInfo.InvariantCulture);

                joint_i = "";
                joint_f = "";

                resizeX = 0;
                resizeY = 0;
                resizeZ = 0;

                nResources = 0;

                fRSDResources = new List<FieldRSDResource>();

                if (loadgeometryQ)
                {
                    if (inRSDLine.Length > 2)
                        nResources = int.Parse(rsdRes[0].Substring(0, inRSDLine.IndexOf(" ")));

                    joint_i = inJointI;
                    joint_f = inJointF;

                    resizeX = 1;
                    resizeY = 1;
                    resizeZ = 1;

                    // Populate resources (RSD files)
                    if (nResources > 0)
                    {
                        if (nResources > rsdRes.Length - 1) nResources = rsdRes.Length - 1;

                        for (i = 0; i < nResources && rsdRes[i + 1] != null; i++)
                        {
                            fRSDResources.Add(new FieldRSDResource(rsdRes[i + 1], ref texturesPool, strFolderName,
                                                                   ignoreMissingPFiles, repairPolys, removeTextureCoords));
                        }
                    }
                }
            }
        }



        //  ---------------------------------------------------------------------------------------------------
        //  ============================================= SAVING ==============================================
        //  ---------------------------------------------------------------------------------------------------
        public static void WriteFieldBone(ref StringBuilder strHRCContent, ref FieldBone fBone,
                                          string strDirectoryPath)
        {
            int ri;
            string strRSDList;

            FieldRSDResource tmpRSDResource;

            strHRCContent.AppendLine("");
            strHRCContent.AppendLine(fBone.joint_i);
            strHRCContent.AppendLine(fBone.joint_f);
            strHRCContent.AppendLine(fBone.len.ToString("0.0######", CultureInfo.InvariantCulture));

            strRSDList = fBone.nResources.ToString();

            if (fBone.nResources > 0)
            {
                // Write resources (if there is number involved it begins with 1, if we use "0" as first number there are issues).
                //                  first resource has no number.
                for (ri = 0; ri < fBone.nResources; ri++)
                {
                    tmpRSDResource = fBone.fRSDResources[ri];

                    strRSDList = strRSDList + " " + tmpRSDResource.res_file.ToUpper();

                    WriteRSDResource(tmpRSDResource, strDirectoryPath + "\\" + fBone.fRSDResources[ri].res_file.ToUpper() + ".RSD", strDirectoryPath);

                    if (tmpRSDResource.Model.Polys != null)
                        WriteGlobalPModel(ref tmpRSDResource.Model, strDirectoryPath + "\\" + fBone.fRSDResources[ri].Model.fileName.ToUpper());

                    fBone.fRSDResources[ri] = tmpRSDResource;
                }
            }
            else strRSDList += " ";

            strHRCContent.AppendLine(strRSDList);
        }

        public static void WriteFieldSkeleton(ref FieldSkeleton fSkeleton, string fileName)
        {
            int bi;
            StringBuilder strHRCContent = new StringBuilder();
            FieldBone tmpfBone;

            strHRCContent.AppendLine(":HEADER_BLOCK 2");
            strHRCContent.AppendLine(":SKELETON " + fSkeleton.name);
            strHRCContent.AppendLine(":BONES " + fSkeleton.nBones);

            //for (bi = 0; bi < fSkeleton.nBones; bi++)
            for (bi = 0; bi < fSkeleton.bones.Count; bi++)
            {
                tmpfBone = fSkeleton.bones[bi];
                WriteFieldBone(ref strHRCContent, ref tmpfBone, (Path.GetDirectoryName(fileName) ?? string.Empty));
                fSkeleton.bones[bi] = tmpfBone;
            }

            File.WriteAllText(fileName.ToUpper(), strHRCContent.ToString());
        }
    }
}
