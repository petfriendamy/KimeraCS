using System;
using System.Collections.Generic;
using System.Windows.Forms;
using KimeraCS.Core;

namespace KimeraCS
{
    using static FF7Skeleton;
    using static FileTools;
    using static Utils;

    public partial class FrmStatistics : Form
    {

        private static string strFileName = string.Empty;

        public FrmStatistics()
        {
            InitializeComponent();
        }

        private void FrmStatistics_Load(object sender, EventArgs e)
        {
            if (skeleton != null)
            {
                int iBoneIdx, iRSDIdx, iGroupIdx, iPolyIdx, iVertIdx, iTexCount, iWpnIdx, iWpnCounted;
                int iTotalVerts, iTotalPolys, iTotalTexCoords, iRSDVertsUsage, iTotalVertsUsage;

                int[] arrVertexUsage;
                HashSet<string> hsTexList = new HashSet<string>();

                iTotalVerts = iTotalPolys = iTotalVertsUsage = iTotalTexCoords = iWpnCounted = 0;

                switch (modelType)
                {
                    // FIELD SKELETON
                    case ModelType.HRCSkeleton:
                        strFileName = skeleton.Name + ".TXT";

                        rtbStats.Text = "Field Model:\t\t" + skeleton.Name + "\n";
                        rtbStats.Text += "Number of Bones:\t" + skeleton.Bones.Count + "\n\n";

                        rtbStats.Text += "BONES\t\t\t\n";

                        for (iBoneIdx = 0; iBoneIdx < skeleton.Bones.Count; iBoneIdx++)
                        {
                            rtbStats.Text += "Bone " + iBoneIdx.ToString("00") + "\t\t\t" +
                                             "Joint:\t" + skeleton.Bones[iBoneIdx].Name + " - " +
                                                          skeleton.Bones[iBoneIdx].ParentName + "\n";
                            rtbStats.Text += "Number of RSD:\t\t" +
                                             skeleton.Bones[iBoneIdx].Models.Count.ToString("00") +
                                             "\n";

                            if (skeleton.Bones[iBoneIdx].Models.Count == 0) rtbStats.Text += "\n";

                            for (iRSDIdx = 0; iRSDIdx < skeleton.Bones[iBoneIdx].Models.Count; iRSDIdx++)
                            {
                                rtbStats.Text += "Resource Name:\t\t" + skeleton.Bones[iBoneIdx].Models[iRSDIdx].ResourceFile + ".RSD" + "\n";
                                rtbStats.Text += "Model Name:\t\t" + skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.fileName + "\n";
                                rtbStats.Text += "Vertices: " + skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.Verts.Length.ToString("00000") + "\t" +
                                                 "Triangles: " + skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.Polys.Length.ToString("00000") + "\t\t" +
                                                 "TexCoords: " + skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.TexCoords.Length.ToString("00000") + "\n";

                                iTotalVerts += skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.Verts.Length;
                                iTotalPolys += skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.Polys.Length;
                                iTotalTexCoords += skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.TexCoords.Length;

                                // Calculate the real usage of the vertices.
                                // Different polys can use the same vertex, so, this is more load for the GPU.
                                arrVertexUsage = new int[skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.Verts.Length];
                                iRSDVertsUsage = 0;

                                for (iGroupIdx = 0;
                                     iGroupIdx < skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.Header.numGroups;
                                     iGroupIdx++)
                                {

                                    for (iPolyIdx = skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.Groups[iGroupIdx].offsetPoly;
                                         iPolyIdx < skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.Groups[iGroupIdx].numPoly;
                                         iPolyIdx++)
                                    {
                                        for (iVertIdx = 0; iVertIdx < 3; iVertIdx++)
                                        {
                                            arrVertexUsage[skeleton.Bones[iBoneIdx].Models[iRSDIdx].
                                                                     Model.Polys[iPolyIdx].Verts[iVertIdx]] += 1;
                                        }
                                    }
                                }

                                for (iVertIdx = 0; iVertIdx < arrVertexUsage.Length; iVertIdx++)
                                {
                                    iRSDVertsUsage += arrVertexUsage[iVertIdx];
                                }

                                iTotalVertsUsage += iRSDVertsUsage;

                                rtbStats.Text += "Vertices Usage:\t" + iRSDVertsUsage.ToString() + "\n";

                                // We need to add the Textures info if any
                                if (skeleton.Bones[iBoneIdx].Models[iRSDIdx].TextureCount > 0)
                                {

                                    rtbStats.Text += "Textures Used:\t\t";

                                    for (iTexCount = 0;
                                         iTexCount < skeleton.Bones[iBoneIdx].Models[iRSDIdx].TextureCount;
                                         iTexCount++)
                                    {
                                        rtbStats.Text += skeleton.Bones[iBoneIdx].Models[iRSDIdx].Textures[iTexCount].TEXfileName;

                                        if (iTexCount < skeleton.Bones[iBoneIdx].Models[iRSDIdx].TextureCount - 1) rtbStats.Text += ", ";

                                        hsTexList.Add(skeleton.Bones[iBoneIdx].Models[iRSDIdx].Textures[iTexCount].TEXfileName + "_" +
                                                      skeleton.Bones[iBoneIdx].Models[iRSDIdx].Textures[iTexCount].width.ToString() + " x " +
                                                      skeleton.Bones[iBoneIdx].Models[iRSDIdx].Textures[iTexCount].height.ToString());
                                    }

                                    rtbStats.Text += "\n\n";
                                }

                            }

                        }

                        // List of textures and its size
                        if (hsTexList.Count > 0)
                        {
                            rtbStats.Text += "\nTEXTURE LIST\n";

                            foreach (string itmTex in hsTexList)
                            {
                                rtbStats.Text += "Name:\t\t" + itmTex.Split('_')[0] + "\t\t\t" +
                                                 "Size:\t" + itmTex.Split('_')[1] + "\n";
                            }

                            rtbStats.Text += "\n";
                        }


                        // TOTALS
                        rtbStats.Text += "TOTAL TRIANGLES:\t\t" + iTotalPolys.ToString() + "\n";
                        rtbStats.Text += "TOTAL VERTICES:\t\t" + iTotalVerts.ToString() + "\n";
                        rtbStats.Text += "TOTAL TEX COORDS (UV):\t" + iTotalTexCoords.ToString() + "\n";
                        rtbStats.Text += "TOTAL VERTICES USAGE:\t" + iTotalVertsUsage.ToString() + "\n";

                        break;

                    // BATTLE/MAGIC SKELETON
                    case ModelType.AASkeleton:
                    case ModelType.MagicSkeleton:
                        strFileName = skeleton.Name;

                        // Put battle model type
                        if (modelType == ModelType.MagicSkeleton)
                            rtbStats.Text = "Magic Model:";
                        else if (skeleton.IsBattleLocation)
                            rtbStats.Text = "Battle Model (Location):";
                        else if (skeleton.Name[0] + skeleton.Name[1] >= 165)
                            rtbStats.Text = "Battle Model (Main):";
                        else
                            rtbStats.Text = "Battle Model (Enemy)";

                        rtbStats.Text += "\t\t" + skeleton.Name + "\n";
                        rtbStats.Text += "Number of Bones:\t" + skeleton.BoneCount.ToString() + "\n\n";

                        rtbStats.Text += "BONES\t\t\t\n";

                        for (iBoneIdx = 0; iBoneIdx < skeleton.BoneCount; iBoneIdx++)
                        {
                            rtbStats.Text += "Bone " + iBoneIdx.ToString("00") + "\t\t\t";
                            rtbStats.Text += "Number of Models:\t\t" +
                                             skeleton.Bones[iBoneIdx].Models.Count.ToString("00") + "\n";

                            if (skeleton.Bones[iBoneIdx].HasModel)
                            {

                                for (iRSDIdx = 0; iRSDIdx < skeleton.Bones[iBoneIdx].Models.Count; iRSDIdx++)
                                {
                                    rtbStats.Text +=
                                        "Model Name:\t\t" + skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.fileName + "\n";

                                    rtbStats.Text += "Vertices: " + skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.Verts.Length.ToString("00000") + "\t" +
                                                     "Triangles: " + skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.Polys.Length.ToString("00000") + "\t\t" +
                                                     "TexCoords: " + skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.TexCoords.Length.ToString("00000") + "\n";

                                    iTotalVerts += skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.Verts.Length;
                                    iTotalPolys += skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.Polys.Length;
                                    iTotalTexCoords += skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.TexCoords.Length;

                                    // Calculate the real usage of the vertices.
                                    // Different polys can use the same vertex, so, this is more load for the GPU.
                                    arrVertexUsage = new int[skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.Verts.Length];
                                    iRSDVertsUsage = 0;

                                    for (iGroupIdx = 0;
                                            iGroupIdx < skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.Header.numGroups;
                                            iGroupIdx++)
                                    {

                                        for (iPolyIdx = skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.Groups[iGroupIdx].offsetPoly;
                                                iPolyIdx < skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.Groups[iGroupIdx].numPoly +
                                                           skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.Groups[iGroupIdx].offsetPoly;
                                                iPolyIdx++)
                                        {
                                            for (iVertIdx = 0; iVertIdx < 3; iVertIdx++)
                                            {
                                                arrVertexUsage[skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.
                                                                            Polys[iPolyIdx].Verts[iVertIdx]] += 1;
                                            }
                                        }
                                    }

                                    for (iVertIdx = 0; iVertIdx < arrVertexUsage.Length; iVertIdx++)
                                    {
                                        iRSDVertsUsage += arrVertexUsage[iVertIdx];
                                    }

                                    iTotalVertsUsage += iRSDVertsUsage;

                                    rtbStats.Text += "Vertices Usage:\t" + iRSDVertsUsage.ToString() + "\n";

                                    // We need to add the Textures info if any
                                    if (skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.TexCoords.Length > 0)
                                    {
                                        rtbStats.Text += "Textures Used:\t\t";

                                        for (iGroupIdx = 0;
                                                iGroupIdx < skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.Groups.Length;
                                                iGroupIdx++)
                                        {
                                            if (skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.Groups[iGroupIdx].texFlag == 1)
                                            {
                                                rtbStats.Text += skeleton.Textures[skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.
                                                                        Groups[iGroupIdx].texID].TEXfileName;

                                                if (iGroupIdx < skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.Groups.Length - 1)
                                                    rtbStats.Text += ", ";

                                                hsTexList.Add(skeleton.Textures[skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.
                                                                        Groups[iGroupIdx].texID].TEXfileName + "_" +
                                                                skeleton.Textures[skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.
                                                                        Groups[iGroupIdx].texID].width.ToString() + " x " +
                                                                skeleton.Textures[skeleton.Bones[iBoneIdx].Models[iRSDIdx].Model.
                                                                        Groups[iGroupIdx].texID].height.ToString());
                                            }
                                        }

                                        rtbStats.Text += "\n\n";
                                    }
                                    else rtbStats.Text += "\n";
                                }
                            }
                            else rtbStats.Text += "\n";

                        }

                        //// WEAPONS
                        // Let add the weapons if there are any
                        if (skeleton.WeaponCount > 0)
                        {
                            rtbStats.Text += "\nWEAPONS BONES\n";
                            rtbStats.Text += "Number of Weapons:" + skeleton.WeaponCount.ToString("00") + "\n\n";

                            for (iWpnIdx = 0; iWpnIdx < skeleton.Weapons.Count; iWpnIdx++)
                            {

                                if (skeleton.Weapons[iWpnIdx].Header.numGroups > 0)
                                {
                                    iWpnCounted++;

                                    rtbStats.Text +=
                                        "Model Name:\t\t" + skeleton.Weapons[iWpnIdx].fileName + "\n";

                                    rtbStats.Text += "Vertices:\t" + skeleton.Weapons[iWpnIdx].Verts.Length.ToString() + "\t\t" +
                                                     "Triangles:\t" + skeleton.Weapons[iWpnIdx].Polys.Length.ToString() + "\t\t" +
                                                     "TexCoords:\t" + skeleton.Weapons[iWpnIdx].TexCoords.Length.ToString() + "\n";

                                    iTotalVerts += skeleton.Weapons[iWpnIdx].Verts.Length;
                                    iTotalPolys += skeleton.Weapons[iWpnIdx].Polys.Length;
                                    iTotalTexCoords += skeleton.Weapons[iWpnIdx].TexCoords.Length;

                                    // Calculate the real usage of the vertices.
                                    // Different polys can use the same vertex, so, this is more load for the GPU.
                                    arrVertexUsage = new int[skeleton.Weapons[iWpnIdx].Verts.Length];
                                    iRSDVertsUsage = 0;

                                    for (iGroupIdx = 0;
                                         iGroupIdx < skeleton.Weapons[iWpnIdx].Header.numGroups;
                                         iGroupIdx++)
                                    {

                                        for (iPolyIdx = skeleton.Weapons[iWpnIdx].Groups[iGroupIdx].offsetPoly;
                                             iPolyIdx < skeleton.Weapons[iWpnIdx].Groups[iGroupIdx].numPoly +
                                                        skeleton.Weapons[iWpnIdx].Groups[iGroupIdx].offsetPoly;
                                             iPolyIdx++)
                                        {
                                            for (iVertIdx = 0; iVertIdx < 3; iVertIdx++)
                                            {
                                                arrVertexUsage[skeleton.Weapons[iWpnIdx].
                                                                            Polys[iPolyIdx].Verts[iVertIdx]] += 1;
                                            }
                                        }
                                    }

                                    for (iVertIdx = 0; iVertIdx < arrVertexUsage.Length; iVertIdx++)
                                    {
                                        iRSDVertsUsage += arrVertexUsage[iVertIdx];
                                    }

                                    iTotalVertsUsage += iRSDVertsUsage;

                                    rtbStats.Text += "Vertices Usage:\t" + iRSDVertsUsage.ToString() + "\n";

                                    // We need to add the Textures info if any
                                    if (skeleton.Weapons[iWpnIdx].TexCoords.Length > 0)
                                    {
                                        rtbStats.Text += "Textures Used:\t\t";

                                        for (iGroupIdx = 0;
                                             iGroupIdx < skeleton.Weapons[iWpnIdx].Groups.Length;
                                             iGroupIdx++)
                                        {
                                            if (skeleton.Weapons[iWpnIdx].Groups[iGroupIdx].texFlag == 1)
                                            {
                                                rtbStats.Text += skeleton.Textures[skeleton.Weapons[iWpnIdx].
                                                                                    Groups[iGroupIdx].texID].TEXfileName;

                                                if (iGroupIdx < skeleton.Weapons[iWpnIdx].Groups.Length - 1)
                                                    rtbStats.Text += ", ";

                                                hsTexList.Add(skeleton.Textures[skeleton.Weapons[iWpnIdx].
                                                                    Groups[iGroupIdx].texID].TEXfileName + "_" +
                                                              skeleton.Textures[skeleton.Weapons[iWpnIdx].
                                                                    Groups[iGroupIdx].texID].width.ToString() + " x " +
                                                              skeleton.Textures[skeleton.Weapons[iWpnIdx].
                                                                    Groups[iGroupIdx].texID].height.ToString());
                                            }
                                        }

                                        rtbStats.Text += "\n\n";
                                    }
                                    else rtbStats.Text += "\n";

                                }
                            }

                            if (iWpnCounted == 0)
                            {
                                rtbStats.Text += "----------------------------------------------------\n";
                                rtbStats.Text += "WARNING: There are NOT weapons loaded and the Battle " +
                                                 "Model should have weapons.\n";
                                rtbStats.Text += "----------------------------------------------------\n\n";
                            }

                            if (iWpnCounted > 0 && iWpnCounted != skeleton.WeaponCount)
                            {
                                rtbStats.Text += "----------------------------------------------------\n";
                                rtbStats.Text += "WARNING: There number of weapons loaded are not the same " +
                                                 "as the number of weapons of the Battle Model.\n";
                                rtbStats.Text += "----------------------------------------------------\n\n";
                            }
                        }


                        // List of textures and its size
                        if (hsTexList.Count > 0)
                        {
                            rtbStats.Text += "\nTEXTURE LIST\n";

                            foreach (string itmTex in hsTexList)
                            {
                                rtbStats.Text += "Name:\t\t" + itmTex.Split('_')[0] + "\t\t\t" +
                                                 "Size:\t" + itmTex.Split('_')[1] + "\n";
                            }

                            rtbStats.Text += "\n";
                        }


                        // TOTALS
                        rtbStats.Text += "TOTAL TRIANGLES:\t\t" + iTotalPolys.ToString() + "\n";
                        rtbStats.Text += "TOTAL VERTICES:\t\t" + iTotalVerts.ToString() + "\n";
                        rtbStats.Text += "TOTAL TEX COORDS (UV):\t" + iTotalTexCoords.ToString() + "\n";
                        rtbStats.Text += "TOTAL VERTICES USAGE:\t" + iTotalVertsUsage.ToString() + "\n";

                        break;

                }

                if (strFileName.Length > 0) Text += " - " + strFileName;
            }
        }

        private void BtnClose_Click(object sender, EventArgs e)
        {
            Close();
        }

        private void BtnSaveStats_Click(object sender, EventArgs e)
        {
            saveFile.Title = "Save Statistics As...";
            saveFile.Filter = "Plain Text file|*.TXT|All files|*.*";
            saveFile.FilterIndex = 1;
            saveFile.InitialDirectory = strGlobalPathFieldSkeletonFolder;

            saveFile.FileName = strFileName;

            try
            {
                // Process input if the user clicked OK.
                if (saveFile.ShowDialog() == DialogResult.OK)
                {
                    rtbStats.SaveFile(saveFile.FileName, RichTextBoxStreamType.PlainText);
                }
            }
            catch (Exception ex)
            {
                strGlobalExceptionMessage = ex.Message;

                MessageBox.Show("Error saving statistics file.", "Error");
                return;
            }
        }

    }
}
