using System;
using System.IO;
using System.Collections.Generic;
using System.Globalization;
using System.Windows.Forms;
using KimeraCS.Core;

namespace KimeraCS
{
    using static FF7Skeleton;
    using static FF7FieldSkeleton;
    using static FF7FieldAnimation;
    using static FF7FieldRSDResource;

    using static Utils;
    using static FileTools;

    public partial class FrmSkeletonJoints : Form
    {

        readonly private FrmSkeletonEditor frmSkelEdit;

        private static string strAttachedTo = string.Empty, strAttachedName = string.Empty;
        private static int iMode;                        // Mode: 0 - Add Joint / 1 - Edit Joint
        private static int iNumBoneSelected;


        public FrmSkeletonJoints(FrmSkeletonEditor frmSkelEdit, int iMenuMode, string strJoint)
        {
            InitializeComponent();

            this.frmSkelEdit = frmSkelEdit;
            Owner = frmSkelEdit;

            iMode = iMenuMode;
            iNumBoneSelected = -1;

            if (iMode == 1)
            {
                strAttachedTo = strJoint.Split('-')[1];
                strAttachedName = strJoint.Split('-')[0];
            }
            else
            {
                strAttachedName = "new_bone";

                if (strJoint != "")
                    strAttachedTo = strJoint.Split('-')[0];
            }
        }

        private void FillJointInfo()
        {
            if (skeleton != null)
            {
                int iBoneCounter;
                float ftmpNum = 0f;

                // Let's fill the combobox with the bones names
                cbBoneAttachedTo.Items.Clear();
                txtLength.Text = ftmpNum.ToString("F1", CultureInfo.InvariantCulture.NumberFormat);

                for (iBoneCounter = 0; iBoneCounter < skeleton.BoneCount; iBoneCounter++)
                {
                    if (skeleton.Bones[iBoneCounter].Name != strAttachedName)
                    {
                        if (iBoneCounter == 0)
                        {
                            cbBoneAttachedTo.Items.Add(skeleton.Bones[iBoneCounter].ParentName);
                        }
                        cbBoneAttachedTo.Items.Add(skeleton.Bones[iBoneCounter].Name);
                    }
                    else iNumBoneSelected = iBoneCounter;
                }

                cbBoneAttachedTo.SelectedItem = strAttachedTo;
                txtBoneName.Text = strAttachedName;

                if (iNumBoneSelected != -1)
                    txtLength.Text = skeleton.Bones[iNumBoneSelected].Length.ToString();
            }
        }

        private void FrmSkeletonJoints_Load(object sender, EventArgs e)
        {
            if (skeleton != null)
            {
                // We will fill the combo boxes with the list of bones
                FillJointInfo();

                // We will select the bones to where the joint is connected (if exist any)
                if (iMode == 1)
                {
                    txtLength.Enabled = false;

                    txtRSDFile.Enabled = false;

                    if (skeleton.Bones[iNumBoneSelected].Models.Count > 0)
                        txtRSDFile.Text = skeleton.Bones[iNumBoneSelected].Models[0].ResourceFile + ".RSD";
                    else
                        txtRSDFile.Text = "NO MODEL ASSIGNED";

                    btnPathRSDFile.Enabled = false;

                    this.Text = "Edit Joint";
                }
                else
                {
                    txtLength.Enabled = true;

                    txtRSDFile.Enabled = true;
                    txtRSDFile.Text = "";
                    btnPathRSDFile.Enabled = true;

                    this.Text = "Add Joint";
                }
            }
        }

        private void ReorderBones()
        {
            if (skeleton != null && animation != null)
            {
                int iBoneCounter, iMoveCounter, iStoppedCounter, iJumpedCounter;

                var lstTmpBones = new List<UnifiedBone>();
                var lstTmpBonesAnim = new List<UnifiedBoneRotation>();

                iBoneCounter = 0;

                while (iBoneCounter < skeleton.BoneCount)
                {
                    if (iBoneCounter == cbBoneAttachedTo.SelectedIndex)
                    {
                        // In this case the bones are moved upward

                        iStoppedCounter = iBoneCounter;
                        iMoveCounter = iNumBoneSelected;
                        while (iMoveCounter < skeleton.BoneCount - 1 &&
                               skeleton.Bones[iMoveCounter].Name == skeleton.Bones[iMoveCounter + 1].ParentName)
                        {
                            lstTmpBones.Add(new UnifiedBone(skeleton.Bones[iMoveCounter]));
                            lstTmpBonesAnim.Add(new UnifiedBoneRotation(animation.Frames[0].BoneRotations[iMoveCounter]));
                            iMoveCounter++;
                            iBoneCounter++;
                        }

                        // Add last bone if needed
                        if (iMoveCounter < skeleton.BoneCount)
                        {
                            lstTmpBones.Add(new UnifiedBone(skeleton.Bones[iMoveCounter]));
                            lstTmpBonesAnim.Add(new UnifiedBoneRotation(animation.Frames[0].BoneRotations[iMoveCounter]));
                        }

                        while (iStoppedCounter < iNumBoneSelected)
                        {
                            lstTmpBones.Add(new UnifiedBone(skeleton.Bones[iStoppedCounter]));
                            lstTmpBonesAnim.Add(new UnifiedBoneRotation(animation.Frames[0].BoneRotations[iStoppedCounter]));

                            iStoppedCounter++;
                            iBoneCounter++;
                        }
                    }
                    else if (iBoneCounter == iNumBoneSelected)
                    {
                        // In this case the bones are moved downward
                        iStoppedCounter = iBoneCounter;
                        iMoveCounter = iBoneCounter + 1;
                        iJumpedCounter = 1;

                        // Let's jump the joints we need to move.
                        while (iMoveCounter < skeleton.BoneCount &&
                               skeleton.Bones[iMoveCounter].ParentName == skeleton.Bones[iMoveCounter - 1].Name)
                        {
                            iMoveCounter++;
                            iJumpedCounter++;
                        }

                        // Add lower joints until selected bone into the new temporary list
                        while (iMoveCounter < skeleton.BoneCount &&
                               skeleton.Bones[iMoveCounter].Name != cbBoneAttachedTo.Text)
                        {
                            lstTmpBones.Add(new UnifiedBone(skeleton.Bones[iMoveCounter]));
                            lstTmpBonesAnim.Add(new UnifiedBoneRotation(animation.Frames[0].BoneRotations[iMoveCounter]));

                            iMoveCounter++;
                            iBoneCounter++;
                        }

                        // Let's add last lower joint
                        lstTmpBones.Add(new UnifiedBone(skeleton.Bones[iMoveCounter]));
                        lstTmpBonesAnim.Add(new UnifiedBoneRotation(animation.Frames[0].BoneRotations[iMoveCounter]));

                        // Add jumped bones
                        while (iJumpedCounter > 0)
                        {
                            lstTmpBones.Add(new UnifiedBone(skeleton.Bones[iStoppedCounter]));
                            lstTmpBonesAnim.Add(new UnifiedBoneRotation(animation.Frames[0].BoneRotations[iStoppedCounter]));

                            iStoppedCounter++;
                            iBoneCounter++;
                            iJumpedCounter--;
                        }

                    }
                    else
                    {
                        // This is the default case where we copy directly the bone
                        lstTmpBones.Add(new UnifiedBone(skeleton.Bones[iBoneCounter]));
                        lstTmpBonesAnim.Add(new UnifiedBoneRotation(animation.Frames[0].BoneRotations[iBoneCounter]));
                    }

                    iBoneCounter++;
                }

                // Assign ordered bones
                skeleton.Bones.Clear();
                animation.Frames[0].BoneRotations.Clear();
                for (iBoneCounter = 0; iBoneCounter < skeleton.BoneCount; iBoneCounter++)
                {
                    skeleton.Bones.Add(new UnifiedBone(lstTmpBones[iBoneCounter]));
                    animation.Frames[0]
                        .BoneRotations.Add(new UnifiedBoneRotation(lstTmpBonesAnim[iBoneCounter]));
                }
                lstTmpBones.Clear();
                lstTmpBonesAnim.Clear();
            }
        }


        //        Name ParentName       ==>     Name ParentName     ==>     Name ParentName
        //          1 	- 	0                      1 	- 	0                   1 	- 	0
        //          2 	- 	1                      2 	- 	1                   2 	- 	1
        //          3 	- 	2                      3 	- 	2                   3 	- 	2
        //          4 	- 	1                      4 	- 	1                   6   -   2
        //          5 	- 	4                      5 	- 	4                   7   -   6
        //          6 	- 	5                      6 	- 	2                   8   -   7
        //          7 	- 	6                      7 	- 	6                   4   -   1
        //          8 	- 	7                      8 	- 	7                   5   -   4
        //          9 	- 	1                      9 	- 	1                   9 	- 	1
        private void ProcessEditJoint()
        {
            if (skeleton != null)
            {
                int iBoneCounter;

                UnifiedBone tmpfBone;

                try
                {
                    // Update the bone data
                    // Assign new updated values
                    tmpfBone = new UnifiedBone(skeleton.Bones[iNumBoneSelected]);
                    tmpfBone.Name = txtBoneName.Text;
                    tmpfBone.ParentName = cbBoneAttachedTo.Text;

                    skeleton.Bones[iNumBoneSelected] = new UnifiedBone(tmpfBone);

                    // Now we need to update the other instances of "attached to" where
                    // the name is strAttached in case it has been changed.
                    for (iBoneCounter = 0; iBoneCounter < skeleton.BoneCount; iBoneCounter++)
                    {
                        if (skeleton.Bones[iBoneCounter].ParentName == strAttachedName)
                        {
                            tmpfBone = new UnifiedBone(skeleton.Bones[iBoneCounter]);
                            tmpfBone.ParentName = txtBoneName.Text;
                            skeleton.Bones[iBoneCounter] = new UnifiedBone(tmpfBone);
                        }
                    }

                    // Let's check and/or reorder if needed the joints for good rendering
                    // We must think also that we need to change the rotations.
                    ReorderBones();

                    // Update panel model and Combobox of Selected Bone in Skeleton window
                    frmSkelEdit.UpdateBones(txtBoneName.Text + "-" + cbBoneAttachedTo.Text);

                    // Refresh local combobox with Bones list
                    strAttachedTo = cbBoneAttachedTo.Text;
                    FillJointInfo();
                }
                catch (Exception ex)
                {
                    strGlobalExceptionMessage = ex.Message;

                    MessageBox.Show("There has been some error in the Edit Joint process.", "Error");
                }
            }
        }

        // Add Joint is different to Edit Joint because here we will add a new
        // bone, and if needed or selected, we will leave the choice to the user
        // to add also the .RSD model.
        private void ProcessAddJoint()
        {
            if (skeleton != null && animation != null)
            {
                FieldBone tmpfBone;
                FieldRSDResource tmpfRSDResource;

                // First let's check things to be sure we can proceed like:
                // - bone attached to selected, text bone name, float len or .rsd file)
                if (cbBoneAttachedTo.SelectedIndex == -1)
                {
                    MessageBox.Show("You must select some bone where to attach the new joint.", "Information");
                    return;
                }

                if (txtBoneName.Text == "")
                {
                    MessageBox.Show("You must write some name for the bone.", "Information");
                    return;
                }

                if (!float.TryParse(txtLength.Text, NumberStyles.Any,
                                    CultureInfo.InvariantCulture.NumberFormat, out float ftmpLen))
                {
                    MessageBox.Show("Is not possible to parse the specified length.", "Information");
                    return;
                }


                // We have to create the Bone and/or the .RSD part structure if needed.
                tmpfBone = new FieldBone()
                {
                    len = ftmpLen,
                    joint_f = cbBoneAttachedTo.Text,
                    joint_i = txtBoneName.Text,

                    nResources = 0,
                    fRSDResources = new List<FieldRSDResource>(),

                    resizeX = 1,
                    resizeY = 1,
                    resizeZ = 1,
                };

                // We have to add the .RSD part structure in the bone list of the Skeleton
                // if the user has selected one
                if (txtRSDFile.Text != "")
                {
                    tmpfBone.nResources = 1;

                    tmpfRSDResource = UserPrompts.FieldRSDLoader(Path.GetFileNameWithoutExtension(txtRSDFile.Text),
                                                               ref skeleton.Textures,
                                                               (Path.GetDirectoryName(txtRSDFile.Text) ?? string.Empty));
                    tmpfBone.fRSDResources.Add(tmpfRSDResource);
                }

                // Finally we add the Bone/Joint to the Skeleton
                iNumBoneSelected = skeleton.BoneCount;
                var bone = new UnifiedBone(tmpfBone, iNumBoneSelected, cbBoneAttachedTo.SelectedIndex);
                skeleton.Bones.Add(bone);
                //skeleton.BoneCount++;

                // We have to add the new rotation of the new Bone/Joint in the Animation
                animation.Frames[0].BoneRotations.Add(new UnifiedBoneRotation(0, 0, 0));
                animation.BoneCount++;

                // Let's check and/or reorder if needed the joints for good rendering
                // We must think also that we need to change the rotations.
                ReorderBones();

                // Update panel model and Combobox of Selected Bone in Skeleton window
                frmSkelEdit.UpdateBones(txtBoneName.Text + "-" + cbBoneAttachedTo.Text);

                // Refresh local combobox with Bones list
                strAttachedTo = cbBoneAttachedTo.Text;
                FillJointInfo();
            }
        }

        private void BtnCommit_Click(object sender, EventArgs e)
        {
            if (iMode == 1)
            {
                // Edit Joint
                ProcessEditJoint();
            }
            else
            {
                // Add Joint
                ProcessAddJoint();
            }
        }

        private void BtnPathRSDFile_Click(object sender, EventArgs e)
        {
            // Set filter options and filter index.
            openFile.Title = "Open RSD Resource";
            openFile.Filter = "FF7 RSD Resource|*.RSD|All files|*.*";
            openFile.FilterIndex = 1;
            openFile.FileName = null;

            // Check Initial Directory
            if (strGlobalPathRSDResourceFolder != null)
            {
                openFile.InitialDirectory = strGlobalPathRSDResourceFolder;
            }
            else
            {
                openFile.InitialDirectory = strGlobalPath;
            }

            try
            {
                // Process input if the user clicked OK.
                if (openFile.ShowDialog() == DialogResult.OK)
                {
                    if (File.Exists(openFile.FileName) && Path.GetExtension(openFile.FileName).ToUpper() == ".RSD")
                    {
                        // Set Global Paths
                        strGlobalPathRSDResourceFolder = Path.GetDirectoryName(openFile.FileName);

                        // Update txtPath variable
                        txtRSDFile.Text = openFile.FileName;
                    }
                }
            }
            catch (Exception ex)
            {
                strGlobalExceptionMessage = ex.Message;

                MessageBox.Show("Error opening Model file " + openFile.FileName.ToUpper() + ".",
                                "Error");
                return;
            }
        }

        private void BtnClose_Click(object sender, EventArgs e)
        {
            Close();
        }

        private void FrmSkeletonJoints_FormClosing(object sender, FormClosingEventArgs e)
        {
            cbBoneAttachedTo.SelectedIndex = -1;
            strAttachedName = "";
            strAttachedTo = "";
        }

    }
}
