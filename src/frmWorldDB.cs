using KimeraCS.Core;
using System;
using System.Reflection;

namespace KimeraCS
{

    using static FileTools;

    public partial class FrmWorldDB : Form
    {

        public static string strWorldFile = "", strAnimFile = "", strLocalModelName = "", strLocalAnimName = "";
        public static bool bSelectedFileFromDB = false;

        public FrmWorldDB()
        {
            InitializeComponent();
        }

        private void FrmWorldDB_Load(object sender, EventArgs e)
        {
            int mi;

            bSelectedFileFromDB = false;

            if (strWorldLGPPathSrc == "") strWorldLGPPathSrc = strGlobalPath;

            txtWorldDataDir.Text = strWorldLGPPathSrc;

            //Set Double buffering on the Grid using reflection and the bindingflags enum.
            typeof(DataGridView).InvokeMember("DoubleBuffered",
                                              BindingFlags.NonPublic | BindingFlags.Instance | BindingFlags.SetProperty, null,
                                              dgvWorld, new object[] { true });

            // World
            if (bWorldDBLoaded)
            {
                dgvWorld.Rows.Clear();
                dgvWorld.Refresh();

                for (mi = 0; mi < lstWorldLGPRegisters.Count; mi++)
                {
                    dgvWorld.Rows.Insert(mi, mi,
                                            lstWorldLGPRegisters[mi].fileName,
                                            lstWorldLGPRegisters[mi].modelName);
                }

                if (strLocalModelName == "")
                {
                    dgvWorld.Rows[0].Selected = true;
                }
                else
                {
                    DataGridViewRow dgvRow = dgvWorld.Rows
                                                .Cast<DataGridViewRow>()
                                                .Where(r => (r.Cells[1].Value?.ToString() ?? string.Empty).Equals(strLocalModelName))
                                                .First();

                    dgvWorld.CurrentCell = dgvWorld.Rows[dgvRow.Index].Cells[1];
                    dgvWorld.Rows[dgvRow.Index].Selected = true;
                }
            }

            UpdatelbAnimation();

            // Select lbAnimation index and value
            if (strLocalAnimName == "")
            {
                lbAnimation.SelectedIndex = 0;
                lbAnimation.SelectedValue = 0;
            }

            if (strLocalAnimName == "") lbAnimation.SetSelected(0, true);
            else lbAnimation.SetSelected(lbAnimation.Items.IndexOf(strLocalAnimName), true);
        }

        private void BtnSaveWorldDataDir_Click(object sender, EventArgs e)
        {
            WriteCFGFile();
        }

        private void BtnSelectDirBrowser_Click(object sender, EventArgs e)
        {
            DialogResult result;
            string path;

            using (var folderDialog = new FolderBrowserDialog())
            {
                folderDialog.Description =
                    "Select the Folder where WORLD_xx.LGP file has been extracted:";
                result = folderDialog.ShowDialog();
                path = folderDialog.SelectedPath;
            }

            if (result == DialogResult.OK)
            {
                if (!string.IsNullOrEmpty(path))
                {
                    // Put Global folder for input unswizzled.
                    strWorldLGPPathSrc = path;
                    txtWorldDataDir.Text = path;
                }
            }
        }

        private void UpdatelbAnimation()
        {
            lbAnimation.Items.Clear();

            if (dgvWorld.SelectedRows.Count > 0)
            {
                int index = dgvWorld.SelectedRows[0].Index;
                for (int i = 0; i < lstWorldLGPRegisters[index].lstAnims.Count; i++)
                {
                    lbAnimation.Items.Add(lstWorldLGPRegisters[index].lstAnims[i]);
                }

                strLocalModelName = lstWorldLGPRegisters[index].fileName;
                lbAnimation.SetSelected(0, true);
            }
        }

        private void dgvWorld_CellClick(object sender, DataGridViewCellEventArgs e)
        {
            UpdatelbAnimation();
        }

        private void FrmWorldDB_FormClosed(object sender, FormClosedEventArgs e)
        {
            int index = dgvWorld.SelectedRows[0].Index;
            strLocalModelName = lstWorldLGPRegisters[index].fileName;
            strLocalAnimName = lbAnimation.Text;

            if (!bSelectedFileFromDB)
            {
                strWorldFile = "";
                strAnimFile = "";
            }
        }

        private void BtnLoadModelAnimation_Click(object sender, EventArgs e)
        {
            bSelectedFileFromDB = false;
            string name = string.Empty;

            if (dgvWorld.SelectedRows.Count > 0)
            {
                int index = dgvWorld.SelectedRows[0].Index;
                name = lstWorldLGPRegisters[index].fileName;
                strWorldFile = txtWorldDataDir.Text + "\\" + name + ".HRC";
            }
            else
            {
                MessageBox.Show("You have not selected any World Model.", "Information", MessageBoxButtons.OK);
                return;
            }

            if (lbAnimation.Text != "")
                strAnimFile = txtWorldDataDir.Text + "\\" + lbAnimation.Text + ".A";
            else
            {
                MessageBox.Show("You have not selected any World Animation for the model.", "Information", MessageBoxButtons.OK);
                return;

            }

            // Check existance of the files.
            if (!File.Exists(strWorldFile))
            {
                MessageBox.Show("The file selected as World Model " + name + " does not exist.",
                                "Error");
                return;
            }

            if (!File.Exists(strAnimFile))
            {
                MessageBox.Show("The file selected as World Animation " + lbAnimation.SelectedItem + " does not exist.",
                                "Error");
                return;
            }

            bSelectedFileFromDB = true;

            Close();
        }

        private void BtnClose_Click(object sender, EventArgs e)
        {
            Close();
        }
    }
}
