using System;
using System.Windows.Forms;

namespace KimeraCS
{
    public partial class frmBattleAnimationImport : Form
    {
        public int AnimationPosition
        {
            get { return cbAnimationPosition.SelectedIndex; }
        }
        public bool InsertNew
        {
            get { return rbInsert.Checked || IsLast; }
        }
        private bool IsLast
        {
            get { return AnimationPosition == numAnims; }
        }
        private int numAnims;

        public frmBattleAnimationImport(int numAnims)
        {
            InitializeComponent();

            this.numAnims = numAnims;
            for (int i = 0; i < numAnims; i++)
            {
                cbAnimationPosition.Items.Add(i);
            }
            cbAnimationPosition.Items.Add("(Insert at end)");
            cbAnimationPosition.SelectedIndex = 0;
        }

        private void cbAnimationPosition_SelectedIndexChanged(object sender, EventArgs e)
        {
            rbInsert.Enabled = rbOverwrite.Enabled = !IsLast;
        }

        private void buttonOK_Click(object sender, EventArgs e)
        {
            DialogResult = DialogResult.OK;
            Close();
        }
    }
}
