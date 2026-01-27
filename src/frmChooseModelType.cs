using KimeraCS.Core;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Forms;

namespace KimeraCS
{
    public partial class frmChooseModelType : Form
    {
        public ModelType ModelType { get; private set; } = ModelType.None;
        private Dictionary<ModelType, string> ModelLookup = new()
        {
            { ModelType.HRCSkeleton, "Field Skeleton" },
            { ModelType.AASkeleton, "Battle Skeleton" },
            { ModelType.MagicSkeleton, "Magic Skeleton" },
            { ModelType.ImportedModel, "P File" }
        };

        public frmChooseModelType()
        {
            InitializeComponent();

            foreach (var model in ModelLookup)
            {
                cbModelTypes.Items.Add(model.Value);
            }
            cbModelTypes.SelectedIndex = 0;
        }

        private void buttonOK_Click(object sender, System.EventArgs e)
        {
            var keys = ModelLookup.Keys.ToArray();
            ModelType = keys[cbModelTypes.SelectedIndex];
            DialogResult = DialogResult.OK;
            Close();
        }
    }
}
