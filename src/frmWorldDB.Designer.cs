
namespace KimeraCS
{
    partial class FrmWorldDB
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            System.ComponentModel.ComponentResourceManager resources = new System.ComponentModel.ComponentResourceManager(typeof(FrmWorldDB));
            lblModel = new Label();
            btnLoadModelAnimation = new Button();
            lblAnimations = new Label();
            lbAnimation = new ListBox();
            lblWorldDataDir = new Label();
            txtWorldDataDir = new TextBox();
            btnSelectDirBrowser = new Button();
            btnSaveWorldDataDir = new Button();
            lblLine = new Label();
            btnClose = new Button();
            dgvWorld = new DataGridView();
            colNumber = new DataGridViewTextBoxColumn();
            colFileName = new DataGridViewTextBoxColumn();
            colName = new DataGridViewTextBoxColumn();
            ((System.ComponentModel.ISupportInitialize)dgvWorld).BeginInit();
            SuspendLayout();
            // 
            // lblModel
            // 
            lblModel.AutoSize = true;
            lblModel.ForeColor = SystemColors.Control;
            lblModel.Location = new Point(18, 17);
            lblModel.Name = "lblModel";
            lblModel.Size = new Size(46, 15);
            lblModel.TabIndex = 0;
            lblModel.Text = "Models";
            // 
            // btnLoadModelAnimation
            // 
            btnLoadModelAnimation.Anchor = AnchorStyles.Top | AnchorStyles.Left | AnchorStyles.Right;
            btnLoadModelAnimation.Location = new Point(22, 224);
            btnLoadModelAnimation.Margin = new Padding(3, 2, 3, 2);
            btnLoadModelAnimation.Name = "btnLoadModelAnimation";
            btnLoadModelAnimation.Size = new Size(541, 45);
            btnLoadModelAnimation.TabIndex = 3;
            btnLoadModelAnimation.Text = "Load Model and Animation";
            btnLoadModelAnimation.UseVisualStyleBackColor = true;
            btnLoadModelAnimation.Click += BtnLoadModelAnimation_Click;
            // 
            // lblAnimations
            // 
            lblAnimations.AutoSize = true;
            lblAnimations.ForeColor = SystemColors.Control;
            lblAnimations.Location = new Point(397, 17);
            lblAnimations.Name = "lblAnimations";
            lblAnimations.Size = new Size(68, 15);
            lblAnimations.TabIndex = 4;
            lblAnimations.Text = "Animations";
            // 
            // lbAnimation
            // 
            lbAnimation.Anchor = AnchorStyles.Top | AnchorStyles.Right;
            lbAnimation.FormattingEnabled = true;
            lbAnimation.Location = new Point(397, 36);
            lbAnimation.Margin = new Padding(3, 2, 3, 2);
            lbAnimation.Name = "lbAnimation";
            lbAnimation.Size = new Size(166, 169);
            lbAnimation.TabIndex = 5;
            // 
            // lblWorldDataDir
            // 
            lblWorldDataDir.AutoSize = true;
            lblWorldDataDir.ForeColor = SystemColors.Control;
            lblWorldDataDir.Location = new Point(20, 300);
            lblWorldDataDir.Name = "lblWorldDataDir";
            lblWorldDataDir.Size = new Size(120, 15);
            lblWorldDataDir.TabIndex = 6;
            lblWorldDataDir.Text = "World Data Directory:";
            // 
            // txtWorldDataDir
            // 
            txtWorldDataDir.Anchor = AnchorStyles.Top | AnchorStyles.Left | AnchorStyles.Right;
            txtWorldDataDir.Location = new Point(22, 319);
            txtWorldDataDir.Margin = new Padding(3, 2, 3, 2);
            txtWorldDataDir.Name = "txtWorldDataDir";
            txtWorldDataDir.Size = new Size(501, 23);
            txtWorldDataDir.TabIndex = 7;
            // 
            // btnSelectDirBrowser
            // 
            btnSelectDirBrowser.Anchor = AnchorStyles.Top | AnchorStyles.Right;
            btnSelectDirBrowser.Location = new Point(528, 319);
            btnSelectDirBrowser.Margin = new Padding(3, 2, 3, 2);
            btnSelectDirBrowser.Name = "btnSelectDirBrowser";
            btnSelectDirBrowser.Size = new Size(35, 23);
            btnSelectDirBrowser.TabIndex = 8;
            btnSelectDirBrowser.Text = "...";
            btnSelectDirBrowser.UseVisualStyleBackColor = true;
            btnSelectDirBrowser.Click += BtnSelectDirBrowser_Click;
            // 
            // btnSaveWorldDataDir
            // 
            btnSaveWorldDataDir.Location = new Point(22, 353);
            btnSaveWorldDataDir.Margin = new Padding(3, 2, 3, 2);
            btnSaveWorldDataDir.Name = "btnSaveWorldDataDir";
            btnSaveWorldDataDir.Size = new Size(314, 39);
            btnSaveWorldDataDir.TabIndex = 9;
            btnSaveWorldDataDir.Text = "Save World Data Directory";
            btnSaveWorldDataDir.UseVisualStyleBackColor = true;
            btnSaveWorldDataDir.Click += BtnSaveWorldDataDir_Click;
            // 
            // lblLine
            // 
            lblLine.Anchor = AnchorStyles.Top | AnchorStyles.Left | AnchorStyles.Right;
            lblLine.BorderStyle = BorderStyle.Fixed3D;
            lblLine.Location = new Point(22, 284);
            lblLine.Margin = new Padding(4, 0, 4, 0);
            lblLine.Name = "lblLine";
            lblLine.Size = new Size(541, 2);
            lblLine.TabIndex = 10;
            // 
            // btnClose
            // 
            btnClose.DialogResult = DialogResult.Cancel;
            btnClose.Location = new Point(343, 355);
            btnClose.Margin = new Padding(4);
            btnClose.Name = "btnClose";
            btnClose.Size = new Size(220, 39);
            btnClose.TabIndex = 11;
            btnClose.Text = "Close";
            btnClose.UseVisualStyleBackColor = true;
            btnClose.Click += BtnClose_Click;
            // 
            // dgvWorld
            // 
            dgvWorld.AllowUserToAddRows = false;
            dgvWorld.AllowUserToDeleteRows = false;
            dgvWorld.AllowUserToResizeColumns = false;
            dgvWorld.AllowUserToResizeRows = false;
            dgvWorld.Anchor = AnchorStyles.Top | AnchorStyles.Left | AnchorStyles.Right;
            dgvWorld.ClipboardCopyMode = DataGridViewClipboardCopyMode.Disable;
            dgvWorld.ColumnHeadersHeight = 20;
            dgvWorld.Columns.AddRange(new DataGridViewColumn[] { colNumber, colFileName, colName });
            dgvWorld.EditMode = DataGridViewEditMode.EditProgrammatically;
            dgvWorld.Location = new Point(18, 36);
            dgvWorld.MultiSelect = false;
            dgvWorld.Name = "dgvWorld";
            dgvWorld.ReadOnly = true;
            dgvWorld.RowHeadersVisible = false;
            dgvWorld.SelectionMode = DataGridViewSelectionMode.FullRowSelect;
            dgvWorld.Size = new Size(373, 169);
            dgvWorld.TabIndex = 12;
            dgvWorld.CellClick += dgvWorld_CellClick;
            // 
            // colNumber
            // 
            colNumber.Frozen = true;
            colNumber.HeaderText = "";
            colNumber.Name = "colNumber";
            colNumber.ReadOnly = true;
            colNumber.Width = 60;
            // 
            // colFileName
            // 
            colFileName.HeaderText = "File Name";
            colFileName.Name = "colFileName";
            colFileName.ReadOnly = true;
            // 
            // colName
            // 
            colName.AutoSizeMode = DataGridViewAutoSizeColumnMode.Fill;
            colName.HeaderText = "Name";
            colName.Name = "colName";
            colName.ReadOnly = true;
            // 
            // FrmWorldDB
            // 
            AcceptButton = btnLoadModelAnimation;
            AutoScaleDimensions = new SizeF(7F, 15F);
            AutoScaleMode = AutoScaleMode.Font;
            BackColor = SystemColors.ControlDarkDark;
            CancelButton = btnClose;
            ClientSize = new Size(584, 407);
            Controls.Add(dgvWorld);
            Controls.Add(btnClose);
            Controls.Add(lblLine);
            Controls.Add(btnSaveWorldDataDir);
            Controls.Add(btnSelectDirBrowser);
            Controls.Add(txtWorldDataDir);
            Controls.Add(lblWorldDataDir);
            Controls.Add(lbAnimation);
            Controls.Add(lblAnimations);
            Controls.Add(btnLoadModelAnimation);
            Controls.Add(lblModel);
            DoubleBuffered = true;
            FormBorderStyle = FormBorderStyle.FixedDialog;
            Icon = (Icon)resources.GetObject("$this.Icon");
            KeyPreview = true;
            Margin = new Padding(3, 2, 3, 2);
            MaximizeBox = false;
            MinimizeBox = false;
            Name = "FrmWorldDB";
            StartPosition = FormStartPosition.CenterParent;
            Text = "FF7 World Database";
            FormClosed += FrmWorldDB_FormClosed;
            Load += FrmWorldDB_Load;
            ((System.ComponentModel.ISupportInitialize)dgvWorld).EndInit();
            ResumeLayout(false);
            PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Label lblModel;
        private System.Windows.Forms.Button btnLoadModelAnimation;
        private System.Windows.Forms.Label lblAnimations;
        private System.Windows.Forms.ListBox lbAnimation;
        private System.Windows.Forms.Label lblWorldDataDir;
        private System.Windows.Forms.TextBox txtWorldDataDir;
        private System.Windows.Forms.Button btnSelectDirBrowser;
        private System.Windows.Forms.Button btnSaveWorldDataDir;
        private System.Windows.Forms.Label lblLine;
        private System.Windows.Forms.Button btnClose;
        private DataGridView dgvWorld;
        private DataGridViewTextBoxColumn colNumber;
        private DataGridViewTextBoxColumn colFileName;
        private DataGridViewTextBoxColumn colName;
    }
}