namespace KimeraCS
{
    partial class frmBattleAnimationImport
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
            cbAnimationPosition = new System.Windows.Forms.ComboBox();
            label1 = new System.Windows.Forms.Label();
            rbInsert = new System.Windows.Forms.RadioButton();
            rbOverwrite = new System.Windows.Forms.RadioButton();
            buttonOK = new System.Windows.Forms.Button();
            buttonCancel = new System.Windows.Forms.Button();
            SuspendLayout();
            // 
            // cbAnimationPosition
            // 
            cbAnimationPosition.Anchor = System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left | System.Windows.Forms.AnchorStyles.Right;
            cbAnimationPosition.DropDownStyle = System.Windows.Forms.ComboBoxStyle.DropDownList;
            cbAnimationPosition.FormattingEnabled = true;
            cbAnimationPosition.Location = new System.Drawing.Point(12, 27);
            cbAnimationPosition.Name = "cbAnimationPosition";
            cbAnimationPosition.Size = new System.Drawing.Size(260, 23);
            cbAnimationPosition.TabIndex = 0;
            cbAnimationPosition.SelectedIndexChanged += cbAnimationPosition_SelectedIndexChanged;
            // 
            // label1
            // 
            label1.AutoSize = true;
            label1.Location = new System.Drawing.Point(12, 9);
            label1.Name = "label1";
            label1.Size = new System.Drawing.Size(115, 15);
            label1.TabIndex = 1;
            label1.Text = "Insert animation at...";
            // 
            // rbInsert
            // 
            rbInsert.AutoSize = true;
            rbInsert.Location = new System.Drawing.Point(94, 56);
            rbInsert.Name = "rbInsert";
            rbInsert.Size = new System.Drawing.Size(54, 19);
            rbInsert.TabIndex = 2;
            rbInsert.Text = "Insert";
            rbInsert.UseVisualStyleBackColor = true;
            // 
            // rbOverwrite
            // 
            rbOverwrite.AutoSize = true;
            rbOverwrite.Checked = true;
            rbOverwrite.Location = new System.Drawing.Point(12, 56);
            rbOverwrite.Name = "rbOverwrite";
            rbOverwrite.Size = new System.Drawing.Size(76, 19);
            rbOverwrite.TabIndex = 3;
            rbOverwrite.TabStop = true;
            rbOverwrite.Text = "Overwrite";
            rbOverwrite.UseVisualStyleBackColor = true;
            // 
            // buttonOK
            // 
            buttonOK.Anchor = System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Right;
            buttonOK.Location = new System.Drawing.Point(197, 86);
            buttonOK.Name = "buttonOK";
            buttonOK.Size = new System.Drawing.Size(75, 23);
            buttonOK.TabIndex = 4;
            buttonOK.Text = "OK";
            buttonOK.UseVisualStyleBackColor = true;
            buttonOK.Click += buttonOK_Click;
            // 
            // buttonCancel
            // 
            buttonCancel.Anchor = System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Right;
            buttonCancel.Location = new System.Drawing.Point(116, 86);
            buttonCancel.Name = "buttonCancel";
            buttonCancel.Size = new System.Drawing.Size(75, 23);
            buttonCancel.TabIndex = 5;
            buttonCancel.Text = "Cancel";
            buttonCancel.UseVisualStyleBackColor = true;
            // 
            // frmBattleAnimationImport
            // 
            AcceptButton = buttonOK;
            AutoScaleDimensions = new System.Drawing.SizeF(7F, 15F);
            AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            CancelButton = buttonCancel;
            ClientSize = new System.Drawing.Size(284, 121);
            ControlBox = false;
            Controls.Add(buttonCancel);
            Controls.Add(buttonOK);
            Controls.Add(rbOverwrite);
            Controls.Add(rbInsert);
            Controls.Add(label1);
            Controls.Add(cbAnimationPosition);
            FormBorderStyle = System.Windows.Forms.FormBorderStyle.FixedSingle;
            Name = "frmBattleAnimationImport";
            StartPosition = System.Windows.Forms.FormStartPosition.CenterParent;
            Text = "Import battle animation";
            ResumeLayout(false);
            PerformLayout();
        }

        #endregion

        private System.Windows.Forms.ComboBox cbAnimationPosition;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.RadioButton rbInsert;
        private System.Windows.Forms.RadioButton rbOverwrite;
        private System.Windows.Forms.Button buttonOK;
        private System.Windows.Forms.Button buttonCancel;
    }
}