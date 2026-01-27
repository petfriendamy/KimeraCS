using System;
using System.Timers;

namespace KimeraCS
{
    partial class FrmSkeletonEditor
    {
        /// <summary>
        /// Variable del diseñador necesaria.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Limpiar los recursos que se estén usando.
        /// </summary>
        /// <param name="disposing">true si los recursos administrados se deben desechar; false en caso contrario.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Código generado por el Diseñador de Windows Forms

        /// <summary>
        /// Método necesario para admitir el Diseñador. No se puede modificar
        /// el contenido de este método con el editor de código.
        /// </summary>
        private void InitializeComponent()
        {
            System.ComponentModel.ComponentResourceManager resources = new System.ComponentModel.ComponentResourceManager(typeof(FrmSkeletonEditor));
            openFile = new OpenFileDialog();
            panel1 = new Panel();
            gbTexturesFrame = new GroupBox();
            nUDMoveTextureUpDown = new NumericUpDown();
            btnRemoveTexture = new Button();
            btnAddTexture = new Button();
            btnChangeTexture = new Button();
            btnFlipHorizontal = new Button();
            btnFlipVertical = new Button();
            btnRotate = new Button();
            chkColorKeyFlag = new CheckBox();
            cbTextureSelect = new ComboBox();
            pbTextureViewer = new PictureBox();
            gbAnimationOptionsFrame = new GroupBox();
            btnInterpolateFrame = new Button();
            btnDuplicateFrame = new Button();
            btnRemoveFrame = new Button();
            chkPropagateChangesForward = new CheckBox();
            gbFrameDataPartOptions = new GroupBox();
            nUDZAnimationFramePart = new NumericUpDown();
            nUDYAnimationFramePart = new NumericUpDown();
            nUDXAnimationFramePart = new NumericUpDown();
            lblZAnimationFramePart = new Label();
            lblYAnimationFramePart = new Label();
            lblXAnimationFramePart = new Label();
            nUDFrameDataPart = new NumericUpDown();
            lblFrameOptionsPart = new Label();
            gbSelectedBoneFrame = new GroupBox();
            btnRemovePiece = new Button();
            btnAddPiece = new Button();
            nUDBoneOptionsLength = new NumericUpDown();
            lblBoneOptionsLength = new Label();
            nUDResizeBoneZ = new NumericUpDown();
            lblZScaleBoneOptions = new Label();
            nUDResizeBoneY = new NumericUpDown();
            lblYScaleBoneOptions = new Label();
            nUDResizeBoneX = new NumericUpDown();
            lblXScaleBoneOptions = new Label();
            lblBoneSelector = new Label();
            cbBoneSelector = new ComboBox();
            panel2 = new Panel();
            groupBox1 = new GroupBox();
            chkInifintyFarLights = new CheckBox();
            chkRightLight = new CheckBox();
            chkLeftLight = new CheckBox();
            chkRearLight = new CheckBox();
            chkFrontLight = new CheckBox();
            lblLightPosZ = new Label();
            hsbLightPosZ = new HScrollBar();
            lblLightPosY = new Label();
            hsbLightPosY = new HScrollBar();
            lblLightPosX = new Label();
            hsbLightPosX = new HScrollBar();
            cbWeapon = new ComboBox();
            cbBattleAnimation = new ComboBox();
            lblWeapon = new Label();
            lblBattleAnimation = new Label();
            btnComputeWeaponPosition = new Button();
            btnComputeGroundHeight = new Button();
            chkShowLastFrameGhost = new CheckBox();
            chkShowGround = new CheckBox();
            gbSelectedPieceFrame = new GroupBox();
            gbRotateFrame = new GroupBox();
            txtRotateGamma = new TextBox();
            hsbRotateGamma = new HScrollBar();
            lblRotateGamma = new Label();
            txtRotateBeta = new TextBox();
            hsbRotateBeta = new HScrollBar();
            lblRotateBeta = new Label();
            lblRotateAlpha = new Label();
            txtRotateAlpha = new TextBox();
            hsbRotateAlpha = new HScrollBar();
            gbRepositionFrame = new GroupBox();
            txtRepositionZ = new TextBox();
            hsbRepositionZ = new HScrollBar();
            lblRepositionZ = new Label();
            txtRepositionY = new TextBox();
            hsbRepositionY = new HScrollBar();
            lblRepositionY = new Label();
            lblRepositionX = new Label();
            txtRepositionX = new TextBox();
            hsbRepositionX = new HScrollBar();
            gbResizeFrame = new GroupBox();
            txtResizePieceZ = new TextBox();
            hsbResizePieceZ = new HScrollBar();
            lblResizePieceZ = new Label();
            txtResizePieceY = new TextBox();
            hsbResizePieceY = new HScrollBar();
            lblResizePieceY = new Label();
            lblResizePieceX = new Label();
            txtResizePieceX = new TextBox();
            hsbResizePieceX = new HScrollBar();
            chkDListEnable = new CheckBox();
            panel3 = new Panel();
            btnFramePrev = new Button();
            btnFrameNext = new Button();
            btnPlayStopAnim = new CheckBox();
            btnFrameEnd = new Button();
            btnFrameBegin = new Button();
            chkShowBones = new CheckBox();
            btnInterpolateAnimation = new Button();
            txtCopyPasteFrame = new TextBox();
            btnPasteFrame = new Button();
            btnCopyFrame = new Button();
            tbCurrentFrameScroll = new TrackBar();
            lblAnimationFrame = new Label();
            txtAnimationFrame = new TextBox();
            menuStrip1 = new MenuStrip();
            fileToolStripMenuItem = new ToolStripMenuItem();
            loadFieldSkeletonToolStripMenuItem = new ToolStripMenuItem();
            loadBattleMagicSkeletonToolStripMenuItem = new ToolStripMenuItem();
            loadRSDResourceToolStripMenuItem = new ToolStripMenuItem();
            loadPModelToolStripMenuItem = new ToolStripMenuItem();
            load3DSToolStripMenuItem = new ToolStripMenuItem();
            loadTMDToolStripMenuItem = new ToolStripMenuItem();
            toolStripSeparator1 = new ToolStripSeparator();
            saveSkeletonToolStripMenuItem = new ToolStripMenuItem();
            saveSkeletonAsToolStripMenuItem = new ToolStripMenuItem();
            toolStripSeparator9 = new ToolStripSeparator();
            Import3DSFixingPositionToolStripMenuItem = new ToolStripMenuItem();
            DontCheckDuplicatedPolysVertsToolStripMenuItem = new ToolStripMenuItem();
            toolStripSeparator2 = new ToolStripSeparator();
            extiToolStripMenuItem = new ToolStripMenuItem();
            editToolStripMenuItem2 = new ToolStripMenuItem();
            undoToolStripMenuItem = new ToolStripMenuItem();
            redoToolStripMenuItem = new ToolStripMenuItem();
            toolStripSeparator6 = new ToolStripSeparator();
            ShowAxesToolStripMenuItem = new ToolStripMenuItem();
            toolStripSeparator10 = new ToolStripSeparator();
            resetCameraToolStripMenuItem = new ToolStripMenuItem();
            toolStripSeparator7 = new ToolStripSeparator();
            uIOpacityToolStripMenuItem = new ToolStripMenuItem();
            tsUIOpacity100 = new ToolStripMenuItem();
            tsUIOpacity90 = new ToolStripMenuItem();
            tsUIOpacity75 = new ToolStripMenuItem();
            tsUIOpacity50 = new ToolStripMenuItem();
            tsUIOpacity25 = new ToolStripMenuItem();
            skeletonToolStripMenuItem = new ToolStripMenuItem();
            addJointToolStripMenuItem = new ToolStripMenuItem();
            editJointToolStripMenuItem = new ToolStripMenuItem();
            toolStripSeparator8 = new ToolStripSeparator();
            showNormalsToolStripMenuItem = new ToolStripMenuItem();
            showVertexNormalsToolStripMenuItem = new ToolStripMenuItem();
            showFaceNormalsToolStripMenuItem = new ToolStripMenuItem();
            normalsColorToolStripMenuItem = new ToolStripMenuItem();
            redToolStripMenuItem = new ToolStripMenuItem();
            greenToolStripMenuItem = new ToolStripMenuItem();
            blueToolStripMenuItem = new ToolStripMenuItem();
            normalsScaleToolStripMenuItem = new ToolStripMenuItem();
            oneftoolStripMenuItem = new ToolStripMenuItem();
            fiveftoolStripMenuItem = new ToolStripMenuItem();
            thirtyftoolStripMenuItem = new ToolStripMenuItem();
            thousandftoolStripMenuItem = new ToolStripMenuItem();
            toolStripSeparator11 = new ToolStripSeparator();
            statisticsToolStripMenuItem = new ToolStripMenuItem();
            textureToolStripMenuItem = new ToolStripMenuItem();
            TEXToPNGBatchConversionToolStripMenuItem = new ToolStripMenuItem();
            animationToolStripMenuItem = new ToolStripMenuItem();
            loadFieldAnimationToolStripMenuItem = new ToolStripMenuItem();
            loadBattleMagicLimitsAnimationStripMenuItem = new ToolStripMenuItem();
            toolStripSeparator3 = new ToolStripSeparator();
            saveAnimationToolStripMenuItem = new ToolStripMenuItem();
            saveAnimationAsToolStripMenuItem = new ToolStripMenuItem();
            outputFramesDataTXTToolStripMenuItem = new ToolStripMenuItem();
            inputFramesDataTXTToolStripMenuItem = new ToolStripMenuItem();
            inputFramesDataTXTToolSelectiveStripMenuItem = new ToolStripMenuItem();
            mergeFramesDataTXTToolStripMenuItem = new ToolStripMenuItem();
            toolStripSeparator4 = new ToolStripSeparator();
            toolStripMenuItem1 = new ToolStripMenuItem();
            toolStripFPS15 = new ToolStripMenuItem();
            toolStripFPS30 = new ToolStripMenuItem();
            toolStripFPS60 = new ToolStripMenuItem();
            toolStripSeparator5 = new ToolStripSeparator();
            interpolateAllAnimationsToolStripMenuItem = new ToolStripMenuItem();
            databaseToolStripMenuItem = new ToolStripMenuItem();
            showCharlgpToolStripMenuItem = new ToolStripMenuItem();
            showWorldLgpToolStripMenuItem = new ToolStripMenuItem();
            showBattlelgpToolStripMenuItem = new ToolStripMenuItem();
            showMagiclgpToolStripMenuItem = new ToolStripMenuItem();
            saveFile = new SaveFileDialog();
            panelModel = new OpenTK.GLControl.GLControl();
            panel1.SuspendLayout();
            gbTexturesFrame.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)nUDMoveTextureUpDown).BeginInit();
            ((System.ComponentModel.ISupportInitialize)pbTextureViewer).BeginInit();
            gbAnimationOptionsFrame.SuspendLayout();
            gbFrameDataPartOptions.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)nUDZAnimationFramePart).BeginInit();
            ((System.ComponentModel.ISupportInitialize)nUDYAnimationFramePart).BeginInit();
            ((System.ComponentModel.ISupportInitialize)nUDXAnimationFramePart).BeginInit();
            ((System.ComponentModel.ISupportInitialize)nUDFrameDataPart).BeginInit();
            gbSelectedBoneFrame.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)nUDBoneOptionsLength).BeginInit();
            ((System.ComponentModel.ISupportInitialize)nUDResizeBoneZ).BeginInit();
            ((System.ComponentModel.ISupportInitialize)nUDResizeBoneY).BeginInit();
            ((System.ComponentModel.ISupportInitialize)nUDResizeBoneX).BeginInit();
            panel2.SuspendLayout();
            groupBox1.SuspendLayout();
            gbSelectedPieceFrame.SuspendLayout();
            gbRotateFrame.SuspendLayout();
            gbRepositionFrame.SuspendLayout();
            gbResizeFrame.SuspendLayout();
            panel3.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)tbCurrentFrameScroll).BeginInit();
            menuStrip1.SuspendLayout();
            SuspendLayout();
            // 
            // panel1
            // 
            panel1.BackColor = SystemColors.ControlDarkDark;
            panel1.Controls.Add(gbTexturesFrame);
            panel1.Controls.Add(gbAnimationOptionsFrame);
            panel1.Controls.Add(gbSelectedBoneFrame);
            panel1.Controls.Add(lblBoneSelector);
            panel1.Controls.Add(cbBoneSelector);
            panel1.Dock = DockStyle.Right;
            panel1.Location = new Point(669, 24);
            panel1.Margin = new Padding(4, 3, 4, 3);
            panel1.Name = "panel1";
            panel1.Size = new Size(187, 725);
            panel1.TabIndex = 6;
            // 
            // gbTexturesFrame
            // 
            gbTexturesFrame.Controls.Add(nUDMoveTextureUpDown);
            gbTexturesFrame.Controls.Add(btnRemoveTexture);
            gbTexturesFrame.Controls.Add(btnAddTexture);
            gbTexturesFrame.Controls.Add(btnChangeTexture);
            gbTexturesFrame.Controls.Add(btnFlipHorizontal);
            gbTexturesFrame.Controls.Add(btnFlipVertical);
            gbTexturesFrame.Controls.Add(btnRotate);
            gbTexturesFrame.Controls.Add(chkColorKeyFlag);
            gbTexturesFrame.Controls.Add(cbTextureSelect);
            gbTexturesFrame.Controls.Add(pbTextureViewer);
            gbTexturesFrame.Enabled = false;
            gbTexturesFrame.ForeColor = SystemColors.Control;
            gbTexturesFrame.Location = new Point(5, 223);
            gbTexturesFrame.Margin = new Padding(2);
            gbTexturesFrame.Name = "gbTexturesFrame";
            gbTexturesFrame.Padding = new Padding(2);
            gbTexturesFrame.Size = new Size(177, 257);
            gbTexturesFrame.TabIndex = 15;
            gbTexturesFrame.TabStop = false;
            gbTexturesFrame.Text = "Textures (Part)";
            gbTexturesFrame.Visible = false;
            // 
            // nUDMoveTextureUpDown
            // 
            nUDMoveTextureUpDown.Location = new Point(9, 141);
            nUDMoveTextureUpDown.Margin = new Padding(2);
            nUDMoveTextureUpDown.Maximum = new decimal(new int[] { 99999999, 0, 0, 0 });
            nUDMoveTextureUpDown.Minimum = new decimal(new int[] { 99999999, 0, 0, int.MinValue });
            nUDMoveTextureUpDown.Name = "nUDMoveTextureUpDown";
            nUDMoveTextureUpDown.Size = new Size(21, 23);
            nUDMoveTextureUpDown.TabIndex = 9;
            nUDMoveTextureUpDown.ValueChanged += NudMoveTextureUpDown_ValueChanged;
            // 
            // btnRemoveTexture
            // 
            btnRemoveTexture.ForeColor = SystemColors.ActiveCaptionText;
            btnRemoveTexture.Location = new Point(4, 230);
            btnRemoveTexture.Margin = new Padding(2);
            btnRemoveTexture.Name = "btnRemoveTexture";
            btnRemoveTexture.Size = new Size(170, 24);
            btnRemoveTexture.TabIndex = 8;
            btnRemoveTexture.Text = "Remove Texture";
            btnRemoveTexture.UseVisualStyleBackColor = true;
            btnRemoveTexture.Click += BtnRemoveTexture_Click;
            // 
            // btnAddTexture
            // 
            btnAddTexture.ForeColor = SystemColors.ActiveCaptionText;
            btnAddTexture.Location = new Point(4, 207);
            btnAddTexture.Margin = new Padding(2);
            btnAddTexture.Name = "btnAddTexture";
            btnAddTexture.Size = new Size(170, 24);
            btnAddTexture.TabIndex = 7;
            btnAddTexture.Text = "Add Texture";
            btnAddTexture.UseVisualStyleBackColor = true;
            btnAddTexture.Click += BtnAddTexture_Click;
            // 
            // btnChangeTexture
            // 
            btnChangeTexture.ForeColor = SystemColors.ActiveCaptionText;
            btnChangeTexture.Location = new Point(4, 183);
            btnChangeTexture.Margin = new Padding(2);
            btnChangeTexture.Name = "btnChangeTexture";
            btnChangeTexture.Size = new Size(170, 24);
            btnChangeTexture.TabIndex = 6;
            btnChangeTexture.Text = "Change Texture";
            btnChangeTexture.UseVisualStyleBackColor = true;
            btnChangeTexture.Click += BtnChangeTexture_Click;
            // 
            // btnFlipHorizontal
            // 
            btnFlipHorizontal.BackgroundImage = (Image)resources.GetObject("btnFlipHorizontal.BackgroundImage");
            btnFlipHorizontal.BackgroundImageLayout = ImageLayout.Stretch;
            btnFlipHorizontal.Location = new Point(132, 99);
            btnFlipHorizontal.Margin = new Padding(2);
            btnFlipHorizontal.Name = "btnFlipHorizontal";
            btnFlipHorizontal.Size = new Size(37, 37);
            btnFlipHorizontal.TabIndex = 5;
            btnFlipHorizontal.UseVisualStyleBackColor = true;
            btnFlipHorizontal.Click += BtnFlipHorizontal_Click;
            // 
            // btnFlipVertical
            // 
            btnFlipVertical.BackgroundImage = (Image)resources.GetObject("btnFlipVertical.BackgroundImage");
            btnFlipVertical.BackgroundImageLayout = ImageLayout.Stretch;
            btnFlipVertical.Location = new Point(132, 60);
            btnFlipVertical.Margin = new Padding(2);
            btnFlipVertical.Name = "btnFlipVertical";
            btnFlipVertical.Size = new Size(37, 37);
            btnFlipVertical.TabIndex = 4;
            btnFlipVertical.UseVisualStyleBackColor = true;
            btnFlipVertical.Click += BtnFlipVertical_Click;
            // 
            // btnRotate
            // 
            btnRotate.BackgroundImage = (Image)resources.GetObject("btnRotate.BackgroundImage");
            btnRotate.BackgroundImageLayout = ImageLayout.Stretch;
            btnRotate.Location = new Point(132, 21);
            btnRotate.Margin = new Padding(2);
            btnRotate.Name = "btnRotate";
            btnRotate.Size = new Size(37, 37);
            btnRotate.TabIndex = 3;
            btnRotate.UseVisualStyleBackColor = true;
            btnRotate.Click += BtnRotate_Click;
            // 
            // chkColorKeyFlag
            // 
            chkColorKeyFlag.AutoSize = true;
            chkColorKeyFlag.Enabled = false;
            chkColorKeyFlag.Location = new Point(19, 166);
            chkColorKeyFlag.Margin = new Padding(2);
            chkColorKeyFlag.Name = "chkColorKeyFlag";
            chkColorKeyFlag.Size = new Size(136, 19);
            chkColorKeyFlag.TabIndex = 2;
            chkColorKeyFlag.Text = "Black is transparency";
            chkColorKeyFlag.UseVisualStyleBackColor = true;
            chkColorKeyFlag.Click += ChkZeroAsTransparent_Click;
            // 
            // cbTextureSelect
            // 
            cbTextureSelect.DropDownStyle = ComboBoxStyle.DropDownList;
            cbTextureSelect.FormattingEnabled = true;
            cbTextureSelect.Location = new Point(33, 141);
            cbTextureSelect.Margin = new Padding(2);
            cbTextureSelect.Name = "cbTextureSelect";
            cbTextureSelect.Size = new Size(134, 23);
            cbTextureSelect.TabIndex = 1;
            cbTextureSelect.SelectedIndexChanged += CbTextureSelect_SelectedIndexChanged;
            // 
            // pbTextureViewer
            // 
            pbTextureViewer.BackColor = Color.Transparent;
            pbTextureViewer.BackgroundImageLayout = ImageLayout.None;
            pbTextureViewer.BorderStyle = BorderStyle.Fixed3D;
            pbTextureViewer.Location = new Point(9, 18);
            pbTextureViewer.Margin = new Padding(2);
            pbTextureViewer.Name = "pbTextureViewer";
            pbTextureViewer.Size = new Size(121, 119);
            pbTextureViewer.TabIndex = 0;
            pbTextureViewer.TabStop = false;
            pbTextureViewer.DoubleClick += PbTextureViewer_DoubleClick;
            // 
            // gbAnimationOptionsFrame
            // 
            gbAnimationOptionsFrame.Controls.Add(btnInterpolateFrame);
            gbAnimationOptionsFrame.Controls.Add(btnDuplicateFrame);
            gbAnimationOptionsFrame.Controls.Add(btnRemoveFrame);
            gbAnimationOptionsFrame.Controls.Add(chkPropagateChangesForward);
            gbAnimationOptionsFrame.Controls.Add(gbFrameDataPartOptions);
            gbAnimationOptionsFrame.Controls.Add(nUDFrameDataPart);
            gbAnimationOptionsFrame.Controls.Add(lblFrameOptionsPart);
            gbAnimationOptionsFrame.ForeColor = SystemColors.Control;
            gbAnimationOptionsFrame.Location = new Point(7, 487);
            gbAnimationOptionsFrame.Margin = new Padding(2);
            gbAnimationOptionsFrame.Name = "gbAnimationOptionsFrame";
            gbAnimationOptionsFrame.Padding = new Padding(2);
            gbAnimationOptionsFrame.Size = new Size(177, 232);
            gbAnimationOptionsFrame.TabIndex = 14;
            gbAnimationOptionsFrame.TabStop = false;
            gbAnimationOptionsFrame.Text = "Frame options";
            gbAnimationOptionsFrame.Visible = false;
            // 
            // btnInterpolateFrame
            // 
            btnInterpolateFrame.ForeColor = SystemColors.ActiveCaptionText;
            btnInterpolateFrame.Location = new Point(4, 204);
            btnInterpolateFrame.Margin = new Padding(2);
            btnInterpolateFrame.Name = "btnInterpolateFrame";
            btnInterpolateFrame.Size = new Size(170, 24);
            btnInterpolateFrame.TabIndex = 9;
            btnInterpolateFrame.Text = "Interpolate Frame";
            btnInterpolateFrame.UseVisualStyleBackColor = true;
            btnInterpolateFrame.Click += BtnInterpolateFrame_Click;
            // 
            // btnDuplicateFrame
            // 
            btnDuplicateFrame.ForeColor = SystemColors.ActiveCaptionText;
            btnDuplicateFrame.Location = new Point(4, 181);
            btnDuplicateFrame.Margin = new Padding(2);
            btnDuplicateFrame.Name = "btnDuplicateFrame";
            btnDuplicateFrame.Size = new Size(170, 24);
            btnDuplicateFrame.TabIndex = 8;
            btnDuplicateFrame.Text = "Duplicate Frame";
            btnDuplicateFrame.UseVisualStyleBackColor = true;
            btnDuplicateFrame.Click += BtnDuplicateFrame_Click;
            // 
            // btnRemoveFrame
            // 
            btnRemoveFrame.ForeColor = SystemColors.ActiveCaptionText;
            btnRemoveFrame.Location = new Point(4, 158);
            btnRemoveFrame.Margin = new Padding(2);
            btnRemoveFrame.Name = "btnRemoveFrame";
            btnRemoveFrame.Size = new Size(170, 24);
            btnRemoveFrame.TabIndex = 7;
            btnRemoveFrame.Text = "Remove Frame";
            btnRemoveFrame.UseVisualStyleBackColor = true;
            btnRemoveFrame.Click += BtnRemoveFrame_Click;
            // 
            // chkPropagateChangesForward
            // 
            chkPropagateChangesForward.AutoSize = true;
            chkPropagateChangesForward.Checked = true;
            chkPropagateChangesForward.CheckState = CheckState.Checked;
            chkPropagateChangesForward.Location = new Point(26, 141);
            chkPropagateChangesForward.Margin = new Padding(4, 3, 4, 3);
            chkPropagateChangesForward.Name = "chkPropagateChangesForward";
            chkPropagateChangesForward.Size = new Size(124, 19);
            chkPropagateChangesForward.TabIndex = 3;
            chkPropagateChangesForward.Text = "Propagate forward";
            chkPropagateChangesForward.UseVisualStyleBackColor = true;
            // 
            // gbFrameDataPartOptions
            // 
            gbFrameDataPartOptions.BackColor = SystemColors.ControlDarkDark;
            gbFrameDataPartOptions.Controls.Add(nUDZAnimationFramePart);
            gbFrameDataPartOptions.Controls.Add(nUDYAnimationFramePart);
            gbFrameDataPartOptions.Controls.Add(nUDXAnimationFramePart);
            gbFrameDataPartOptions.Controls.Add(lblZAnimationFramePart);
            gbFrameDataPartOptions.Controls.Add(lblYAnimationFramePart);
            gbFrameDataPartOptions.Controls.Add(lblXAnimationFramePart);
            gbFrameDataPartOptions.ForeColor = SystemColors.ButtonFace;
            gbFrameDataPartOptions.Location = new Point(10, 39);
            gbFrameDataPartOptions.Margin = new Padding(2);
            gbFrameDataPartOptions.Name = "gbFrameDataPartOptions";
            gbFrameDataPartOptions.Padding = new Padding(2);
            gbFrameDataPartOptions.Size = new Size(155, 100);
            gbFrameDataPartOptions.TabIndex = 2;
            gbFrameDataPartOptions.TabStop = false;
            gbFrameDataPartOptions.Text = "Bone rotation";
            // 
            // nUDZAnimationFramePart
            // 
            nUDZAnimationFramePart.DecimalPlaces = 6;
            nUDZAnimationFramePart.Location = new Point(44, 70);
            nUDZAnimationFramePart.Margin = new Padding(2);
            nUDZAnimationFramePart.Maximum = new decimal(new int[] { 999999999, 0, 0, 0 });
            nUDZAnimationFramePart.Minimum = new decimal(new int[] { 999999999, 0, 0, int.MinValue });
            nUDZAnimationFramePart.Name = "nUDZAnimationFramePart";
            nUDZAnimationFramePart.Size = new Size(103, 23);
            nUDZAnimationFramePart.TabIndex = 10;
            nUDZAnimationFramePart.ValueChanged += NudZAnimationFramePart_ValueChanged;
            // 
            // nUDYAnimationFramePart
            // 
            nUDYAnimationFramePart.DecimalPlaces = 6;
            nUDYAnimationFramePart.Location = new Point(44, 45);
            nUDYAnimationFramePart.Margin = new Padding(2);
            nUDYAnimationFramePart.Maximum = new decimal(new int[] { 999999999, 0, 0, 0 });
            nUDYAnimationFramePart.Minimum = new decimal(new int[] { 999999999, 0, 0, int.MinValue });
            nUDYAnimationFramePart.Name = "nUDYAnimationFramePart";
            nUDYAnimationFramePart.Size = new Size(103, 23);
            nUDYAnimationFramePart.TabIndex = 9;
            nUDYAnimationFramePart.ValueChanged += NudYAnimationFramePart_ValueChanged;
            // 
            // nUDXAnimationFramePart
            // 
            nUDXAnimationFramePart.DecimalPlaces = 6;
            nUDXAnimationFramePart.Location = new Point(44, 20);
            nUDXAnimationFramePart.Margin = new Padding(2);
            nUDXAnimationFramePart.Maximum = new decimal(new int[] { 999999999, 0, 0, 0 });
            nUDXAnimationFramePart.Minimum = new decimal(new int[] { 999999999, 0, 0, int.MinValue });
            nUDXAnimationFramePart.Name = "nUDXAnimationFramePart";
            nUDXAnimationFramePart.Size = new Size(103, 23);
            nUDXAnimationFramePart.TabIndex = 8;
            nUDXAnimationFramePart.ValueChanged += NudXAnimationFramePart_ValueChanged;
            // 
            // lblZAnimationFramePart
            // 
            lblZAnimationFramePart.AutoSize = true;
            lblZAnimationFramePart.Location = new Point(19, 73);
            lblZAnimationFramePart.Margin = new Padding(2, 0, 2, 0);
            lblZAnimationFramePart.Name = "lblZAnimationFramePart";
            lblZAnimationFramePart.Size = new Size(14, 15);
            lblZAnimationFramePart.TabIndex = 7;
            lblZAnimationFramePart.Text = "Z";
            // 
            // lblYAnimationFramePart
            // 
            lblYAnimationFramePart.AutoSize = true;
            lblYAnimationFramePart.Location = new Point(19, 47);
            lblYAnimationFramePart.Margin = new Padding(2, 0, 2, 0);
            lblYAnimationFramePart.Name = "lblYAnimationFramePart";
            lblYAnimationFramePart.Size = new Size(14, 15);
            lblYAnimationFramePart.TabIndex = 6;
            lblYAnimationFramePart.Text = "Y";
            // 
            // lblXAnimationFramePart
            // 
            lblXAnimationFramePart.AutoSize = true;
            lblXAnimationFramePart.Location = new Point(19, 22);
            lblXAnimationFramePart.Margin = new Padding(2, 0, 2, 0);
            lblXAnimationFramePart.Name = "lblXAnimationFramePart";
            lblXAnimationFramePart.Size = new Size(14, 15);
            lblXAnimationFramePart.TabIndex = 5;
            lblXAnimationFramePart.Text = "X";
            // 
            // nUDFrameDataPart
            // 
            nUDFrameDataPart.Location = new Point(124, 18);
            nUDFrameDataPart.Margin = new Padding(2);
            nUDFrameDataPart.Maximum = new decimal(new int[] { 99999999, 0, 0, 0 });
            nUDFrameDataPart.Minimum = new decimal(new int[] { 99999999, 0, 0, int.MinValue });
            nUDFrameDataPart.Name = "nUDFrameDataPart";
            nUDFrameDataPart.Size = new Size(21, 23);
            nUDFrameDataPart.TabIndex = 1;
            nUDFrameDataPart.ValueChanged += NudFrameDataPart_ValueChanged;
            // 
            // lblFrameOptionsPart
            // 
            lblFrameOptionsPart.AutoSize = true;
            lblFrameOptionsPart.Location = new Point(27, 21);
            lblFrameOptionsPart.Margin = new Padding(2, 0, 2, 0);
            lblFrameOptionsPart.Name = "lblFrameOptionsPart";
            lblFrameOptionsPart.Size = new Size(90, 15);
            lblFrameOptionsPart.TabIndex = 0;
            lblFrameOptionsPart.Text = "Frame data part";
            // 
            // gbSelectedBoneFrame
            // 
            gbSelectedBoneFrame.Controls.Add(btnRemovePiece);
            gbSelectedBoneFrame.Controls.Add(btnAddPiece);
            gbSelectedBoneFrame.Controls.Add(nUDBoneOptionsLength);
            gbSelectedBoneFrame.Controls.Add(lblBoneOptionsLength);
            gbSelectedBoneFrame.Controls.Add(nUDResizeBoneZ);
            gbSelectedBoneFrame.Controls.Add(lblZScaleBoneOptions);
            gbSelectedBoneFrame.Controls.Add(nUDResizeBoneY);
            gbSelectedBoneFrame.Controls.Add(lblYScaleBoneOptions);
            gbSelectedBoneFrame.Controls.Add(nUDResizeBoneX);
            gbSelectedBoneFrame.Controls.Add(lblXScaleBoneOptions);
            gbSelectedBoneFrame.Enabled = false;
            gbSelectedBoneFrame.ForeColor = SystemColors.Control;
            gbSelectedBoneFrame.Location = new Point(5, 44);
            gbSelectedBoneFrame.Margin = new Padding(2);
            gbSelectedBoneFrame.Name = "gbSelectedBoneFrame";
            gbSelectedBoneFrame.Padding = new Padding(2);
            gbSelectedBoneFrame.Size = new Size(177, 173);
            gbSelectedBoneFrame.TabIndex = 13;
            gbSelectedBoneFrame.TabStop = false;
            gbSelectedBoneFrame.Text = "Joint options";
            gbSelectedBoneFrame.Visible = false;
            // 
            // btnRemovePiece
            // 
            btnRemovePiece.ForeColor = SystemColors.ActiveCaptionText;
            btnRemovePiece.Location = new Point(4, 145);
            btnRemovePiece.Margin = new Padding(2);
            btnRemovePiece.Name = "btnRemovePiece";
            btnRemovePiece.RightToLeft = RightToLeft.No;
            btnRemovePiece.Size = new Size(170, 24);
            btnRemovePiece.TabIndex = 9;
            btnRemovePiece.Text = "Remove part from the Joint";
            btnRemovePiece.UseVisualStyleBackColor = true;
            btnRemovePiece.Click += BtnRemovePiece_Click;
            // 
            // btnAddPiece
            // 
            btnAddPiece.ForeColor = SystemColors.ActiveCaptionText;
            btnAddPiece.Location = new Point(4, 122);
            btnAddPiece.Margin = new Padding(2);
            btnAddPiece.Name = "btnAddPiece";
            btnAddPiece.RightToLeft = RightToLeft.No;
            btnAddPiece.Size = new Size(170, 24);
            btnAddPiece.TabIndex = 8;
            btnAddPiece.Text = "Add part to the Joint";
            btnAddPiece.UseVisualStyleBackColor = true;
            btnAddPiece.Click += BtnAddPiece_Click;
            // 
            // nUDBoneOptionsLength
            // 
            nUDBoneOptionsLength.DecimalPlaces = 6;
            nUDBoneOptionsLength.Location = new Point(55, 97);
            nUDBoneOptionsLength.Margin = new Padding(2);
            nUDBoneOptionsLength.Maximum = new decimal(new int[] { 65536, 0, 0, 0 });
            nUDBoneOptionsLength.Minimum = new decimal(new int[] { 65536, 0, 0, int.MinValue });
            nUDBoneOptionsLength.Name = "nUDBoneOptionsLength";
            nUDBoneOptionsLength.Size = new Size(119, 23);
            nUDBoneOptionsLength.TabIndex = 7;
            nUDBoneOptionsLength.TextChanged += NudBoneLength_TextChanged;
            nUDBoneOptionsLength.ValueChanged += NudBoneLength_ValueChanged;
            // 
            // lblBoneOptionsLength
            // 
            lblBoneOptionsLength.AutoSize = true;
            lblBoneOptionsLength.Location = new Point(5, 99);
            lblBoneOptionsLength.Margin = new Padding(2, 0, 2, 0);
            lblBoneOptionsLength.Name = "lblBoneOptionsLength";
            lblBoneOptionsLength.Size = new Size(44, 15);
            lblBoneOptionsLength.TabIndex = 6;
            lblBoneOptionsLength.Text = "Length";
            // 
            // nUDResizeBoneZ
            // 
            nUDResizeBoneZ.Location = new Point(55, 70);
            nUDResizeBoneZ.Margin = new Padding(2);
            nUDResizeBoneZ.Maximum = new decimal(new int[] { 999999999, 0, 0, 0 });
            nUDResizeBoneZ.Minimum = new decimal(new int[] { 999999999, 0, 0, int.MinValue });
            nUDResizeBoneZ.Name = "nUDResizeBoneZ";
            nUDResizeBoneZ.Size = new Size(119, 23);
            nUDResizeBoneZ.TabIndex = 5;
            nUDResizeBoneZ.ValueChanged += NudResizeBoneZ_ValueChanged;
            // 
            // lblZScaleBoneOptions
            // 
            lblZScaleBoneOptions.AutoSize = true;
            lblZScaleBoneOptions.Location = new Point(5, 73);
            lblZScaleBoneOptions.Margin = new Padding(2, 0, 2, 0);
            lblZScaleBoneOptions.Name = "lblZScaleBoneOptions";
            lblZScaleBoneOptions.Size = new Size(44, 15);
            lblZScaleBoneOptions.TabIndex = 4;
            lblZScaleBoneOptions.Text = "Z Scale";
            // 
            // nUDResizeBoneY
            // 
            nUDResizeBoneY.Location = new Point(55, 44);
            nUDResizeBoneY.Margin = new Padding(2);
            nUDResizeBoneY.Maximum = new decimal(new int[] { 999999999, 0, 0, 0 });
            nUDResizeBoneY.Minimum = new decimal(new int[] { 999999999, 0, 0, int.MinValue });
            nUDResizeBoneY.Name = "nUDResizeBoneY";
            nUDResizeBoneY.Size = new Size(119, 23);
            nUDResizeBoneY.TabIndex = 3;
            nUDResizeBoneY.ValueChanged += NudResizeBoneY_ValueChanged;
            // 
            // lblYScaleBoneOptions
            // 
            lblYScaleBoneOptions.AutoSize = true;
            lblYScaleBoneOptions.Location = new Point(5, 46);
            lblYScaleBoneOptions.Margin = new Padding(2, 0, 2, 0);
            lblYScaleBoneOptions.Name = "lblYScaleBoneOptions";
            lblYScaleBoneOptions.Size = new Size(44, 15);
            lblYScaleBoneOptions.TabIndex = 2;
            lblYScaleBoneOptions.Text = "Y Scale";
            // 
            // nUDResizeBoneX
            // 
            nUDResizeBoneX.Location = new Point(55, 17);
            nUDResizeBoneX.Margin = new Padding(2);
            nUDResizeBoneX.Maximum = new decimal(new int[] { 999999999, 0, 0, 0 });
            nUDResizeBoneX.Minimum = new decimal(new int[] { 999999999, 0, 0, int.MinValue });
            nUDResizeBoneX.Name = "nUDResizeBoneX";
            nUDResizeBoneX.Size = new Size(119, 23);
            nUDResizeBoneX.TabIndex = 1;
            nUDResizeBoneX.ValueChanged += NudResizeBoneX_ValueChanged;
            // 
            // lblXScaleBoneOptions
            // 
            lblXScaleBoneOptions.AutoSize = true;
            lblXScaleBoneOptions.Location = new Point(5, 20);
            lblXScaleBoneOptions.Margin = new Padding(2, 0, 2, 0);
            lblXScaleBoneOptions.Name = "lblXScaleBoneOptions";
            lblXScaleBoneOptions.Size = new Size(44, 15);
            lblXScaleBoneOptions.TabIndex = 0;
            lblXScaleBoneOptions.Text = "X Scale";
            // 
            // lblBoneSelector
            // 
            lblBoneSelector.AutoSize = true;
            lblBoneSelector.ForeColor = SystemColors.Control;
            lblBoneSelector.Location = new Point(4, 2);
            lblBoneSelector.Margin = new Padding(2, 0, 2, 0);
            lblBoneSelector.Name = "lblBoneSelector";
            lblBoneSelector.Size = new Size(82, 15);
            lblBoneSelector.TabIndex = 12;
            lblBoneSelector.Text = "Selected Joint:";
            lblBoneSelector.Visible = false;
            // 
            // cbBoneSelector
            // 
            cbBoneSelector.DropDownStyle = ComboBoxStyle.DropDownList;
            cbBoneSelector.FormattingEnabled = true;
            cbBoneSelector.Location = new Point(5, 18);
            cbBoneSelector.Margin = new Padding(2);
            cbBoneSelector.Name = "cbBoneSelector";
            cbBoneSelector.Size = new Size(177, 23);
            cbBoneSelector.TabIndex = 11;
            cbBoneSelector.Visible = false;
            cbBoneSelector.SelectedIndexChanged += CbBoneSelector_SelectedIndexChanged;
            // 
            // panel2
            // 
            panel2.BackColor = SystemColors.ControlDarkDark;
            panel2.Controls.Add(groupBox1);
            panel2.Controls.Add(cbWeapon);
            panel2.Controls.Add(cbBattleAnimation);
            panel2.Controls.Add(lblWeapon);
            panel2.Controls.Add(lblBattleAnimation);
            panel2.Controls.Add(btnComputeWeaponPosition);
            panel2.Controls.Add(btnComputeGroundHeight);
            panel2.Controls.Add(chkShowLastFrameGhost);
            panel2.Controls.Add(chkShowGround);
            panel2.Controls.Add(gbSelectedPieceFrame);
            panel2.Controls.Add(chkDListEnable);
            panel2.Dock = DockStyle.Left;
            panel2.Location = new Point(0, 24);
            panel2.Margin = new Padding(4, 3, 4, 3);
            panel2.Name = "panel2";
            panel2.Size = new Size(184, 725);
            panel2.TabIndex = 7;
            // 
            // groupBox1
            // 
            groupBox1.Controls.Add(chkInifintyFarLights);
            groupBox1.Controls.Add(chkRightLight);
            groupBox1.Controls.Add(chkLeftLight);
            groupBox1.Controls.Add(chkRearLight);
            groupBox1.Controls.Add(chkFrontLight);
            groupBox1.Controls.Add(lblLightPosZ);
            groupBox1.Controls.Add(hsbLightPosZ);
            groupBox1.Controls.Add(lblLightPosY);
            groupBox1.Controls.Add(hsbLightPosY);
            groupBox1.Controls.Add(lblLightPosX);
            groupBox1.Controls.Add(hsbLightPosX);
            groupBox1.ForeColor = SystemColors.Control;
            groupBox1.Location = new Point(6, 564);
            groupBox1.Margin = new Padding(4, 3, 4, 3);
            groupBox1.Name = "groupBox1";
            groupBox1.Padding = new Padding(4, 3, 4, 3);
            groupBox1.Size = new Size(173, 155);
            groupBox1.TabIndex = 13;
            groupBox1.TabStop = false;
            groupBox1.Text = "General Lighting";
            // 
            // chkInifintyFarLights
            // 
            chkInifintyFarLights.AutoSize = true;
            chkInifintyFarLights.Location = new Point(15, 132);
            chkInifintyFarLights.Margin = new Padding(4, 3, 4, 3);
            chkInifintyFarLights.Name = "chkInifintyFarLights";
            chkInifintyFarLights.Size = new Size(124, 19);
            chkInifintyFarLights.TabIndex = 14;
            chkInifintyFarLights.Text = "Inifinitely far lights";
            chkInifintyFarLights.UseVisualStyleBackColor = true;
            chkInifintyFarLights.CheckedChanged += ChkInifintyFarLights_CheckedChanged;
            // 
            // chkRightLight
            // 
            chkRightLight.AutoSize = true;
            chkRightLight.Location = new Point(92, 112);
            chkRightLight.Margin = new Padding(4, 3, 4, 3);
            chkRightLight.Name = "chkRightLight";
            chkRightLight.Size = new Size(54, 19);
            chkRightLight.TabIndex = 13;
            chkRightLight.Text = "Right";
            chkRightLight.UseVisualStyleBackColor = true;
            chkRightLight.CheckedChanged += ChkRight_CheckedChanged;
            // 
            // chkLeftLight
            // 
            chkLeftLight.AutoSize = true;
            chkLeftLight.Location = new Point(92, 93);
            chkLeftLight.Margin = new Padding(4, 3, 4, 3);
            chkLeftLight.Name = "chkLeftLight";
            chkLeftLight.Size = new Size(46, 19);
            chkLeftLight.TabIndex = 12;
            chkLeftLight.Text = "Left";
            chkLeftLight.UseVisualStyleBackColor = true;
            chkLeftLight.CheckedChanged += ChkLeftLight_CheckedChanged;
            // 
            // chkRearLight
            // 
            chkRearLight.AutoSize = true;
            chkRearLight.Location = new Point(15, 112);
            chkRearLight.Margin = new Padding(4, 3, 4, 3);
            chkRearLight.Name = "chkRearLight";
            chkRearLight.Size = new Size(49, 19);
            chkRearLight.TabIndex = 11;
            chkRearLight.Text = "Rear";
            chkRearLight.UseVisualStyleBackColor = true;
            chkRearLight.CheckedChanged += ChkRearLight_CheckedChanged;
            // 
            // chkFrontLight
            // 
            chkFrontLight.AutoSize = true;
            chkFrontLight.Location = new Point(15, 93);
            chkFrontLight.Margin = new Padding(4, 3, 4, 3);
            chkFrontLight.Name = "chkFrontLight";
            chkFrontLight.Size = new Size(54, 19);
            chkFrontLight.TabIndex = 10;
            chkFrontLight.Text = "Front";
            chkFrontLight.UseVisualStyleBackColor = true;
            chkFrontLight.CheckedChanged += ChkFrontLight_CheckedChanged;
            // 
            // lblLightPosZ
            // 
            lblLightPosZ.AutoSize = true;
            lblLightPosZ.Location = new Point(8, 72);
            lblLightPosZ.Margin = new Padding(2, 0, 2, 0);
            lblLightPosZ.Name = "lblLightPosZ";
            lblLightPosZ.Size = new Size(14, 15);
            lblLightPosZ.TabIndex = 8;
            lblLightPosZ.Text = "Z";
            // 
            // hsbLightPosZ
            // 
            hsbLightPosZ.LargeChange = 1;
            hsbLightPosZ.Location = new Point(33, 68);
            hsbLightPosZ.Maximum = 32767;
            hsbLightPosZ.Name = "hsbLightPosZ";
            hsbLightPosZ.Size = new Size(128, 19);
            hsbLightPosZ.TabIndex = 7;
            hsbLightPosZ.ValueChanged += HsbLightPosZ_ValueChanged;
            // 
            // lblLightPosY
            // 
            lblLightPosY.AutoSize = true;
            lblLightPosY.Location = new Point(8, 48);
            lblLightPosY.Margin = new Padding(2, 0, 2, 0);
            lblLightPosY.Name = "lblLightPosY";
            lblLightPosY.Size = new Size(14, 15);
            lblLightPosY.TabIndex = 6;
            lblLightPosY.Text = "Y";
            // 
            // hsbLightPosY
            // 
            hsbLightPosY.LargeChange = 1;
            hsbLightPosY.Location = new Point(33, 45);
            hsbLightPosY.Maximum = 32767;
            hsbLightPosY.Name = "hsbLightPosY";
            hsbLightPosY.Size = new Size(128, 19);
            hsbLightPosY.TabIndex = 5;
            hsbLightPosY.ValueChanged += HsbLightPosY_ValueChanged;
            // 
            // lblLightPosX
            // 
            lblLightPosX.AutoSize = true;
            lblLightPosX.Location = new Point(8, 25);
            lblLightPosX.Margin = new Padding(2, 0, 2, 0);
            lblLightPosX.Name = "lblLightPosX";
            lblLightPosX.Size = new Size(14, 15);
            lblLightPosX.TabIndex = 4;
            lblLightPosX.Text = "X";
            // 
            // hsbLightPosX
            // 
            hsbLightPosX.LargeChange = 1;
            hsbLightPosX.Location = new Point(33, 22);
            hsbLightPosX.Maximum = 32767;
            hsbLightPosX.Name = "hsbLightPosX";
            hsbLightPosX.Size = new Size(128, 19);
            hsbLightPosX.TabIndex = 3;
            hsbLightPosX.ValueChanged += HsbLightPosX_ValueChanged;
            // 
            // cbWeapon
            // 
            cbWeapon.DropDownStyle = ComboBoxStyle.DropDownList;
            cbWeapon.FormattingEnabled = true;
            cbWeapon.Location = new Point(105, 537);
            cbWeapon.Margin = new Padding(4, 3, 4, 3);
            cbWeapon.Name = "cbWeapon";
            cbWeapon.Size = new Size(74, 23);
            cbWeapon.TabIndex = 12;
            cbWeapon.Visible = false;
            cbWeapon.SelectedIndexChanged += CbWeapon_SelectedIndexChanged;
            // 
            // cbBattleAnimation
            // 
            cbBattleAnimation.DropDownStyle = ComboBoxStyle.DropDownList;
            cbBattleAnimation.FormattingEnabled = true;
            cbBattleAnimation.Location = new Point(105, 512);
            cbBattleAnimation.Margin = new Padding(4, 3, 4, 3);
            cbBattleAnimation.Name = "cbBattleAnimation";
            cbBattleAnimation.Size = new Size(74, 23);
            cbBattleAnimation.TabIndex = 11;
            cbBattleAnimation.Visible = false;
            cbBattleAnimation.SelectedIndexChanged += CbBattleAnimation_SelectedIndexChanged;
            // 
            // lblWeapon
            // 
            lblWeapon.ForeColor = SystemColors.Control;
            lblWeapon.Location = new Point(46, 539);
            lblWeapon.Margin = new Padding(4, 0, 4, 0);
            lblWeapon.Name = "lblWeapon";
            lblWeapon.Size = new Size(59, 16);
            lblWeapon.TabIndex = 10;
            lblWeapon.Text = "Weapon:";
            lblWeapon.TextAlign = ContentAlignment.MiddleRight;
            lblWeapon.Visible = false;
            // 
            // lblBattleAnimation
            // 
            lblBattleAnimation.ForeColor = SystemColors.Control;
            lblBattleAnimation.Location = new Point(5, 515);
            lblBattleAnimation.Margin = new Padding(4, 0, 4, 0);
            lblBattleAnimation.Name = "lblBattleAnimation";
            lblBattleAnimation.RightToLeft = RightToLeft.No;
            lblBattleAnimation.Size = new Size(100, 16);
            lblBattleAnimation.TabIndex = 9;
            lblBattleAnimation.Text = "Battle Animation:";
            lblBattleAnimation.TextAlign = ContentAlignment.MiddleRight;
            lblBattleAnimation.Visible = false;
            // 
            // btnComputeWeaponPosition
            // 
            btnComputeWeaponPosition.Location = new Point(6, 470);
            btnComputeWeaponPosition.Margin = new Padding(4, 3, 4, 3);
            btnComputeWeaponPosition.Name = "btnComputeWeaponPosition";
            btnComputeWeaponPosition.Size = new Size(173, 40);
            btnComputeWeaponPosition.TabIndex = 8;
            btnComputeWeaponPosition.Text = "Compute attached weapon position";
            btnComputeWeaponPosition.UseVisualStyleBackColor = true;
            btnComputeWeaponPosition.Visible = false;
            btnComputeWeaponPosition.Click += BtnComputeWeaponPosition_Click;
            // 
            // btnComputeGroundHeight
            // 
            btnComputeGroundHeight.Location = new Point(6, 445);
            btnComputeGroundHeight.Margin = new Padding(4, 3, 4, 3);
            btnComputeGroundHeight.Name = "btnComputeGroundHeight";
            btnComputeGroundHeight.Size = new Size(173, 25);
            btnComputeGroundHeight.TabIndex = 7;
            btnComputeGroundHeight.Text = "Compute ground height";
            btnComputeGroundHeight.UseVisualStyleBackColor = true;
            btnComputeGroundHeight.Visible = false;
            btnComputeGroundHeight.Click += BtnComputeGroundHeight_Click;
            // 
            // chkShowLastFrameGhost
            // 
            chkShowLastFrameGhost.AutoSize = true;
            chkShowLastFrameGhost.ForeColor = SystemColors.Control;
            chkShowLastFrameGhost.Location = new Point(12, 390);
            chkShowLastFrameGhost.Margin = new Padding(2);
            chkShowLastFrameGhost.Name = "chkShowLastFrameGhost";
            chkShowLastFrameGhost.Size = new Size(155, 19);
            chkShowLastFrameGhost.TabIndex = 6;
            chkShowLastFrameGhost.Text = "Overlap last frame ghost";
            chkShowLastFrameGhost.UseVisualStyleBackColor = true;
            chkShowLastFrameGhost.CheckedChanged += ChkShowLastFrameGhost_CheckedChanged;
            // 
            // chkShowGround
            // 
            chkShowGround.AutoSize = true;
            chkShowGround.ForeColor = SystemColors.Control;
            chkShowGround.Location = new Point(12, 407);
            chkShowGround.Margin = new Padding(2);
            chkShowGround.Name = "chkShowGround";
            chkShowGround.Size = new Size(97, 19);
            chkShowGround.TabIndex = 5;
            chkShowGround.Text = "Show ground";
            chkShowGround.UseVisualStyleBackColor = true;
            chkShowGround.CheckedChanged += ChkShowGround_CheckedChanged;
            // 
            // gbSelectedPieceFrame
            // 
            gbSelectedPieceFrame.Controls.Add(gbRotateFrame);
            gbSelectedPieceFrame.Controls.Add(gbRepositionFrame);
            gbSelectedPieceFrame.Controls.Add(gbResizeFrame);
            gbSelectedPieceFrame.Enabled = false;
            gbSelectedPieceFrame.ForeColor = SystemColors.Control;
            gbSelectedPieceFrame.Location = new Point(6, 3);
            gbSelectedPieceFrame.Margin = new Padding(2);
            gbSelectedPieceFrame.Name = "gbSelectedPieceFrame";
            gbSelectedPieceFrame.Padding = new Padding(2);
            gbSelectedPieceFrame.Size = new Size(173, 384);
            gbSelectedPieceFrame.TabIndex = 4;
            gbSelectedPieceFrame.TabStop = false;
            gbSelectedPieceFrame.Text = "Selected Piece";
            // 
            // gbRotateFrame
            // 
            gbRotateFrame.Controls.Add(txtRotateGamma);
            gbRotateFrame.Controls.Add(hsbRotateGamma);
            gbRotateFrame.Controls.Add(lblRotateGamma);
            gbRotateFrame.Controls.Add(txtRotateBeta);
            gbRotateFrame.Controls.Add(hsbRotateBeta);
            gbRotateFrame.Controls.Add(lblRotateBeta);
            gbRotateFrame.Controls.Add(lblRotateAlpha);
            gbRotateFrame.Controls.Add(txtRotateAlpha);
            gbRotateFrame.Controls.Add(hsbRotateAlpha);
            gbRotateFrame.ForeColor = SystemColors.ButtonFace;
            gbRotateFrame.Location = new Point(6, 222);
            gbRotateFrame.Margin = new Padding(4, 3, 4, 3);
            gbRotateFrame.Name = "gbRotateFrame";
            gbRotateFrame.Padding = new Padding(4, 3, 4, 3);
            gbRotateFrame.Size = new Size(161, 158);
            gbRotateFrame.TabIndex = 2;
            gbRotateFrame.TabStop = false;
            gbRotateFrame.Text = "Rotation";
            // 
            // txtRotateGamma
            // 
            txtRotateGamma.Location = new Point(118, 128);
            txtRotateGamma.Margin = new Padding(4, 3, 4, 3);
            txtRotateGamma.Name = "txtRotateGamma";
            txtRotateGamma.Size = new Size(32, 23);
            txtRotateGamma.TabIndex = 8;
            txtRotateGamma.Text = "0";
            txtRotateGamma.TextChanged += TxtRotateGamma_TextChanged;
            // 
            // hsbRotateGamma
            // 
            hsbRotateGamma.LargeChange = 1;
            hsbRotateGamma.Location = new Point(8, 128);
            hsbRotateGamma.Maximum = 360;
            hsbRotateGamma.Name = "hsbRotateGamma";
            hsbRotateGamma.Size = new Size(105, 19);
            hsbRotateGamma.TabIndex = 7;
            hsbRotateGamma.ValueChanged += HsbRotateGamma_ValueChanged;
            // 
            // lblRotateGamma
            // 
            lblRotateGamma.AutoSize = true;
            lblRotateGamma.Location = new Point(6, 111);
            lblRotateGamma.Margin = new Padding(2, 0, 2, 0);
            lblRotateGamma.Name = "lblRotateGamma";
            lblRotateGamma.Size = new Size(139, 15);
            lblRotateGamma.TabIndex = 6;
            lblRotateGamma.Text = "Gamma rotation (Z-Axis)";
            // 
            // txtRotateBeta
            // 
            txtRotateBeta.Location = new Point(119, 83);
            txtRotateBeta.Margin = new Padding(4, 3, 4, 3);
            txtRotateBeta.Name = "txtRotateBeta";
            txtRotateBeta.Size = new Size(32, 23);
            txtRotateBeta.TabIndex = 5;
            txtRotateBeta.Text = "0";
            txtRotateBeta.TextChanged += TxtRotateBeta_TextChanged;
            // 
            // hsbRotateBeta
            // 
            hsbRotateBeta.LargeChange = 1;
            hsbRotateBeta.Location = new Point(9, 83);
            hsbRotateBeta.Maximum = 360;
            hsbRotateBeta.Name = "hsbRotateBeta";
            hsbRotateBeta.Size = new Size(105, 19);
            hsbRotateBeta.TabIndex = 4;
            hsbRotateBeta.ValueChanged += HsbRotateBeta_ValueChanged;
            // 
            // lblRotateBeta
            // 
            lblRotateBeta.AutoSize = true;
            lblRotateBeta.Location = new Point(6, 66);
            lblRotateBeta.Margin = new Padding(2, 0, 2, 0);
            lblRotateBeta.Name = "lblRotateBeta";
            lblRotateBeta.Size = new Size(120, 15);
            lblRotateBeta.TabIndex = 3;
            lblRotateBeta.Text = "Beta rotation (Y-Axis)";
            // 
            // lblRotateAlpha
            // 
            lblRotateAlpha.AutoSize = true;
            lblRotateAlpha.Location = new Point(6, 21);
            lblRotateAlpha.Margin = new Padding(2, 0, 2, 0);
            lblRotateAlpha.Name = "lblRotateAlpha";
            lblRotateAlpha.Size = new Size(128, 15);
            lblRotateAlpha.TabIndex = 2;
            lblRotateAlpha.Text = "Alpha rotation (X-Axis)";
            // 
            // txtRotateAlpha
            // 
            txtRotateAlpha.Location = new Point(119, 38);
            txtRotateAlpha.Margin = new Padding(4, 3, 4, 3);
            txtRotateAlpha.Name = "txtRotateAlpha";
            txtRotateAlpha.Size = new Size(32, 23);
            txtRotateAlpha.TabIndex = 1;
            txtRotateAlpha.Text = "0";
            txtRotateAlpha.TextChanged += TxtRotateAlpha_TextChanged;
            // 
            // hsbRotateAlpha
            // 
            hsbRotateAlpha.LargeChange = 1;
            hsbRotateAlpha.Location = new Point(9, 38);
            hsbRotateAlpha.Maximum = 360;
            hsbRotateAlpha.Name = "hsbRotateAlpha";
            hsbRotateAlpha.Size = new Size(105, 19);
            hsbRotateAlpha.TabIndex = 0;
            hsbRotateAlpha.ValueChanged += HsbRotateAlpha_ValueChanged;
            // 
            // gbRepositionFrame
            // 
            gbRepositionFrame.Controls.Add(txtRepositionZ);
            gbRepositionFrame.Controls.Add(hsbRepositionZ);
            gbRepositionFrame.Controls.Add(lblRepositionZ);
            gbRepositionFrame.Controls.Add(txtRepositionY);
            gbRepositionFrame.Controls.Add(hsbRepositionY);
            gbRepositionFrame.Controls.Add(lblRepositionY);
            gbRepositionFrame.Controls.Add(lblRepositionX);
            gbRepositionFrame.Controls.Add(txtRepositionX);
            gbRepositionFrame.Controls.Add(hsbRepositionX);
            gbRepositionFrame.ForeColor = SystemColors.ButtonFace;
            gbRepositionFrame.Location = new Point(6, 119);
            gbRepositionFrame.Margin = new Padding(4, 3, 4, 3);
            gbRepositionFrame.Name = "gbRepositionFrame";
            gbRepositionFrame.Padding = new Padding(4, 3, 4, 3);
            gbRepositionFrame.Size = new Size(161, 99);
            gbRepositionFrame.TabIndex = 1;
            gbRepositionFrame.TabStop = false;
            gbRepositionFrame.Text = "Reposition";
            // 
            // txtRepositionZ
            // 
            txtRepositionZ.Location = new Point(121, 69);
            txtRepositionZ.Margin = new Padding(4, 3, 4, 3);
            txtRepositionZ.Name = "txtRepositionZ";
            txtRepositionZ.Size = new Size(32, 23);
            txtRepositionZ.TabIndex = 8;
            txtRepositionZ.Text = "0";
            txtRepositionZ.TextChanged += TxtRepositionZ_TextChanged;
            // 
            // hsbRepositionZ
            // 
            hsbRepositionZ.LargeChange = 1;
            hsbRepositionZ.Location = new Point(26, 69);
            hsbRepositionZ.Maximum = 500;
            hsbRepositionZ.Minimum = -500;
            hsbRepositionZ.Name = "hsbRepositionZ";
            hsbRepositionZ.Size = new Size(91, 19);
            hsbRepositionZ.TabIndex = 7;
            hsbRepositionZ.ValueChanged += HsbRepositionZ_ValueChanged;
            // 
            // lblRepositionZ
            // 
            lblRepositionZ.AutoSize = true;
            lblRepositionZ.Location = new Point(7, 70);
            lblRepositionZ.Margin = new Padding(2, 0, 2, 0);
            lblRepositionZ.Name = "lblRepositionZ";
            lblRepositionZ.Size = new Size(14, 15);
            lblRepositionZ.TabIndex = 6;
            lblRepositionZ.Text = "Z";
            // 
            // txtRepositionY
            // 
            txtRepositionY.Location = new Point(121, 46);
            txtRepositionY.Margin = new Padding(4, 3, 4, 3);
            txtRepositionY.Name = "txtRepositionY";
            txtRepositionY.Size = new Size(32, 23);
            txtRepositionY.TabIndex = 5;
            txtRepositionY.Text = "0";
            txtRepositionY.TextChanged += TxtRepositionY_TextChanged;
            // 
            // hsbRepositionY
            // 
            hsbRepositionY.LargeChange = 1;
            hsbRepositionY.Location = new Point(26, 46);
            hsbRepositionY.Maximum = 500;
            hsbRepositionY.Minimum = -500;
            hsbRepositionY.Name = "hsbRepositionY";
            hsbRepositionY.Size = new Size(91, 19);
            hsbRepositionY.TabIndex = 4;
            hsbRepositionY.ValueChanged += HsbRepositionY_ValueChanged;
            // 
            // lblRepositionY
            // 
            lblRepositionY.AutoSize = true;
            lblRepositionY.Location = new Point(7, 48);
            lblRepositionY.Margin = new Padding(2, 0, 2, 0);
            lblRepositionY.Name = "lblRepositionY";
            lblRepositionY.Size = new Size(14, 15);
            lblRepositionY.TabIndex = 3;
            lblRepositionY.Text = "Y";
            // 
            // lblRepositionX
            // 
            lblRepositionX.AutoSize = true;
            lblRepositionX.Location = new Point(7, 27);
            lblRepositionX.Margin = new Padding(2, 0, 2, 0);
            lblRepositionX.Name = "lblRepositionX";
            lblRepositionX.Size = new Size(14, 15);
            lblRepositionX.TabIndex = 2;
            lblRepositionX.Text = "X";
            // 
            // txtRepositionX
            // 
            txtRepositionX.Location = new Point(121, 23);
            txtRepositionX.Margin = new Padding(4, 3, 4, 3);
            txtRepositionX.Name = "txtRepositionX";
            txtRepositionX.Size = new Size(32, 23);
            txtRepositionX.TabIndex = 1;
            txtRepositionX.Text = "0";
            txtRepositionX.TextChanged += TxtRepositionX_TextChanged;
            // 
            // hsbRepositionX
            // 
            hsbRepositionX.LargeChange = 1;
            hsbRepositionX.Location = new Point(26, 23);
            hsbRepositionX.Maximum = 500;
            hsbRepositionX.Minimum = -500;
            hsbRepositionX.Name = "hsbRepositionX";
            hsbRepositionX.Size = new Size(91, 19);
            hsbRepositionX.TabIndex = 0;
            hsbRepositionX.ValueChanged += HsbRepositionX_ValueChanged;
            // 
            // gbResizeFrame
            // 
            gbResizeFrame.Controls.Add(txtResizePieceZ);
            gbResizeFrame.Controls.Add(hsbResizePieceZ);
            gbResizeFrame.Controls.Add(lblResizePieceZ);
            gbResizeFrame.Controls.Add(txtResizePieceY);
            gbResizeFrame.Controls.Add(hsbResizePieceY);
            gbResizeFrame.Controls.Add(lblResizePieceY);
            gbResizeFrame.Controls.Add(lblResizePieceX);
            gbResizeFrame.Controls.Add(txtResizePieceX);
            gbResizeFrame.Controls.Add(hsbResizePieceX);
            gbResizeFrame.ForeColor = SystemColors.ButtonFace;
            gbResizeFrame.Location = new Point(6, 16);
            gbResizeFrame.Margin = new Padding(4, 3, 4, 3);
            gbResizeFrame.Name = "gbResizeFrame";
            gbResizeFrame.Padding = new Padding(4, 3, 4, 3);
            gbResizeFrame.Size = new Size(161, 99);
            gbResizeFrame.TabIndex = 0;
            gbResizeFrame.TabStop = false;
            gbResizeFrame.Text = "Resize";
            // 
            // txtResizePieceZ
            // 
            txtResizePieceZ.Location = new Point(121, 69);
            txtResizePieceZ.Margin = new Padding(4, 3, 4, 3);
            txtResizePieceZ.Name = "txtResizePieceZ";
            txtResizePieceZ.Size = new Size(32, 23);
            txtResizePieceZ.TabIndex = 8;
            txtResizePieceZ.Text = "100";
            txtResizePieceZ.TextChanged += TxtResizePieceZ_TextChanged;
            // 
            // hsbResizePieceZ
            // 
            hsbResizePieceZ.LargeChange = 1;
            hsbResizePieceZ.Location = new Point(26, 69);
            hsbResizePieceZ.Maximum = 500;
            hsbResizePieceZ.Name = "hsbResizePieceZ";
            hsbResizePieceZ.Size = new Size(91, 19);
            hsbResizePieceZ.TabIndex = 7;
            hsbResizePieceZ.Value = 100;
            hsbResizePieceZ.ValueChanged += HsbResizePieceZ_ValueChanged;
            // 
            // lblResizePieceZ
            // 
            lblResizePieceZ.AutoSize = true;
            lblResizePieceZ.Location = new Point(7, 70);
            lblResizePieceZ.Margin = new Padding(2, 0, 2, 0);
            lblResizePieceZ.Name = "lblResizePieceZ";
            lblResizePieceZ.Size = new Size(14, 15);
            lblResizePieceZ.TabIndex = 6;
            lblResizePieceZ.Text = "Z";
            // 
            // txtResizePieceY
            // 
            txtResizePieceY.Location = new Point(121, 46);
            txtResizePieceY.Margin = new Padding(4, 3, 4, 3);
            txtResizePieceY.Name = "txtResizePieceY";
            txtResizePieceY.Size = new Size(32, 23);
            txtResizePieceY.TabIndex = 5;
            txtResizePieceY.Text = "100";
            txtResizePieceY.TextChanged += TxtResizePieceY_TextChanged;
            // 
            // hsbResizePieceY
            // 
            hsbResizePieceY.LargeChange = 1;
            hsbResizePieceY.Location = new Point(26, 46);
            hsbResizePieceY.Maximum = 500;
            hsbResizePieceY.Name = "hsbResizePieceY";
            hsbResizePieceY.Size = new Size(91, 19);
            hsbResizePieceY.TabIndex = 4;
            hsbResizePieceY.Value = 100;
            hsbResizePieceY.ValueChanged += HsbResizePieceY_ValueChanged;
            // 
            // lblResizePieceY
            // 
            lblResizePieceY.AutoSize = true;
            lblResizePieceY.Location = new Point(7, 48);
            lblResizePieceY.Margin = new Padding(2, 0, 2, 0);
            lblResizePieceY.Name = "lblResizePieceY";
            lblResizePieceY.Size = new Size(14, 15);
            lblResizePieceY.TabIndex = 3;
            lblResizePieceY.Text = "Y";
            // 
            // lblResizePieceX
            // 
            lblResizePieceX.AutoSize = true;
            lblResizePieceX.Location = new Point(7, 27);
            lblResizePieceX.Margin = new Padding(2, 0, 2, 0);
            lblResizePieceX.Name = "lblResizePieceX";
            lblResizePieceX.Size = new Size(14, 15);
            lblResizePieceX.TabIndex = 2;
            lblResizePieceX.Text = "X";
            // 
            // txtResizePieceX
            // 
            txtResizePieceX.Location = new Point(121, 23);
            txtResizePieceX.Margin = new Padding(4, 3, 4, 3);
            txtResizePieceX.Name = "txtResizePieceX";
            txtResizePieceX.Size = new Size(32, 23);
            txtResizePieceX.TabIndex = 1;
            txtResizePieceX.Text = "100";
            txtResizePieceX.TextChanged += TxtResizePieceX_TextChanged;
            // 
            // hsbResizePieceX
            // 
            hsbResizePieceX.LargeChange = 1;
            hsbResizePieceX.Location = new Point(26, 23);
            hsbResizePieceX.Maximum = 500;
            hsbResizePieceX.Name = "hsbResizePieceX";
            hsbResizePieceX.Size = new Size(91, 19);
            hsbResizePieceX.TabIndex = 0;
            hsbResizePieceX.Value = 100;
            hsbResizePieceX.ValueChanged += HsbResizePieceX_ValueChanged;
            // 
            // chkDListEnable
            // 
            chkDListEnable.AutoSize = true;
            chkDListEnable.ForeColor = SystemColors.Control;
            chkDListEnable.Location = new Point(12, 426);
            chkDListEnable.Margin = new Padding(2);
            chkDListEnable.Name = "chkDListEnable";
            chkDListEnable.Size = new Size(129, 19);
            chkDListEnable.TabIndex = 3;
            chkDListEnable.Text = "Render using DLists";
            chkDListEnable.UseVisualStyleBackColor = true;
            chkDListEnable.CheckedChanged += CheckBDListEnable_CheckedChanged;
            // 
            // panel3
            // 
            panel3.BackColor = SystemColors.WindowFrame;
            panel3.Controls.Add(btnFramePrev);
            panel3.Controls.Add(btnFrameNext);
            panel3.Controls.Add(btnPlayStopAnim);
            panel3.Controls.Add(btnFrameEnd);
            panel3.Controls.Add(btnFrameBegin);
            panel3.Controls.Add(chkShowBones);
            panel3.Controls.Add(btnInterpolateAnimation);
            panel3.Controls.Add(txtCopyPasteFrame);
            panel3.Controls.Add(btnPasteFrame);
            panel3.Controls.Add(btnCopyFrame);
            panel3.Controls.Add(tbCurrentFrameScroll);
            panel3.Controls.Add(lblAnimationFrame);
            panel3.Controls.Add(txtAnimationFrame);
            panel3.Dock = DockStyle.Bottom;
            panel3.Location = new Point(184, 667);
            panel3.Margin = new Padding(4, 3, 4, 3);
            panel3.Name = "panel3";
            panel3.Size = new Size(485, 82);
            panel3.TabIndex = 8;
            // 
            // btnFramePrev
            // 
            btnFramePrev.BackgroundImage = (Image)resources.GetObject("btnFramePrev.BackgroundImage");
            btnFramePrev.BackgroundImageLayout = ImageLayout.Stretch;
            btnFramePrev.Location = new Point(334, 8);
            btnFramePrev.Margin = new Padding(4, 3, 4, 3);
            btnFramePrev.Name = "btnFramePrev";
            btnFramePrev.Size = new Size(37, 37);
            btnFramePrev.TabIndex = 25;
            btnFramePrev.UseVisualStyleBackColor = true;
            btnFramePrev.Visible = false;
            btnFramePrev.Click += BtnFramePrev_Click;
            // 
            // btnFrameNext
            // 
            btnFrameNext.BackgroundImage = (Image)resources.GetObject("btnFrameNext.BackgroundImage");
            btnFrameNext.BackgroundImageLayout = ImageLayout.Stretch;
            btnFrameNext.Location = new Point(406, 8);
            btnFrameNext.Margin = new Padding(4, 3, 4, 3);
            btnFrameNext.Name = "btnFrameNext";
            btnFrameNext.Size = new Size(37, 37);
            btnFrameNext.TabIndex = 24;
            btnFrameNext.UseVisualStyleBackColor = true;
            btnFrameNext.Visible = false;
            btnFrameNext.Click += BtnFrameNext_Click;
            // 
            // btnPlayStopAnim
            // 
            btnPlayStopAnim.Appearance = Appearance.Button;
            btnPlayStopAnim.BackgroundImage = (Image)resources.GetObject("btnPlayStopAnim.BackgroundImage");
            btnPlayStopAnim.BackgroundImageLayout = ImageLayout.Stretch;
            btnPlayStopAnim.Location = new Point(370, 8);
            btnPlayStopAnim.Margin = new Padding(4, 3, 4, 3);
            btnPlayStopAnim.Name = "btnPlayStopAnim";
            btnPlayStopAnim.Size = new Size(37, 37);
            btnPlayStopAnim.TabIndex = 23;
            btnPlayStopAnim.UseVisualStyleBackColor = true;
            btnPlayStopAnim.Visible = false;
            btnPlayStopAnim.CheckedChanged += BtnPlayStopAnm_CheckedChanged;
            // 
            // btnFrameEnd
            // 
            btnFrameEnd.BackgroundImage = (Image)resources.GetObject("btnFrameEnd.BackgroundImage");
            btnFrameEnd.BackgroundImageLayout = ImageLayout.Stretch;
            btnFrameEnd.Location = new Point(442, 8);
            btnFrameEnd.Margin = new Padding(4, 3, 4, 3);
            btnFrameEnd.Name = "btnFrameEnd";
            btnFrameEnd.Size = new Size(37, 37);
            btnFrameEnd.TabIndex = 22;
            btnFrameEnd.UseVisualStyleBackColor = true;
            btnFrameEnd.Visible = false;
            btnFrameEnd.Click += BtnFrameEnd_Click;
            // 
            // btnFrameBegin
            // 
            btnFrameBegin.BackgroundImage = (Image)resources.GetObject("btnFrameBegin.BackgroundImage");
            btnFrameBegin.BackgroundImageLayout = ImageLayout.Stretch;
            btnFrameBegin.Location = new Point(298, 8);
            btnFrameBegin.Margin = new Padding(4, 3, 4, 3);
            btnFrameBegin.Name = "btnFrameBegin";
            btnFrameBegin.Size = new Size(37, 37);
            btnFrameBegin.TabIndex = 21;
            btnFrameBegin.UseVisualStyleBackColor = true;
            btnFrameBegin.Visible = false;
            btnFrameBegin.Click += BtnFrameBegin_Click;
            // 
            // chkShowBones
            // 
            chkShowBones.AutoSize = true;
            chkShowBones.ForeColor = SystemColors.Control;
            chkShowBones.Location = new Point(224, 59);
            chkShowBones.Margin = new Padding(2);
            chkShowBones.Name = "chkShowBones";
            chkShowBones.Size = new Size(90, 19);
            chkShowBones.TabIndex = 19;
            chkShowBones.Text = "Show Bones";
            chkShowBones.UseVisualStyleBackColor = true;
            chkShowBones.Visible = false;
            chkShowBones.CheckedChanged += ChkShowBones_CheckedChanged;
            // 
            // btnInterpolateAnimation
            // 
            btnInterpolateAnimation.Location = new Point(336, 54);
            btnInterpolateAnimation.Margin = new Padding(2);
            btnInterpolateAnimation.Name = "btnInterpolateAnimation";
            btnInterpolateAnimation.Size = new Size(141, 24);
            btnInterpolateAnimation.TabIndex = 18;
            btnInterpolateAnimation.Text = "Interpolate Animation";
            btnInterpolateAnimation.UseVisualStyleBackColor = true;
            btnInterpolateAnimation.Visible = false;
            btnInterpolateAnimation.Click += BtnInterpolateAnimation_Click;
            // 
            // txtCopyPasteFrame
            // 
            txtCopyPasteFrame.BackColor = SystemColors.ControlLight;
            txtCopyPasteFrame.CausesValidation = false;
            txtCopyPasteFrame.Location = new Point(82, 54);
            txtCopyPasteFrame.Margin = new Padding(2);
            txtCopyPasteFrame.Name = "txtCopyPasteFrame";
            txtCopyPasteFrame.ReadOnly = true;
            txtCopyPasteFrame.ShortcutsEnabled = false;
            txtCopyPasteFrame.Size = new Size(122, 23);
            txtCopyPasteFrame.TabIndex = 17;
            txtCopyPasteFrame.TabStop = false;
            txtCopyPasteFrame.Visible = false;
            // 
            // btnPasteFrame
            // 
            btnPasteFrame.BackgroundImage = (Image)resources.GetObject("btnPasteFrame.BackgroundImage");
            btnPasteFrame.BackgroundImageLayout = ImageLayout.Stretch;
            btnPasteFrame.Enabled = false;
            btnPasteFrame.Location = new Point(44, 42);
            btnPasteFrame.Margin = new Padding(2);
            btnPasteFrame.Name = "btnPasteFrame";
            btnPasteFrame.Size = new Size(37, 37);
            btnPasteFrame.TabIndex = 16;
            btnPasteFrame.UseVisualStyleBackColor = true;
            btnPasteFrame.Visible = false;
            btnPasteFrame.Click += BtnPasteFrame_Click;
            // 
            // btnCopyFrame
            // 
            btnCopyFrame.BackgroundImage = (Image)resources.GetObject("btnCopyFrame.BackgroundImage");
            btnCopyFrame.BackgroundImageLayout = ImageLayout.Stretch;
            btnCopyFrame.Location = new Point(8, 42);
            btnCopyFrame.Margin = new Padding(2);
            btnCopyFrame.Name = "btnCopyFrame";
            btnCopyFrame.Size = new Size(37, 37);
            btnCopyFrame.TabIndex = 15;
            btnCopyFrame.UseVisualStyleBackColor = true;
            btnCopyFrame.Visible = false;
            btnCopyFrame.Click += BtnCopyFrame_Click;
            // 
            // tbCurrentFrameScroll
            // 
            tbCurrentFrameScroll.Enabled = false;
            tbCurrentFrameScroll.Location = new Point(106, 2);
            tbCurrentFrameScroll.Margin = new Padding(2);
            tbCurrentFrameScroll.Maximum = 5;
            tbCurrentFrameScroll.Name = "tbCurrentFrameScroll";
            tbCurrentFrameScroll.Size = new Size(187, 45);
            tbCurrentFrameScroll.TabIndex = 14;
            tbCurrentFrameScroll.Visible = false;
            tbCurrentFrameScroll.Scroll += TbCurrentFrameScroll_Scroll;
            tbCurrentFrameScroll.ValueChanged += TbCurrentFrameScroll_ValueChanged;
            // 
            // lblAnimationFrame
            // 
            lblAnimationFrame.AutoSize = true;
            lblAnimationFrame.ForeColor = SystemColors.Control;
            lblAnimationFrame.Location = new Point(8, 14);
            lblAnimationFrame.Margin = new Padding(2, 0, 2, 0);
            lblAnimationFrame.Name = "lblAnimationFrame";
            lblAnimationFrame.Size = new Size(43, 15);
            lblAnimationFrame.TabIndex = 13;
            lblAnimationFrame.Text = "Frame:";
            lblAnimationFrame.Visible = false;
            // 
            // txtAnimationFrame
            // 
            txtAnimationFrame.Location = new Point(54, 9);
            txtAnimationFrame.Margin = new Padding(2);
            txtAnimationFrame.Name = "txtAnimationFrame";
            txtAnimationFrame.ReadOnly = true;
            txtAnimationFrame.Size = new Size(48, 23);
            txtAnimationFrame.TabIndex = 1;
            txtAnimationFrame.TextAlign = HorizontalAlignment.Center;
            txtAnimationFrame.Visible = false;
            // 
            // menuStrip1
            // 
            menuStrip1.ImageScalingSize = new Size(20, 20);
            menuStrip1.Items.AddRange(new ToolStripItem[] { fileToolStripMenuItem, editToolStripMenuItem2, skeletonToolStripMenuItem, textureToolStripMenuItem, animationToolStripMenuItem, databaseToolStripMenuItem });
            menuStrip1.Location = new Point(0, 0);
            menuStrip1.Name = "menuStrip1";
            menuStrip1.Padding = new Padding(5, 2, 0, 2);
            menuStrip1.Size = new Size(856, 24);
            menuStrip1.TabIndex = 10;
            menuStrip1.Text = "menuStrip1";
            // 
            // fileToolStripMenuItem
            // 
            fileToolStripMenuItem.DropDownItems.AddRange(new ToolStripItem[] { loadFieldSkeletonToolStripMenuItem, loadBattleMagicSkeletonToolStripMenuItem, loadRSDResourceToolStripMenuItem, loadPModelToolStripMenuItem, load3DSToolStripMenuItem, loadTMDToolStripMenuItem, toolStripSeparator1, saveSkeletonToolStripMenuItem, saveSkeletonAsToolStripMenuItem, toolStripSeparator9, Import3DSFixingPositionToolStripMenuItem, DontCheckDuplicatedPolysVertsToolStripMenuItem, toolStripSeparator2, extiToolStripMenuItem });
            fileToolStripMenuItem.Name = "fileToolStripMenuItem";
            fileToolStripMenuItem.Size = new Size(37, 20);
            fileToolStripMenuItem.Text = "&File";
            // 
            // loadFieldSkeletonToolStripMenuItem
            // 
            loadFieldSkeletonToolStripMenuItem.Name = "loadFieldSkeletonToolStripMenuItem";
            loadFieldSkeletonToolStripMenuItem.Size = new Size(292, 22);
            loadFieldSkeletonToolStripMenuItem.Text = "Load Field/World Skeleton";
            loadFieldSkeletonToolStripMenuItem.Click += LoadFieldSkeletonToolStripMenuItem_Click;
            // 
            // loadBattleMagicSkeletonToolStripMenuItem
            // 
            loadBattleMagicSkeletonToolStripMenuItem.Name = "loadBattleMagicSkeletonToolStripMenuItem";
            loadBattleMagicSkeletonToolStripMenuItem.Size = new Size(292, 22);
            loadBattleMagicSkeletonToolStripMenuItem.Text = "Load Battle/Magic Skeleton";
            loadBattleMagicSkeletonToolStripMenuItem.Click += LoadBattleMagicSkeletonToolStripMenuItem_Click;
            // 
            // loadRSDResourceToolStripMenuItem
            // 
            loadRSDResourceToolStripMenuItem.Name = "loadRSDResourceToolStripMenuItem";
            loadRSDResourceToolStripMenuItem.Size = new Size(292, 22);
            loadRSDResourceToolStripMenuItem.Text = "Load RSD Resource";
            loadRSDResourceToolStripMenuItem.Click += LoadRSDResourceToolStripMenuItem_Click;
            // 
            // loadPModelToolStripMenuItem
            // 
            loadPModelToolStripMenuItem.Name = "loadPModelToolStripMenuItem";
            loadPModelToolStripMenuItem.Size = new Size(292, 22);
            loadPModelToolStripMenuItem.Text = "Load P Model";
            loadPModelToolStripMenuItem.Click += LoadPModelToolStripMenuItem_Click;
            // 
            // load3DSToolStripMenuItem
            // 
            load3DSToolStripMenuItem.Name = "load3DSToolStripMenuItem";
            load3DSToolStripMenuItem.Size = new Size(292, 22);
            load3DSToolStripMenuItem.Text = "Load External Model";
            load3DSToolStripMenuItem.Click += Load3DSToolStripMenuItem_Click;
            // 
            // loadTMDToolStripMenuItem
            // 
            loadTMDToolStripMenuItem.Name = "loadTMDToolStripMenuItem";
            loadTMDToolStripMenuItem.Size = new Size(292, 22);
            loadTMDToolStripMenuItem.Text = "Load TMD";
            loadTMDToolStripMenuItem.Click += LoadTMDToolStripMenuItem_Click;
            // 
            // toolStripSeparator1
            // 
            toolStripSeparator1.Name = "toolStripSeparator1";
            toolStripSeparator1.Size = new Size(289, 6);
            // 
            // saveSkeletonToolStripMenuItem
            // 
            saveSkeletonToolStripMenuItem.Enabled = false;
            saveSkeletonToolStripMenuItem.Name = "saveSkeletonToolStripMenuItem";
            saveSkeletonToolStripMenuItem.ShortcutKeys = Keys.Control | Keys.S;
            saveSkeletonToolStripMenuItem.Size = new Size(292, 22);
            saveSkeletonToolStripMenuItem.Text = "Save Skeleton";
            saveSkeletonToolStripMenuItem.Click += SaveSkeletonToolStripMenuItem_Click;
            // 
            // saveSkeletonAsToolStripMenuItem
            // 
            saveSkeletonAsToolStripMenuItem.Enabled = false;
            saveSkeletonAsToolStripMenuItem.Name = "saveSkeletonAsToolStripMenuItem";
            saveSkeletonAsToolStripMenuItem.Size = new Size(292, 22);
            saveSkeletonAsToolStripMenuItem.Text = "Save Skeleton/RSD/Model As...";
            saveSkeletonAsToolStripMenuItem.Click += SaveSkeletonAsToolStripMenuItem_Click;
            // 
            // toolStripSeparator9
            // 
            toolStripSeparator9.Name = "toolStripSeparator9";
            toolStripSeparator9.Size = new Size(289, 6);
            // 
            // Import3DSFixingPositionToolStripMenuItem
            // 
            Import3DSFixingPositionToolStripMenuItem.CheckOnClick = true;
            Import3DSFixingPositionToolStripMenuItem.Name = "Import3DSFixingPositionToolStripMenuItem";
            Import3DSFixingPositionToolStripMenuItem.Size = new Size(292, 22);
            Import3DSFixingPositionToolStripMenuItem.Text = "Fix position of external models on import";
            Import3DSFixingPositionToolStripMenuItem.Click += Import3DSFixingPositionToolStripMenuItem_Click;
            // 
            // DontCheckDuplicatedPolysVertsToolStripMenuItem
            // 
            DontCheckDuplicatedPolysVertsToolStripMenuItem.CheckOnClick = true;
            DontCheckDuplicatedPolysVertsToolStripMenuItem.Name = "DontCheckDuplicatedPolysVertsToolStripMenuItem";
            DontCheckDuplicatedPolysVertsToolStripMenuItem.Size = new Size(292, 22);
            DontCheckDuplicatedPolysVertsToolStripMenuItem.Text = "Don't Check Duplicated Polys/Verts";
            DontCheckDuplicatedPolysVertsToolStripMenuItem.Click += DontCheckRepeatedPolysVertsToolStripMenuItem_Click;
            // 
            // toolStripSeparator2
            // 
            toolStripSeparator2.Name = "toolStripSeparator2";
            toolStripSeparator2.Size = new Size(289, 6);
            // 
            // extiToolStripMenuItem
            // 
            extiToolStripMenuItem.Name = "extiToolStripMenuItem";
            extiToolStripMenuItem.Size = new Size(292, 22);
            extiToolStripMenuItem.Text = "Exit";
            extiToolStripMenuItem.Click += ExitToolStripMenuItem_Click;
            // 
            // editToolStripMenuItem2
            // 
            editToolStripMenuItem2.DropDownItems.AddRange(new ToolStripItem[] { undoToolStripMenuItem, redoToolStripMenuItem, toolStripSeparator6, ShowAxesToolStripMenuItem, toolStripSeparator10, resetCameraToolStripMenuItem, toolStripSeparator7, uIOpacityToolStripMenuItem });
            editToolStripMenuItem2.Name = "editToolStripMenuItem2";
            editToolStripMenuItem2.Size = new Size(39, 20);
            editToolStripMenuItem2.Text = "&Edit";
            // 
            // undoToolStripMenuItem
            // 
            undoToolStripMenuItem.Enabled = false;
            undoToolStripMenuItem.Name = "undoToolStripMenuItem";
            undoToolStripMenuItem.ShortcutKeys = Keys.Control | Keys.U;
            undoToolStripMenuItem.Size = new Size(221, 22);
            undoToolStripMenuItem.Text = "Undo";
            undoToolStripMenuItem.Click += UndoToolStripMenuItem_Click;
            // 
            // redoToolStripMenuItem
            // 
            redoToolStripMenuItem.Enabled = false;
            redoToolStripMenuItem.Name = "redoToolStripMenuItem";
            redoToolStripMenuItem.ShortcutKeys = Keys.Control | Keys.R;
            redoToolStripMenuItem.Size = new Size(221, 22);
            redoToolStripMenuItem.Text = "Redo";
            redoToolStripMenuItem.Click += RedoToolStripMenuItem_Click;
            // 
            // toolStripSeparator6
            // 
            toolStripSeparator6.Name = "toolStripSeparator6";
            toolStripSeparator6.Size = new Size(218, 6);
            // 
            // ShowAxesToolStripMenuItem
            // 
            ShowAxesToolStripMenuItem.CheckOnClick = true;
            ShowAxesToolStripMenuItem.Name = "ShowAxesToolStripMenuItem";
            ShowAxesToolStripMenuItem.Size = new Size(221, 22);
            ShowAxesToolStripMenuItem.Text = "Show Axes";
            ShowAxesToolStripMenuItem.Click += DrawAxesToolStripMenuItem_Click;
            // 
            // toolStripSeparator10
            // 
            toolStripSeparator10.Name = "toolStripSeparator10";
            toolStripSeparator10.Size = new Size(218, 6);
            // 
            // resetCameraToolStripMenuItem
            // 
            resetCameraToolStripMenuItem.Name = "resetCameraToolStripMenuItem";
            resetCameraToolStripMenuItem.ShortcutKeyDisplayString = "CTRL+Home";
            resetCameraToolStripMenuItem.Size = new Size(221, 22);
            resetCameraToolStripMenuItem.Text = "Reset Camera";
            resetCameraToolStripMenuItem.Click += ResetCameraToolStripMenuItem_Click;
            // 
            // toolStripSeparator7
            // 
            toolStripSeparator7.Name = "toolStripSeparator7";
            toolStripSeparator7.Size = new Size(218, 6);
            // 
            // uIOpacityToolStripMenuItem
            // 
            uIOpacityToolStripMenuItem.DropDownItems.AddRange(new ToolStripItem[] { tsUIOpacity100, tsUIOpacity90, tsUIOpacity75, tsUIOpacity50, tsUIOpacity25 });
            uIOpacityToolStripMenuItem.Name = "uIOpacityToolStripMenuItem";
            uIOpacityToolStripMenuItem.Size = new Size(221, 22);
            uIOpacityToolStripMenuItem.Text = "UI Opacity";
            // 
            // tsUIOpacity100
            // 
            tsUIOpacity100.Checked = true;
            tsUIOpacity100.CheckState = CheckState.Checked;
            tsUIOpacity100.Name = "tsUIOpacity100";
            tsUIOpacity100.Size = new Size(102, 22);
            tsUIOpacity100.Text = "100%";
            tsUIOpacity100.Click += TsUIOpacity100_Click;
            // 
            // tsUIOpacity90
            // 
            tsUIOpacity90.Name = "tsUIOpacity90";
            tsUIOpacity90.Size = new Size(102, 22);
            tsUIOpacity90.Text = "90%";
            tsUIOpacity90.Click += TsUIOpacity90_Click;
            // 
            // tsUIOpacity75
            // 
            tsUIOpacity75.Name = "tsUIOpacity75";
            tsUIOpacity75.Size = new Size(102, 22);
            tsUIOpacity75.Text = "75%";
            tsUIOpacity75.Click += TsUIOpacity75_Click;
            // 
            // tsUIOpacity50
            // 
            tsUIOpacity50.Name = "tsUIOpacity50";
            tsUIOpacity50.Size = new Size(102, 22);
            tsUIOpacity50.Text = "50%";
            tsUIOpacity50.Click += TsUIOpacity50_Click;
            // 
            // tsUIOpacity25
            // 
            tsUIOpacity25.Name = "tsUIOpacity25";
            tsUIOpacity25.Size = new Size(102, 22);
            tsUIOpacity25.Text = "25%";
            tsUIOpacity25.Click += TsUIOpacity25_Click;
            // 
            // skeletonToolStripMenuItem
            // 
            skeletonToolStripMenuItem.DropDownItems.AddRange(new ToolStripItem[] { addJointToolStripMenuItem, editJointToolStripMenuItem, toolStripSeparator8, showNormalsToolStripMenuItem, normalsColorToolStripMenuItem, normalsScaleToolStripMenuItem, toolStripSeparator11, statisticsToolStripMenuItem });
            skeletonToolStripMenuItem.Name = "skeletonToolStripMenuItem";
            skeletonToolStripMenuItem.Size = new Size(64, 20);
            skeletonToolStripMenuItem.Text = "&Skeleton";
            // 
            // addJointToolStripMenuItem
            // 
            addJointToolStripMenuItem.Enabled = false;
            addJointToolStripMenuItem.Name = "addJointToolStripMenuItem";
            addJointToolStripMenuItem.Size = new Size(151, 22);
            addJointToolStripMenuItem.Text = "Add Joint";
            addJointToolStripMenuItem.Click += AddJointToolStripMenuItem_Click;
            // 
            // editJointToolStripMenuItem
            // 
            editJointToolStripMenuItem.Enabled = false;
            editJointToolStripMenuItem.Name = "editJointToolStripMenuItem";
            editJointToolStripMenuItem.Size = new Size(151, 22);
            editJointToolStripMenuItem.Text = "Edit Joint";
            editJointToolStripMenuItem.Click += EditJointToolStripMenuItem_Click;
            // 
            // toolStripSeparator8
            // 
            toolStripSeparator8.Name = "toolStripSeparator8";
            toolStripSeparator8.Size = new Size(148, 6);
            // 
            // showNormalsToolStripMenuItem
            // 
            showNormalsToolStripMenuItem.DropDownItems.AddRange(new ToolStripItem[] { showVertexNormalsToolStripMenuItem, showFaceNormalsToolStripMenuItem });
            showNormalsToolStripMenuItem.Enabled = false;
            showNormalsToolStripMenuItem.Name = "showNormalsToolStripMenuItem";
            showNormalsToolStripMenuItem.Size = new Size(151, 22);
            showNormalsToolStripMenuItem.Text = "Show Normals";
            // 
            // showVertexNormalsToolStripMenuItem
            // 
            showVertexNormalsToolStripMenuItem.CheckOnClick = true;
            showVertexNormalsToolStripMenuItem.Name = "showVertexNormalsToolStripMenuItem";
            showVertexNormalsToolStripMenuItem.ShortcutKeyDisplayString = "V";
            showVertexNormalsToolStripMenuItem.Size = new Size(120, 22);
            showVertexNormalsToolStripMenuItem.Text = "Vertex";
            showVertexNormalsToolStripMenuItem.Click += ShowVertexNormalsToolStripMenuItem_Click;
            // 
            // showFaceNormalsToolStripMenuItem
            // 
            showFaceNormalsToolStripMenuItem.CheckOnClick = true;
            showFaceNormalsToolStripMenuItem.Name = "showFaceNormalsToolStripMenuItem";
            showFaceNormalsToolStripMenuItem.ShortcutKeyDisplayString = "F";
            showFaceNormalsToolStripMenuItem.Size = new Size(120, 22);
            showFaceNormalsToolStripMenuItem.Text = "Face";
            showFaceNormalsToolStripMenuItem.Click += ShowFaceNormalsToolStripMenuItem_Click;
            // 
            // normalsColorToolStripMenuItem
            // 
            normalsColorToolStripMenuItem.DropDownItems.AddRange(new ToolStripItem[] { redToolStripMenuItem, greenToolStripMenuItem, blueToolStripMenuItem });
            normalsColorToolStripMenuItem.Name = "normalsColorToolStripMenuItem";
            normalsColorToolStripMenuItem.Size = new Size(151, 22);
            normalsColorToolStripMenuItem.Text = "Normals Color";
            // 
            // redToolStripMenuItem
            // 
            redToolStripMenuItem.CheckOnClick = true;
            redToolStripMenuItem.Name = "redToolStripMenuItem";
            redToolStripMenuItem.Size = new Size(105, 22);
            redToolStripMenuItem.Text = "Red";
            redToolStripMenuItem.Click += RedToolStripMenuItem_Click;
            // 
            // greenToolStripMenuItem
            // 
            greenToolStripMenuItem.Checked = true;
            greenToolStripMenuItem.CheckOnClick = true;
            greenToolStripMenuItem.CheckState = CheckState.Checked;
            greenToolStripMenuItem.Name = "greenToolStripMenuItem";
            greenToolStripMenuItem.Size = new Size(105, 22);
            greenToolStripMenuItem.Text = "Green";
            greenToolStripMenuItem.Click += GreenToolStripMenuItem_Click;
            // 
            // blueToolStripMenuItem
            // 
            blueToolStripMenuItem.CheckOnClick = true;
            blueToolStripMenuItem.Name = "blueToolStripMenuItem";
            blueToolStripMenuItem.Size = new Size(105, 22);
            blueToolStripMenuItem.Text = "Blue";
            blueToolStripMenuItem.Click += BlueToolStripMenuItem_Click;
            // 
            // normalsScaleToolStripMenuItem
            // 
            normalsScaleToolStripMenuItem.DropDownItems.AddRange(new ToolStripItem[] { oneftoolStripMenuItem, fiveftoolStripMenuItem, thirtyftoolStripMenuItem, thousandftoolStripMenuItem });
            normalsScaleToolStripMenuItem.Name = "normalsScaleToolStripMenuItem";
            normalsScaleToolStripMenuItem.Size = new Size(151, 22);
            normalsScaleToolStripMenuItem.Text = "Normals Scale";
            // 
            // oneftoolStripMenuItem
            // 
            oneftoolStripMenuItem.Checked = true;
            oneftoolStripMenuItem.CheckOnClick = true;
            oneftoolStripMenuItem.CheckState = CheckState.Checked;
            oneftoolStripMenuItem.Name = "oneftoolStripMenuItem";
            oneftoolStripMenuItem.Size = new Size(216, 22);
            oneftoolStripMenuItem.Text = "1.0 (Field Models)";
            oneftoolStripMenuItem.Click += OneftoolStripMenuItem_Click;
            // 
            // fiveftoolStripMenuItem
            // 
            fiveftoolStripMenuItem.CheckOnClick = true;
            fiveftoolStripMenuItem.Name = "fiveftoolStripMenuItem";
            fiveftoolStripMenuItem.Size = new Size(216, 22);
            fiveftoolStripMenuItem.Text = "5.0";
            fiveftoolStripMenuItem.Click += FiveftoolStripMenuItem_Click;
            // 
            // thirtyftoolStripMenuItem
            // 
            thirtyftoolStripMenuItem.CheckOnClick = true;
            thirtyftoolStripMenuItem.Name = "thirtyftoolStripMenuItem";
            thirtyftoolStripMenuItem.Size = new Size(216, 22);
            thirtyftoolStripMenuItem.Text = "30.0 (Battle/Magic Models)";
            thirtyftoolStripMenuItem.Click += ThirtyftoolStripMenuItem_Click;
            // 
            // thousandftoolStripMenuItem
            // 
            thousandftoolStripMenuItem.CheckOnClick = true;
            thousandftoolStripMenuItem.Name = "thousandftoolStripMenuItem";
            thousandftoolStripMenuItem.Size = new Size(216, 22);
            thousandftoolStripMenuItem.Text = "1000.0 (Locations)";
            thousandftoolStripMenuItem.Click += ThousandftoolStripMenuItem_Click;
            // 
            // toolStripSeparator11
            // 
            toolStripSeparator11.Name = "toolStripSeparator11";
            toolStripSeparator11.Size = new Size(148, 6);
            // 
            // statisticsToolStripMenuItem
            // 
            statisticsToolStripMenuItem.Enabled = false;
            statisticsToolStripMenuItem.Name = "statisticsToolStripMenuItem";
            statisticsToolStripMenuItem.Size = new Size(151, 22);
            statisticsToolStripMenuItem.Text = "Statistics";
            statisticsToolStripMenuItem.Click += StatisticsToolStripMenuItem_Click;
            // 
            // textureToolStripMenuItem
            // 
            textureToolStripMenuItem.DropDownItems.AddRange(new ToolStripItem[] { TEXToPNGBatchConversionToolStripMenuItem });
            textureToolStripMenuItem.Name = "textureToolStripMenuItem";
            textureToolStripMenuItem.Size = new Size(57, 20);
            textureToolStripMenuItem.Text = "&Texture";
            // 
            // TEXToPNGBatchConversionToolStripMenuItem
            // 
            TEXToPNGBatchConversionToolStripMenuItem.Name = "TEXToPNGBatchConversionToolStripMenuItem";
            TEXToPNGBatchConversionToolStripMenuItem.Size = new Size(236, 22);
            TEXToPNGBatchConversionToolStripMenuItem.Text = ".TEX to .PNG Batch Conversion";
            TEXToPNGBatchConversionToolStripMenuItem.Click += TEXToPNGBatchConversionToolStripMenuItem_Click;
            // 
            // animationToolStripMenuItem
            // 
            animationToolStripMenuItem.DropDownItems.AddRange(new ToolStripItem[] { loadFieldAnimationToolStripMenuItem, loadBattleMagicLimitsAnimationStripMenuItem, toolStripSeparator3, saveAnimationToolStripMenuItem, saveAnimationAsToolStripMenuItem, outputFramesDataTXTToolStripMenuItem, inputFramesDataTXTToolStripMenuItem, inputFramesDataTXTToolSelectiveStripMenuItem, mergeFramesDataTXTToolStripMenuItem, toolStripSeparator4, toolStripMenuItem1, toolStripSeparator5, interpolateAllAnimationsToolStripMenuItem });
            animationToolStripMenuItem.Name = "animationToolStripMenuItem";
            animationToolStripMenuItem.Size = new Size(75, 20);
            animationToolStripMenuItem.Text = "&Animation";
            // 
            // loadFieldAnimationToolStripMenuItem
            // 
            loadFieldAnimationToolStripMenuItem.Enabled = false;
            loadFieldAnimationToolStripMenuItem.Name = "loadFieldAnimationToolStripMenuItem";
            loadFieldAnimationToolStripMenuItem.Size = new Size(382, 22);
            loadFieldAnimationToolStripMenuItem.Text = "Load Field Animation";
            loadFieldAnimationToolStripMenuItem.Click += LoadFieldAnimationToolStripMenuItem_Click;
            // 
            // loadBattleMagicLimitsAnimationStripMenuItem
            // 
            loadBattleMagicLimitsAnimationStripMenuItem.Enabled = false;
            loadBattleMagicLimitsAnimationStripMenuItem.Name = "loadBattleMagicLimitsAnimationStripMenuItem";
            loadBattleMagicLimitsAnimationStripMenuItem.Size = new Size(382, 22);
            loadBattleMagicLimitsAnimationStripMenuItem.Text = "Load Battle/Magic/Limits Animations";
            loadBattleMagicLimitsAnimationStripMenuItem.Click += LoadBattleMagicLimitAnimationsStripMenuItem_Click;
            // 
            // toolStripSeparator3
            // 
            toolStripSeparator3.Name = "toolStripSeparator3";
            toolStripSeparator3.Size = new Size(379, 6);
            // 
            // saveAnimationToolStripMenuItem
            // 
            saveAnimationToolStripMenuItem.Enabled = false;
            saveAnimationToolStripMenuItem.Name = "saveAnimationToolStripMenuItem";
            saveAnimationToolStripMenuItem.ShortcutKeys = Keys.Control | Keys.W;
            saveAnimationToolStripMenuItem.Size = new Size(382, 22);
            saveAnimationToolStripMenuItem.Text = "Save Animation";
            saveAnimationToolStripMenuItem.Click += SaveAnimationToolStripMenuItem_Click;
            // 
            // saveAnimationAsToolStripMenuItem
            // 
            saveAnimationAsToolStripMenuItem.Enabled = false;
            saveAnimationAsToolStripMenuItem.Name = "saveAnimationAsToolStripMenuItem";
            saveAnimationAsToolStripMenuItem.Size = new Size(382, 22);
            saveAnimationAsToolStripMenuItem.Text = "Save Animation As...";
            saveAnimationAsToolStripMenuItem.Click += SaveAnimationAsToolStripMenuItem_Click;
            // 
            // outputFramesDataTXTToolStripMenuItem
            // 
            outputFramesDataTXTToolStripMenuItem.Enabled = false;
            outputFramesDataTXTToolStripMenuItem.Name = "outputFramesDataTXTToolStripMenuItem";
            outputFramesDataTXTToolStripMenuItem.Size = new Size(382, 22);
            outputFramesDataTXTToolStripMenuItem.Text = "Output Frames Data to .TXT (Only Field Models)";
            outputFramesDataTXTToolStripMenuItem.Click += OutputFramesDataAsToolStripMenuItem_Click;
            // 
            // inputFramesDataTXTToolStripMenuItem
            // 
            inputFramesDataTXTToolStripMenuItem.Enabled = false;
            inputFramesDataTXTToolStripMenuItem.Name = "inputFramesDataTXTToolStripMenuItem";
            inputFramesDataTXTToolStripMenuItem.Size = new Size(382, 22);
            inputFramesDataTXTToolStripMenuItem.Text = "Input Frames Data from .TXT (Only Field Models)";
            inputFramesDataTXTToolStripMenuItem.Click += InputFramesDataTXTToolStripMenuItem_Click;
            // 
            // inputFramesDataTXTToolSelectiveStripMenuItem
            // 
            inputFramesDataTXTToolSelectiveStripMenuItem.Enabled = false;
            inputFramesDataTXTToolSelectiveStripMenuItem.Name = "inputFramesDataTXTToolSelectiveStripMenuItem";
            inputFramesDataTXTToolSelectiveStripMenuItem.Size = new Size(382, 22);
            inputFramesDataTXTToolSelectiveStripMenuItem.Text = "Input Frames Data from .TXT (Only Field Models, Selective)";
            inputFramesDataTXTToolSelectiveStripMenuItem.Click += InputFramesDataFromTXTOnlyFieldModelsSelectiveToolStripMenuItem_Click;
            // 
            // mergeFramesDataTXTToolStripMenuItem
            // 
            mergeFramesDataTXTToolStripMenuItem.Enabled = false;
            mergeFramesDataTXTToolStripMenuItem.Name = "mergeFramesDataTXTToolStripMenuItem";
            mergeFramesDataTXTToolStripMenuItem.Size = new Size(382, 22);
            mergeFramesDataTXTToolStripMenuItem.Text = "Merge Frames Data from .TXT (Only Field Models)";
            mergeFramesDataTXTToolStripMenuItem.Click += MergeFramesDataTXTOnlyFieldModelsToolStripMenuItem_Click;
            // 
            // toolStripSeparator4
            // 
            toolStripSeparator4.Name = "toolStripSeparator4";
            toolStripSeparator4.Size = new Size(379, 6);
            // 
            // toolStripMenuItem1
            // 
            toolStripMenuItem1.DropDownItems.AddRange(new ToolStripItem[] { toolStripFPS15, toolStripFPS30, toolStripFPS60 });
            toolStripMenuItem1.Name = "toolStripMenuItem1";
            toolStripMenuItem1.Size = new Size(382, 22);
            toolStripMenuItem1.Text = "FPS for Play Animation";
            // 
            // toolStripFPS15
            // 
            toolStripFPS15.Checked = true;
            toolStripFPS15.CheckState = CheckState.Checked;
            toolStripFPS15.Name = "toolStripFPS15";
            toolStripFPS15.Size = new Size(86, 22);
            toolStripFPS15.Text = "15";
            toolStripFPS15.Click += ToolStripFPS15_Click;
            // 
            // toolStripFPS30
            // 
            toolStripFPS30.Name = "toolStripFPS30";
            toolStripFPS30.Size = new Size(86, 22);
            toolStripFPS30.Text = "30";
            toolStripFPS30.Click += ToolStripFPS30_Click;
            // 
            // toolStripFPS60
            // 
            toolStripFPS60.Name = "toolStripFPS60";
            toolStripFPS60.Size = new Size(86, 22);
            toolStripFPS60.Text = "60";
            toolStripFPS60.Click += ToolStripFPS60_Click;
            // 
            // toolStripSeparator5
            // 
            toolStripSeparator5.Name = "toolStripSeparator5";
            toolStripSeparator5.Size = new Size(379, 6);
            // 
            // interpolateAllAnimationsToolStripMenuItem
            // 
            interpolateAllAnimationsToolStripMenuItem.Name = "interpolateAllAnimationsToolStripMenuItem";
            interpolateAllAnimationsToolStripMenuItem.Size = new Size(382, 22);
            interpolateAllAnimationsToolStripMenuItem.Text = "Interpolate All Animations";
            interpolateAllAnimationsToolStripMenuItem.Click += InterpolateAllAnimationsToolStripMenuItem_Click;
            // 
            // databaseToolStripMenuItem
            // 
            databaseToolStripMenuItem.DropDownItems.AddRange(new ToolStripItem[] { showCharlgpToolStripMenuItem, showWorldLgpToolStripMenuItem, showBattlelgpToolStripMenuItem, showMagiclgpToolStripMenuItem });
            databaseToolStripMenuItem.Name = "databaseToolStripMenuItem";
            databaseToolStripMenuItem.Size = new Size(67, 20);
            databaseToolStripMenuItem.Text = "&Database";
            // 
            // showCharlgpToolStripMenuItem
            // 
            showCharlgpToolStripMenuItem.Name = "showCharlgpToolStripMenuItem";
            showCharlgpToolStripMenuItem.ShortcutKeys = Keys.Control | Keys.D;
            showCharlgpToolStripMenuItem.Size = new Size(204, 22);
            showCharlgpToolStripMenuItem.Text = "Show CHAR.LGP";
            showCharlgpToolStripMenuItem.Click += ShowCharlgpToolStripMenuItem_Click;
            // 
            // showWorldLgpToolStripMenuItem
            // 
            showWorldLgpToolStripMenuItem.Name = "showWorldLgpToolStripMenuItem";
            showWorldLgpToolStripMenuItem.Size = new Size(204, 22);
            showWorldLgpToolStripMenuItem.Text = "Show WORLD_xx.LGP";
            showWorldLgpToolStripMenuItem.Click += showWorldLgpToolStripMenuItem_Click;
            // 
            // showBattlelgpToolStripMenuItem
            // 
            showBattlelgpToolStripMenuItem.Name = "showBattlelgpToolStripMenuItem";
            showBattlelgpToolStripMenuItem.Size = new Size(204, 22);
            showBattlelgpToolStripMenuItem.Text = "Show BATTLE.LGP";
            showBattlelgpToolStripMenuItem.Click += ShowBattlelgpToolStripMenuItem_Click;
            // 
            // showMagiclgpToolStripMenuItem
            // 
            showMagiclgpToolStripMenuItem.Name = "showMagiclgpToolStripMenuItem";
            showMagiclgpToolStripMenuItem.Size = new Size(204, 22);
            showMagiclgpToolStripMenuItem.Text = "Show MAGIC.LGP";
            showMagiclgpToolStripMenuItem.Click += ShowMagiclgpToolStripMenuItem_Click;
            // 
            // panelModel
            // 
            panelModel.Anchor = AnchorStyles.Top | AnchorStyles.Bottom | AnchorStyles.Left | AnchorStyles.Right;
            panelModel.API = OpenTK.Windowing.Common.ContextAPI.OpenGL;
            panelModel.APIVersion = new Version(3, 3, 0, 0);
            panelModel.BackColor = Color.Black;
            panelModel.Flags = OpenTK.Windowing.Common.ContextFlags.Default;
            panelModel.IsEventDriven = true;
            panelModel.Location = new Point(186, 29);
            panelModel.Margin = new Padding(4, 3, 4, 3);
            panelModel.Name = "panelModel";
            panelModel.Profile = OpenTK.Windowing.Common.ContextProfile.Compatability;
            panelModel.SharedContext = null;
            panelModel.Size = new Size(484, 638);
            panelModel.TabIndex = 9;
            panelModel.TabStop = false;
            panelModel.Paint += PanelModel_Paint;
            panelModel.DoubleClick += PanelModel_DoubleClick;
            panelModel.MouseDown += PanelModel_MouseDown;
            panelModel.MouseMove += PanelModel_MouseMove;
            panelModel.MouseUp += PanelModel_MouseUp;
            // 
            // FrmSkeletonEditor
            // 
            AutoScaleDimensions = new SizeF(7F, 15F);
            AutoScaleMode = AutoScaleMode.Font;
            AutoValidate = AutoValidate.EnablePreventFocusChange;
            ClientSize = new Size(856, 749);
            Controls.Add(panelModel);
            Controls.Add(panel3);
            Controls.Add(panel2);
            Controls.Add(panel1);
            Controls.Add(menuStrip1);
            DoubleBuffered = true;
            Icon = (Icon)resources.GetObject("$this.Icon");
            KeyPreview = true;
            MainMenuStrip = menuStrip1;
            Margin = new Padding(4, 3, 4, 3);
            MinimumSize = new Size(872, 771);
            Name = "FrmSkeletonEditor";
            StartPosition = FormStartPosition.CenterScreen;
            Text = "KimeraCS";
            Activated += FrmSkeletonEditor_Activated;
            FormClosed += FrmSkeletonEditor_FormClosed;
            Load += FrmSkeletonEditor_Load;
            KeyDown += FrmSkeletonEditor_KeyDown;
            KeyUp += FrmSkeletonEditor_KeyUp;
            Move += FrmSkeletonEditor_Move;
            Resize += FrmSkeletonEditor_Resize;
            panel1.ResumeLayout(false);
            panel1.PerformLayout();
            gbTexturesFrame.ResumeLayout(false);
            gbTexturesFrame.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)nUDMoveTextureUpDown).EndInit();
            ((System.ComponentModel.ISupportInitialize)pbTextureViewer).EndInit();
            gbAnimationOptionsFrame.ResumeLayout(false);
            gbAnimationOptionsFrame.PerformLayout();
            gbFrameDataPartOptions.ResumeLayout(false);
            gbFrameDataPartOptions.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)nUDZAnimationFramePart).EndInit();
            ((System.ComponentModel.ISupportInitialize)nUDYAnimationFramePart).EndInit();
            ((System.ComponentModel.ISupportInitialize)nUDXAnimationFramePart).EndInit();
            ((System.ComponentModel.ISupportInitialize)nUDFrameDataPart).EndInit();
            gbSelectedBoneFrame.ResumeLayout(false);
            gbSelectedBoneFrame.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)nUDBoneOptionsLength).EndInit();
            ((System.ComponentModel.ISupportInitialize)nUDResizeBoneZ).EndInit();
            ((System.ComponentModel.ISupportInitialize)nUDResizeBoneY).EndInit();
            ((System.ComponentModel.ISupportInitialize)nUDResizeBoneX).EndInit();
            panel2.ResumeLayout(false);
            panel2.PerformLayout();
            groupBox1.ResumeLayout(false);
            groupBox1.PerformLayout();
            gbSelectedPieceFrame.ResumeLayout(false);
            gbRotateFrame.ResumeLayout(false);
            gbRotateFrame.PerformLayout();
            gbRepositionFrame.ResumeLayout(false);
            gbRepositionFrame.PerformLayout();
            gbResizeFrame.ResumeLayout(false);
            gbResizeFrame.PerformLayout();
            panel3.ResumeLayout(false);
            panel3.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)tbCurrentFrameScroll).EndInit();
            menuStrip1.ResumeLayout(false);
            menuStrip1.PerformLayout();
            ResumeLayout(false);
            PerformLayout();

        }

        #endregion
        private System.Windows.Forms.Panel panel1;
        private System.Windows.Forms.Panel panel2;
        private System.Windows.Forms.Panel panel3;
        public OpenTK.GLControl.GLControl panelModel;
        private System.Windows.Forms.Label lblBoneSelector;
        private System.Windows.Forms.ComboBox cbBoneSelector;
        private System.Windows.Forms.Label lblAnimationFrame;
        private System.Windows.Forms.TextBox txtAnimationFrame;
        public System.Windows.Forms.CheckBox chkDListEnable;
        private System.Windows.Forms.GroupBox gbSelectedBoneFrame;
        private System.Windows.Forms.MenuStrip menuStrip1;
        private System.Windows.Forms.ToolStripMenuItem fileToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem loadFieldSkeletonToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem loadPModelToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem load3DSToolStripMenuItem;
        private System.Windows.Forms.ToolStripSeparator toolStripSeparator1;
        private System.Windows.Forms.ToolStripMenuItem saveSkeletonToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem saveSkeletonAsToolStripMenuItem;
        private System.Windows.Forms.ToolStripSeparator toolStripSeparator2;
        private System.Windows.Forms.ToolStripMenuItem extiToolStripMenuItem;
        private System.Windows.Forms.GroupBox gbSelectedPieceFrame;
        private System.Windows.Forms.ToolStripMenuItem animationToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem loadFieldAnimationToolStripMenuItem;
        private System.Windows.Forms.ToolStripSeparator toolStripSeparator3;
        private System.Windows.Forms.ToolStripMenuItem saveAnimationToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem saveAnimationAsToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem loadBattleMagicLimitsAnimationStripMenuItem;
        private System.Windows.Forms.NumericUpDown nUDResizeBoneZ;
        private System.Windows.Forms.Label lblZScaleBoneOptions;
        private System.Windows.Forms.NumericUpDown nUDResizeBoneY;
        private System.Windows.Forms.Label lblYScaleBoneOptions;
        private System.Windows.Forms.Label lblXScaleBoneOptions;
        private System.Windows.Forms.NumericUpDown nUDBoneOptionsLength;
        private System.Windows.Forms.Label lblBoneOptionsLength;
        private System.Windows.Forms.GroupBox gbAnimationOptionsFrame;
        private System.Windows.Forms.NumericUpDown nUDFrameDataPart;
        private System.Windows.Forms.Label lblFrameOptionsPart;
        private System.Windows.Forms.GroupBox gbFrameDataPartOptions;
        private System.Windows.Forms.Label lblZAnimationFramePart;
        private System.Windows.Forms.Label lblYAnimationFramePart;
        private System.Windows.Forms.Label lblXAnimationFramePart;
        private System.Windows.Forms.NumericUpDown nUDZAnimationFramePart;
        private System.Windows.Forms.NumericUpDown nUDYAnimationFramePart;
        private System.Windows.Forms.NumericUpDown nUDXAnimationFramePart;
        private System.Windows.Forms.Button btnRemovePiece;
        private System.Windows.Forms.Button btnAddPiece;
        private System.Windows.Forms.CheckBox chkColorKeyFlag;
        private System.Windows.Forms.Button btnChangeTexture;
        private System.Windows.Forms.Button btnFlipHorizontal;
        private System.Windows.Forms.Button btnFlipVertical;
        private System.Windows.Forms.Button btnRotate;
        private System.Windows.Forms.Button btnRemoveTexture;
        private System.Windows.Forms.Button btnAddTexture;
        private System.Windows.Forms.CheckBox chkPropagateChangesForward;
        private System.Windows.Forms.Button btnDuplicateFrame;
        private System.Windows.Forms.Button btnRemoveFrame;
        private System.Windows.Forms.Button btnInterpolateFrame;
        private System.Windows.Forms.GroupBox gbResizeFrame;
        private System.Windows.Forms.Label lblResizePieceX;
        private System.Windows.Forms.TextBox txtResizePieceX;
        private System.Windows.Forms.HScrollBar hsbResizePieceX;
        private System.Windows.Forms.TextBox txtResizePieceZ;
        private System.Windows.Forms.HScrollBar hsbResizePieceZ;
        private System.Windows.Forms.Label lblResizePieceZ;
        private System.Windows.Forms.TextBox txtResizePieceY;
        private System.Windows.Forms.HScrollBar hsbResizePieceY;
        private System.Windows.Forms.Label lblResizePieceY;
        private System.Windows.Forms.GroupBox gbRepositionFrame;
        private System.Windows.Forms.TextBox txtRepositionZ;
        private System.Windows.Forms.HScrollBar hsbRepositionZ;
        private System.Windows.Forms.Label lblRepositionZ;
        private System.Windows.Forms.TextBox txtRepositionY;
        private System.Windows.Forms.HScrollBar hsbRepositionY;
        private System.Windows.Forms.Label lblRepositionY;
        private System.Windows.Forms.Label lblRepositionX;
        private System.Windows.Forms.TextBox txtRepositionX;
        private System.Windows.Forms.HScrollBar hsbRepositionX;
        private System.Windows.Forms.GroupBox gbRotateFrame;
        private System.Windows.Forms.TextBox txtRotateGamma;
        private System.Windows.Forms.HScrollBar hsbRotateGamma;
        private System.Windows.Forms.Label lblRotateGamma;
        private System.Windows.Forms.TextBox txtRotateBeta;
        private System.Windows.Forms.HScrollBar hsbRotateBeta;
        private System.Windows.Forms.Label lblRotateBeta;
        private System.Windows.Forms.Label lblRotateAlpha;
        private System.Windows.Forms.TextBox txtRotateAlpha;
        private System.Windows.Forms.HScrollBar hsbRotateAlpha;
        public System.Windows.Forms.CheckBox chkShowLastFrameGhost;
        public System.Windows.Forms.CheckBox chkShowGround;
        private System.Windows.Forms.Button btnComputeWeaponPosition;
        private System.Windows.Forms.Button btnComputeGroundHeight;
        private System.Windows.Forms.Label lblWeapon;
        private System.Windows.Forms.Label lblBattleAnimation;
        private System.Windows.Forms.GroupBox groupBox1;
        private System.Windows.Forms.Label lblLightPosX;
        private System.Windows.Forms.HScrollBar hsbLightPosX;
        private System.Windows.Forms.Label lblLightPosZ;
        private System.Windows.Forms.HScrollBar hsbLightPosZ;
        private System.Windows.Forms.Label lblLightPosY;
        private System.Windows.Forms.HScrollBar hsbLightPosY;
        private System.Windows.Forms.CheckBox chkInifintyFarLights;
        private System.Windows.Forms.CheckBox chkRightLight;
        private System.Windows.Forms.CheckBox chkLeftLight;
        private System.Windows.Forms.CheckBox chkRearLight;
        private System.Windows.Forms.NumericUpDown nUDMoveTextureUpDown;
        private System.Windows.Forms.TextBox txtCopyPasteFrame;
        private System.Windows.Forms.Button btnPasteFrame;
        private System.Windows.Forms.Button btnCopyFrame;
        private System.Windows.Forms.ToolStripSeparator toolStripSeparator4;
        private System.Windows.Forms.ToolStripMenuItem interpolateAllAnimationsToolStripMenuItem;
        private System.Windows.Forms.CheckBox chkShowBones;
        private System.Windows.Forms.Button btnInterpolateAnimation;
        private System.Windows.Forms.ToolStripMenuItem databaseToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem showCharlgpToolStripMenuItem;
        public System.Windows.Forms.CheckBox chkFrontLight;
        public System.Windows.Forms.PictureBox pbTextureViewer;
        private System.Windows.Forms.Button btnFrameEnd;
        private System.Windows.Forms.Button btnFrameBegin;
        private System.Windows.Forms.CheckBox btnPlayStopAnim;
        private System.Windows.Forms.Button btnFramePrev;
        private System.Windows.Forms.Button btnFrameNext;
        public System.Windows.Forms.TrackBar tbCurrentFrameScroll;
        private System.Windows.Forms.ToolStripMenuItem toolStripMenuItem1;
        private System.Windows.Forms.ToolStripMenuItem toolStripFPS15;
        private System.Windows.Forms.ToolStripMenuItem toolStripFPS30;
        private System.Windows.Forms.ToolStripMenuItem toolStripFPS60;
        private System.Windows.Forms.ToolStripSeparator toolStripSeparator5;
        private System.Windows.Forms.ToolStripMenuItem loadBattleMagicSkeletonToolStripMenuItem;
        private System.Windows.Forms.OpenFileDialog openFile;
        private System.Windows.Forms.SaveFileDialog saveFile;
        private System.Windows.Forms.NumericUpDown nUDResizeBoneX;
        public System.Windows.Forms.ComboBox cbTextureSelect;
        public System.Windows.Forms.ComboBox cbWeapon;
        public System.Windows.Forms.ComboBox cbBattleAnimation;
        private System.Windows.Forms.ToolStripMenuItem editToolStripMenuItem2;
        private System.Windows.Forms.ToolStripSeparator toolStripSeparator6;
        private System.Windows.Forms.ToolStripMenuItem resetCameraToolStripMenuItem;
        public System.Windows.Forms.ToolStripMenuItem undoToolStripMenuItem;
        public System.Windows.Forms.ToolStripMenuItem redoToolStripMenuItem;
        private System.Windows.Forms.ToolStripSeparator toolStripSeparator7;
        private System.Windows.Forms.ToolStripMenuItem uIOpacityToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem tsUIOpacity100;
        private System.Windows.Forms.ToolStripMenuItem tsUIOpacity90;
        private System.Windows.Forms.ToolStripMenuItem tsUIOpacity75;
        private System.Windows.Forms.ToolStripMenuItem tsUIOpacity50;
        private System.Windows.Forms.ToolStripMenuItem tsUIOpacity25;
        private System.Windows.Forms.ToolStripMenuItem outputFramesDataTXTToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem inputFramesDataTXTToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem showBattlelgpToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem showMagiclgpToolStripMenuItem;
        public System.Windows.Forms.GroupBox gbTexturesFrame;
        private System.Windows.Forms.ToolStripMenuItem textureToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem TEXToPNGBatchConversionToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem loadRSDResourceToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem inputFramesDataTXTToolSelectiveStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem skeletonToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem addJointToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem editJointToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem loadTMDToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem mergeFramesDataTXTToolStripMenuItem;
        private System.Windows.Forms.ToolStripSeparator toolStripSeparator8;
        private System.Windows.Forms.ToolStripMenuItem showNormalsToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem normalsColorToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem redToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem greenToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem blueToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem normalsScaleToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem oneftoolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem fiveftoolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem thousandftoolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem thirtyftoolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem showVertexNormalsToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem showFaceNormalsToolStripMenuItem;
        private System.Windows.Forms.ToolStripSeparator toolStripSeparator9;
        private System.Windows.Forms.ToolStripMenuItem Import3DSFixingPositionToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem statisticsToolStripMenuItem;
        private System.Windows.Forms.ToolStripSeparator toolStripSeparator11;
        private System.Windows.Forms.ToolStripMenuItem ShowAxesToolStripMenuItem;
        private System.Windows.Forms.ToolStripSeparator toolStripSeparator10;
        private System.Windows.Forms.ToolStripMenuItem DontCheckDuplicatedPolysVertsToolStripMenuItem;
        private ToolStripMenuItem showWorldLgpToolStripMenuItem;
    }
}

