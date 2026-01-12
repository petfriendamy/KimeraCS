using KimeraCS.Core;

namespace KimeraCS.Rendering
{
    using static FF7PModel;
    using static FF7FieldSkeleton;
    using static FF7FieldAnimation;
    using static FF7BattleSkeleton;
    using static FF7BattleAnimationsPack;
    using static ModelDrawing;

    /// <summary>
    /// Camera state for 3D viewport rendering.
    /// </summary>
    public struct CameraState
    {
        public float Alpha;     // X-axis rotation
        public float Beta;      // Y-axis rotation
        public float Gamma;     // Z-axis rotation
        public float Distance;  // Camera distance from origin
        public float PanX;      // Horizontal pan offset
        public float PanY;      // Vertical pan offset
        public float PanZ;      // Depth pan offset

        public static CameraState Default => new CameraState
        {
            Alpha = 0,
            Beta = 0,
            Gamma = 0,
            Distance = -10,
            PanX = 0,
            PanY = 0,
            PanZ = 0
        };
    }

    /// <summary>
    /// Model transformation state (resize/reposition).
    /// </summary>
    public struct ModelTransform
    {
        public float ResizeX;
        public float ResizeY;
        public float ResizeZ;
        public float RepositionX;
        public float RepositionY;
        public float RepositionZ;

        public static ModelTransform Default => new ModelTransform
        {
            ResizeX = 1,
            ResizeY = 1,
            ResizeZ = 1,
            RepositionX = 0,
            RepositionY = 0,
            RepositionZ = 0
        };
    }

    /// <summary>
    /// Rendering options for skeleton/model display.
    /// </summary>
    public struct RenderingOptions
    {
        public bool ShowBones;
        public bool ShowGround;
        public NormalsDisplayMode NormalsDisplayMode;
        public bool ShowLastFrameGhost;
        public bool EnableDisplayLists;
        public float NormalsScale;
        public int NormalsColor;

        public static RenderingOptions Default => new RenderingOptions
        {
            ShowBones = false,
            ShowGround = false,
            NormalsDisplayMode = NormalsDisplayMode.None,
            ShowLastFrameGhost = false,
            EnableDisplayLists = true,
            NormalsScale = DEFAULT_NORMAL_SCALE,
            NormalsColor = DEFAULT_NORMAL_COLOR
        };
    }

    /// <summary>
    /// Lighting configuration for the renderer.
    /// </summary>
    public struct LightingConfig
    {
        public bool FrontLightEnabled;
        public bool RearLightEnabled;
        public bool RightLightEnabled;
        public bool LeftLightEnabled;
        public float PosXScroll;
        public float PosYScroll;
        public float PosZScroll;

        public static LightingConfig Default => new LightingConfig
        {
            FrontLightEnabled = true,
            RearLightEnabled = false,
            RightLightEnabled = false,
            LeftLightEnabled = false,
            PosXScroll = 0,
            PosYScroll = 0,
            PosZScroll = 0
        };

        public bool AnyLightEnabled => FrontLightEnabled || RearLightEnabled || RightLightEnabled || LeftLightEnabled;
    }

    /// <summary>
    /// Animation playback state.
    /// </summary>
    public struct AnimationState
    {
        public int CurrentFrame;
        public int AnimationIndex;
        public int WeaponAnimationIndex;

        public static AnimationState Default => new AnimationState
        {
            CurrentFrame = 0,
            AnimationIndex = 0,
            WeaponAnimationIndex = -1
        };
    }

    /// <summary>
    /// Selection state for skeleton editing.
    /// </summary>
    public struct SelectionState
    {
        public int SelectedBone;
        public int SelectedBonePiece;

        public static SelectionState None => new SelectionState
        {
            SelectedBone = -1,
            SelectedBonePiece = -1
        };
    }

    /// <summary>
    /// Model data for rendering. Contains the actual models and animations
    /// needed for different model types.
    /// </summary>
    public class SkeletonModelData
    {
        // P-Model data (for K_P_FIELD_MODEL, K_P_BATTLE_MODEL, K_P_MAGIC_MODEL, K_3DS_MODEL)
        public PModel PModel { get; set; }
        public uint[] TextureIds { get; set; }

        // Field skeleton data (for K_HRC_SKELETON)
        public FieldSkeleton FieldSkeleton { get; set; }
        public FieldAnimation FieldAnimation { get; set; }

        // Battle skeleton data (for K_AA_SKELETON, K_MAGIC_SKELETON)
        public BattleSkeleton BattleSkeleton { get; set; }
        public BattleAnimationsPack BattleAnimations { get; set; }
    }

    /// <summary>
    /// Complete rendering context combining all state needed for rendering.
    /// This allows Core rendering functions to be independent of form state.
    /// </summary>
    public class RenderingContext
    {
        public ModelType ModelType { get; set; }
        public bool IsLoaded { get; set; }
        public CameraState Camera { get; set; }
        public ModelTransform Transform { get; set; }
        public RenderingOptions Options { get; set; }
        public LightingConfig Lighting { get; set; }
        public AnimationState Animation { get; set; }
        public SelectionState Selection { get; set; }
        public SkeletonModelData ModelData { get; set; }

        public RenderingContext()
        {
            Camera = CameraState.Default;
            Transform = ModelTransform.Default;
            Options = RenderingOptions.Default;
            Lighting = LightingConfig.Default;
            Animation = AnimationState.Default;
            Selection = SelectionState.None;
            ModelData = null; // Must be set by caller if using context-based rendering
        }

        /// <summary>
        /// Creates a context from the current FrmSkeletonEditor state.
        /// Call this from the form to get a snapshot of current state.
        /// </summary>
        public static RenderingContext FromSkeletonEditor(
            ModelType modelType,
            bool isLoaded,
            float alpha, float beta, float gamma,
            float distance,
            float panX, float panY, float panZ,
            int currentFrame,
            int animIndex,
            int weaponAnimIndex,
            int selectedBone,
            int selectedBonePiece,
            bool showBones,
            bool showGround,
            NormalsDisplayMode normalsViewMode,
            bool showLastFrameGhost,
            bool enableDisplayLists,
            float normalsScale,
            int normalsColor,
            bool frontLight,
            bool rearLight,
            bool rightLight,
            bool leftLight,
            float lightPosX,
            float lightPosY,
            float lightPosZ)
        {
            return new RenderingContext
            {
                ModelType = modelType,
                IsLoaded = isLoaded,
                Camera = new CameraState
                {
                    Alpha = alpha,
                    Beta = beta,
                    Gamma = gamma,
                    Distance = distance,
                    PanX = panX,
                    PanY = panY,
                    PanZ = panZ
                },
                Animation = new AnimationState
                {
                    CurrentFrame = currentFrame,
                    AnimationIndex = animIndex,
                    WeaponAnimationIndex = weaponAnimIndex
                },
                Selection = new SelectionState
                {
                    SelectedBone = selectedBone,
                    SelectedBonePiece = selectedBonePiece
                },
                Options = new RenderingOptions
                {
                    ShowBones = showBones,
                    ShowGround = showGround,
                    NormalsDisplayMode = normalsViewMode,
                    ShowLastFrameGhost = showLastFrameGhost,
                    EnableDisplayLists = enableDisplayLists,
                    NormalsScale = normalsScale,
                    NormalsColor = normalsColor
                },
                Lighting = new LightingConfig
                {
                    FrontLightEnabled = frontLight,
                    RearLightEnabled = rearLight,
                    RightLightEnabled = rightLight,
                    LeftLightEnabled = leftLight,
                    PosXScroll = lightPosX,
                    PosYScroll = lightPosY,
                    PosZScroll = lightPosZ
                }
            };
        }

        /// <summary>
        /// Creates a context with model data for fully decoupled rendering.
        /// Use this from external applications that have their own model data.
        /// </summary>
        public static RenderingContext CreateWithModelData(
            ModelType modelType,
            SkeletonModelData modelData,
            CameraState camera,
            AnimationState animation,
            SelectionState selection,
            RenderingOptions options,
            LightingConfig lighting)
        {
            return new RenderingContext
            {
                ModelType = modelType,
                IsLoaded = true,
                Camera = camera,
                Animation = animation,
                Selection = selection,
                Options = options,
                Lighting = lighting,
                ModelData = modelData
            };
        }
    }

    /// <summary>
    /// Context for P Editor rendering, separate from skeleton editor.
    /// </summary>
    public class PEditorContext
    {
        public PModel EditedModel { get; set; }
        public uint[] TextureIds { get; set; }
        public CameraState Camera { get; set; }
        public ModelTransform Transform { get; set; }
        public DrawMode DrawMode { get; set; }
        public bool EnableLighting { get; set; }
        public int LightX { get; set; }
        public int LightY { get; set; }
        public int LightZ { get; set; }

        public PEditorContext()
        {
            Camera = CameraState.Default;
            Transform = ModelTransform.Default;
            DrawMode = DrawMode.K_VCOLORS;
            EnableLighting = false;
            LightX = 0;
            LightY = 1;
            LightZ = 0;
        }

        /// <summary>
        /// Creates a context from FrmPEditor state.
        /// </summary>
        public static PEditorContext FromPEditor(
            PModel editedModel,
            uint[] texIds,
            float alpha, float beta, float gamma,
            float distance,
            float panX, float panY, float panZ,
            float resizeX, float resizeY, float resizeZ,
            float reposX, float reposY, float reposZ,
            DrawMode drawMode,
            bool enableLighting,
            int lightX, int lightY, int lightZ)
        {
            return new PEditorContext
            {
                EditedModel = editedModel,
                TextureIds = texIds,
                Camera = new CameraState
                {
                    Alpha = alpha,
                    Beta = beta,
                    Gamma = gamma,
                    Distance = distance,
                    PanX = panX,
                    PanY = panY,
                    PanZ = panZ
                },
                Transform = new ModelTransform
                {
                    ResizeX = resizeX,
                    ResizeY = resizeY,
                    ResizeZ = resizeZ,
                    RepositionX = reposX,
                    RepositionY = reposY,
                    RepositionZ = reposZ
                },
                DrawMode = drawMode,
                EnableLighting = enableLighting,
                LightX = lightX,
                LightY = lightY,
                LightZ = lightZ
            };
        }
    }
}
