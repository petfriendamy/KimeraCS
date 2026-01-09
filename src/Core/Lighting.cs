using OpenTK.Graphics.OpenGL.Compatibility;
using OpenTK.Mathematics;
using KimeraCS.Rendering;

namespace KimeraCS
{
    using static FrmSkeletonEditor;

    using static FF7Skeleton;
    using static FF7FieldSkeleton;
    using static FF7PModel;

    using static FF7BattleSkeleton;

    using static Utils;

    class Lighting
    {
        public const int LIGHT_STEPS = 20;

        // Light indices (matches GLRenderer arrays)
        public const int LIGHT_RIGHT = 0;
        public const int LIGHT_LEFT = 1;
        public const int LIGHT_FRONT = 2;
        public const int LIGHT_REAR = 3;

        /// <summary>
        /// Modern lighting setup - sets GLRenderer light properties instead of legacy GL.Light calls.
        /// Call this for modern shader-based rendering.
        /// </summary>
        public static void SetLights()
        {
            Point3D p_min = new Point3D();
            Point3D p_max = new Point3D();

            // Check if any lights are enabled
            bool anyLightEnabled = bchkFrontLight || bchkRearLight || bchkRightLight || bchkLeftLight;
            GLRenderer.LightingEnabled = anyLightEnabled;

            if (!anyLightEnabled)
                return;

            // Compute scene bounds for light positioning
            switch (modelType)
            {
                case K_P_FIELD_MODEL:
                case K_P_BATTLE_MODEL:
                case K_P_MAGIC_MODEL:
                case K_3DS_MODEL:
                    ComputePModelBoundingBox(fPModel, ref p_min, ref p_max);
                    break;

                case K_HRC_SKELETON:
                    ComputeFieldBoundingBox(fSkeleton, fAnimation.frames[iCurrentFrameScroll], ref p_min, ref p_max);
                    break;

                case K_AA_SKELETON:
                case K_MAGIC_SKELETON:
                    ComputeBattleBoundingBox(bSkeleton, bAnimationsPack.SkeletonAnimations[ianimIndex].frames[iCurrentFrameScroll],
                                             ref p_min, ref p_max);
                    break;
            }

            float scene_diameter = (float)(-2 * ComputeSceneRadius(p_min, p_max));

            float light_x = scene_diameter / LIGHT_STEPS * fLightPosXScroll;
            float light_y = scene_diameter / LIGHT_STEPS * fLightPosYScroll;
            float light_z = scene_diameter / LIGHT_STEPS * fLightPosZScroll;

            // Set light positions and enable flags
            // Right light
            GLRenderer.LightEnabled[LIGHT_RIGHT] = bchkRightLight;
            if (bchkRightLight)
            {
                GLRenderer.LightPositions[LIGHT_RIGHT] = new Vector3(light_z, light_y, light_x);
                GLRenderer.LightColors[LIGHT_RIGHT] = new Vector3(0.5f, 0.5f, 0.5f);
            }

            // Left light
            GLRenderer.LightEnabled[LIGHT_LEFT] = bchkLeftLight;
            if (bchkLeftLight)
            {
                GLRenderer.LightPositions[LIGHT_LEFT] = new Vector3(-light_z, light_y, light_x);
                GLRenderer.LightColors[LIGHT_LEFT] = new Vector3(0.5f, 0.5f, 0.5f);
            }

            // Front light
            GLRenderer.LightEnabled[LIGHT_FRONT] = bchkFrontLight;
            if (bchkFrontLight)
            {
                GLRenderer.LightPositions[LIGHT_FRONT] = new Vector3(light_x, light_y, light_z);
                GLRenderer.LightColors[LIGHT_FRONT] = new Vector3(1f, 1f, 1f);
            }

            // Rear light
            GLRenderer.LightEnabled[LIGHT_REAR] = bchkRearLight;
            if (bchkRearLight)
            {
                GLRenderer.LightPositions[LIGHT_REAR] = new Vector3(light_x, light_y, -light_z);
                GLRenderer.LightColors[LIGHT_REAR] = new Vector3(0.75f, 0.75f, 0.75f);
            }
        }

        public static void SetLighting(OpenTK.Graphics.OpenGL.Compatibility.LightName lightNumber, float x, float y, float z, float red, float green, float blue, bool infinityFarQ)
        {
            float[] l_color = new float[4];
            float[] l_pos = new float[4];

            l_pos[0] = x;
            l_pos[1] = y;
            l_pos[2] = z;

            if (infinityFarQ) l_pos[3] = 0;
            else l_pos[3] = 1;

            l_color[0] = red;
            l_color[1] = green;
            l_color[2] = blue;
            l_color[3] = 1;

            GL.Enable(EnableCap.Lighting);
            // Cast LightName to EnableCap - they share the same underlying values
            GL.Disable((EnableCap)lightNumber);

            GL.Lightf(lightNumber, LightParameter.Position, l_pos);
            GL.Lightf(lightNumber, LightParameter.Diffuse, l_color);

            GL.Enable((EnableCap)lightNumber);
        }
    }
}
