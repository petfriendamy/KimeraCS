using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;

using OpenTK.Graphics.OpenGL;
using OpenTK.Mathematics;

namespace KimeraCS
{
    using Rendering;

    using static FF7FieldAnimation;
    using static FF7FieldRSDResource;

    using static FF7TEXTexture;
    using static FF7PModel;

    using static ModelDrawing;

    using static Utils;

    public class FF7FieldSkeleton
    {

        //
        // Field Skeleton Structure
        //
        public struct FieldSkeleton
        {
            public string fileName;
            public string name;
            public int nBones;

            public List<FieldBone> bones;
            public List<TEX> textures_pool;

            public FieldSkeleton(string strfileName, bool loadGeometryQ)
            {
                string strFileDirectoryName = Path.GetDirectoryName(strfileName);

                textures_pool = new List<TEX>();

                fileName = Path.GetFileName(strfileName).ToUpper();

                // Let's read HRC file into memory.
                string[] hrcString = File.ReadAllLines(strfileName);

                // Let's process the lines.
                // Skeleton name
                name = hrcString[1].Substring(10);

                // Skeleton number of bones
                nBones = Int32.Parse(hrcString[2].Substring(7));
                // Check model without skeleton

                bones = new List<FieldBone>();

                // Populate bones.
                // There is always "root" even there are no bones in some models.

                //  We will get the rest of bones rows checking if we are reading a correct line from
                //  .HRC file. There are some Field Models that can have '#' lines not useful among other things maybe.
                int i = 4;
                int numLine = 0;
                string rowOne = "", rowTwo = "", rowThree = "", rowFour = "";

                while (i < hrcString.Length)
                {
                    hrcString[i] = hrcString[i].Trim();

                    if (hrcString[i] != "" && hrcString[i][0] != '#' && hrcString[i][0] != ' ')
                    {
                        switch (numLine)
                        {
                            case 0:
                                rowOne = hrcString[i];
                                numLine++;
                                break;
                            case 1:
                                rowTwo = hrcString[i];
                                numLine++;
                                break;
                            case 2:
                                rowThree = hrcString[i];
                                numLine++;
                                break;
                            case 3:
                                rowFour = hrcString[i];
                                numLine = 0;
                                break;
                        }
                        
                        if (numLine == 0)
                        {
                            bones.Add(new FieldBone(rowOne,
                                                    rowTwo,
                                                    rowThree,
                                                    rowFour,
                                                    ref textures_pool,
                                                    loadGeometryQ,
                                                    strFileDirectoryName));
                        }
                    }

                    i++;
                }
            }
        }

        //
        // Field Skeleton Bone Structure
        //
        public struct FieldBone
        {
            public int nResources;
            public List<FieldRSDResource> fRSDResources;
            public double len;
            public string joint_i;
            public string joint_f;
            // added attributes
            public float resizeX;
            public float resizeY;
            public float resizeZ;

            public FieldBone(string inJointI, string inJointF, string inLen, string inRSDLine,
                             ref List<TEX> texturesPool, bool loadgeometryQ, string strFolderName)
            {
                string[] rsdRes = inRSDLine.Split();
                int i;

                len = Double.Parse(inLen, CultureInfo.InvariantCulture);

                joint_i = "";
                joint_f = "";

                resizeX = 0;
                resizeY = 0;
                resizeZ = 0;

                nResources = 0;

                fRSDResources = new List<FieldRSDResource>();

                if (loadgeometryQ)
                {
                    if (inRSDLine.Length > 2)
                        nResources = Int32.Parse(rsdRes[0].Substring(0, inRSDLine.IndexOf(" ")));

                    joint_i = inJointI;
                    joint_f = inJointF;

                    resizeX = 1;
                    resizeY = 1;
                    resizeZ = 1;

                    // Populate resources (RSD files)
                    if (nResources > 0)
                    {
                        if (nResources > rsdRes.Length - 1) nResources = rsdRes.Length - 1;

                        for (i = 0; i < nResources && rsdRes[i + 1] != null; i++)
                        {
                            fRSDResources.Add(new FieldRSDResource(rsdRes[i + 1], ref texturesPool, strFolderName));
                        }
                    }
                }
            }
        }


        //
        // Field Skeleton functions
        //
        public static void ComputeFieldBoneBoundingBox(FieldBone bone, ref Point3D p_min, ref Point3D p_max)
        {
            short ri;
            Point3D p_min_part = new Point3D();
            Point3D p_max_part = new Point3D();

            if (bone.nResources == 0)
            {
                p_max.x = 0;
                p_max.y = 0;
                p_max.z = 0;

                p_min.x = 0;
                p_min.y = 0;
                p_min.z = 0;
            }
            else
            {
                p_max.x = -(float)INFINITY_SINGLE;
                p_max.y = -(float)INFINITY_SINGLE;
                p_max.z = -(float)INFINITY_SINGLE;

                p_min.x = (float)INFINITY_SINGLE;
                p_min.y = (float)INFINITY_SINGLE;
                p_min.z = (float)INFINITY_SINGLE;

                for (ri = 0; ri < bone.nResources; ri++)
                {
                    ComputePModelBoundingBox(bone.fRSDResources[ri].Model, ref p_min_part, ref p_max_part);

                    if (p_max.x < p_max_part.x) p_max.x = p_max_part.x;
                    if (p_max.y < p_max_part.y) p_max.y = p_max_part.y;
                    if (p_max.z < p_max_part.z) p_max.z = p_max_part.z;

                    if (p_min.x > p_min_part.x) p_min.x = p_min_part.x;
                    if (p_min.y > p_min_part.y) p_min.y = p_min_part.y;
                    if (p_min.z > p_min_part.z) p_min.z = p_min_part.z;

                }
            }
        }

        public static void ComputeFieldBoundingBox(FieldSkeleton fSkeleton, FieldFrame fFrame,
                                                   ref Point3D p_min_field, ref Point3D p_max_field)
        {
            string[] joint_stack;
            int jsp;
            double[] rot_mat = new double[16];
            double[] MV_matrix = new double[16];
            int bi;

            Point3D p_max_bone = new Point3D();
            Point3D p_min_bone = new Point3D();

            Point3D p_max_bone_trans = new Point3D();
            Point3D p_min_bone_trans = new Point3D();

            joint_stack = new string[fSkeleton.bones.Count + 1];
            Matrix4[] matrixStack = new Matrix4[fSkeleton.bones.Count + 2];
            int matrixStackPtr = 0;
            jsp = 0;

            joint_stack[jsp] = fSkeleton.bones[0].joint_f;

            p_max_field.x = -(float)INFINITY_SINGLE;
            p_max_field.y = -(float)INFINITY_SINGLE;
            p_max_field.z = -(float)INFINITY_SINGLE;

            p_min_field.x = (float)INFINITY_SINGLE;
            p_min_field.y = (float)INFINITY_SINGLE;
            p_min_field.z = (float)INFINITY_SINGLE;

            // Build initial transform
            Matrix4 currentMatrix = Matrix4.Identity;
            currentMatrix *= Matrix4.CreateTranslation((float)fFrame.rootTranslationX, 0, 0);
            currentMatrix *= Matrix4.CreateTranslation(0, (float)(-fFrame.rootTranslationY), 0);
            currentMatrix *= Matrix4.CreateTranslation(0, 0, (float)fFrame.rootTranslationZ);

            BuildRotationMatrixWithQuaternions(fFrame.rootRotationAlpha, fFrame.rootRotationBeta, fFrame.rootRotationGamma, ref rot_mat);
            Matrix4 rotMatrix = DoubleArrayToMatrix4(rot_mat);
            currentMatrix *= rotMatrix;

            matrixStack[matrixStackPtr++] = currentMatrix;

            for (bi = 0; bi < fSkeleton.bones.Count; bi++)
            {
                while (!(fSkeleton.bones[bi].joint_f == joint_stack[jsp]) && jsp > 0)
                {
                    matrixStackPtr--;
                    currentMatrix = matrixStack[matrixStackPtr - 1];
                    jsp--;
                }
                matrixStack[matrixStackPtr++] = currentMatrix;

                BuildRotationMatrixWithQuaternions(fFrame.rotations[bi].alpha, fFrame.rotations[bi].beta,
                                                         fFrame.rotations[bi].gamma, ref rot_mat);
                rotMatrix = DoubleArrayToMatrix4(rot_mat);
                currentMatrix *= rotMatrix;

                ComputeFieldBoneBoundingBox(fSkeleton.bones[bi], ref p_min_bone, ref p_max_bone);

                MV_matrix = Matrix4ToDoubleArray(currentMatrix);

                ComputeTransformedBoxBoundingBox(MV_matrix, ref p_min_bone, ref p_max_bone, ref p_min_bone_trans, ref p_max_bone_trans);

                if (p_max_field.x < p_max_bone_trans.x) p_max_field.x = p_max_bone_trans.x;
                if (p_max_field.y < p_max_bone_trans.y) p_max_field.y = p_max_bone_trans.y;
                if (p_max_field.z < p_max_bone_trans.z) p_max_field.z = p_max_bone_trans.z;

                if (p_min_field.x > p_min_bone_trans.x) p_min_field.x = p_min_bone_trans.x;
                if (p_min_field.y > p_min_bone_trans.y) p_min_field.y = p_min_bone_trans.y;
                if (p_min_field.z > p_min_bone_trans.z) p_min_field.z = p_min_bone_trans.z;

                currentMatrix *= Matrix4.CreateTranslation(0, 0, (float)(-fSkeleton.bones[bi].len));

                jsp++;
                joint_stack[jsp] = fSkeleton.bones[bi].joint_i;
            }
        }

        public static float ComputeFieldBoneDiameter(FieldBone bone)
        {
            int iResourceIdx;
            float computeDiameter = 0;

            Point3D p_max = new Point3D();
            Point3D p_min = new Point3D();

            if (bone.nResources > 0)
            {
                p_max.x = -(float)INFINITY_SINGLE;
                p_max.y = -(float)INFINITY_SINGLE;
                p_max.z = -(float)INFINITY_SINGLE;

                p_min.x = (float)INFINITY_SINGLE;
                p_min.y = (float)INFINITY_SINGLE;
                p_min.z = (float)INFINITY_SINGLE;

                for (iResourceIdx = 0; iResourceIdx < bone.nResources; iResourceIdx++)
                {
                    if (p_max.x < bone.fRSDResources[iResourceIdx].Model.BoundingBox.max_x)
                        p_max.x = bone.fRSDResources[iResourceIdx].Model.BoundingBox.max_x;
                    if (p_max.y < bone.fRSDResources[iResourceIdx].Model.BoundingBox.max_y)
                        p_max.y = bone.fRSDResources[iResourceIdx].Model.BoundingBox.max_y;
                    if (p_max.z < bone.fRSDResources[iResourceIdx].Model.BoundingBox.max_z)
                        p_max.z = bone.fRSDResources[iResourceIdx].Model.BoundingBox.max_z;

                    if (p_min.x > bone.fRSDResources[iResourceIdx].Model.BoundingBox.min_x)
                        p_min.x = bone.fRSDResources[iResourceIdx].Model.BoundingBox.min_x;
                    if (p_min.y > bone.fRSDResources[iResourceIdx].Model.BoundingBox.min_y)
                        p_min.y = bone.fRSDResources[iResourceIdx].Model.BoundingBox.min_y;
                    if (p_min.z > bone.fRSDResources[iResourceIdx].Model.BoundingBox.min_z)
                        p_min.z = bone.fRSDResources[iResourceIdx].Model.BoundingBox.min_z;
                }

                computeDiameter = CalculateDistance(p_max, p_min);
            }

            return computeDiameter;
        }

        public static float ComputeFieldDiameter(FieldSkeleton Skeleton)
        {
            int iBoneIdx;
            float aux_diam;
            float fcomputeFieldDiameterResult = 0;

            //for (bi = 0; bi < Skeleton.nBones; bi++)
            for (iBoneIdx = 0; iBoneIdx < Skeleton.bones.Count; iBoneIdx++)
            {
                aux_diam = ComputeFieldBoneDiameter(Skeleton.bones[iBoneIdx]);

                fcomputeFieldDiameterResult += aux_diam;
            }

            return fcomputeFieldDiameterResult;
        }

        public static void SelectFieldBoneAndPiece(FieldSkeleton fSkeleton, FieldFrame fFrame, int b_index, int p_index)
        {
            int i, jsp;
            double[] rot_mat = new double[16];

            GL.MatrixMode(MatrixMode.Modelview);
            GL.PushMatrix();

            GL.Translate(fFrame.rootTranslationX, 0, 0);
            GL.Translate(0, -fFrame.rootTranslationY, 0);
            GL.Translate(0, 0, fFrame.rootTranslationZ);

            BuildRotationMatrixWithQuaternions(fFrame.rootRotationAlpha, fFrame.rootRotationBeta, fFrame.rootRotationGamma, ref rot_mat);
            GL.MultMatrix(rot_mat);

            if (b_index > -1)
            {
                jsp = MoveToFieldBone(fSkeleton, fFrame, b_index);
                DrawFieldBoneBoundingBox(fSkeleton.bones[b_index]);

                if (p_index > -1)
                    DrawFieldBonePieceBoundingBox(fSkeleton.bones[b_index], p_index);

                for (i = 0; i <= jsp; i++) GL.PopMatrix();
            }

            GL.PopMatrix();
        }

        /// <summary>
        /// Möller–Trumbore ray-triangle intersection algorithm.
        /// </summary>
        private static bool RayTriangleIntersect(Vector3 rayOrigin, Vector3 rayDir,
                                                  Vector3 v0, Vector3 v1, Vector3 v2,
                                                  out float distance)
        {
            distance = 0;
            const float EPSILON = 0.0000001f;

            Vector3 edge1 = v1 - v0;
            Vector3 edge2 = v2 - v0;
            Vector3 h = Vector3.Cross(rayDir, edge2);
            float a = Vector3.Dot(edge1, h);

            if (a > -EPSILON && a < EPSILON)
                return false; // Ray is parallel to triangle

            float f = 1.0f / a;
            Vector3 s = rayOrigin - v0;
            float u = f * Vector3.Dot(s, h);

            if (u < 0.0f || u > 1.0f)
                return false;

            Vector3 q = Vector3.Cross(s, edge1);
            float v = f * Vector3.Dot(rayDir, q);

            if (v < 0.0f || u + v > 1.0f)
                return false;

            // Compute distance to intersection point
            distance = f * Vector3.Dot(edge2, q);
            return distance > EPSILON;
        }

        /// <summary>
        /// Tests ray intersection with a bone's geometry.
        /// </summary>
        private static bool RayIntersectsBone(Vector3 rayOrigin, Vector3 rayDir,
                                               FieldBone bone, Matrix4 boneTransform,
                                               out float minDist)
        {
            minDist = float.MaxValue;
            bool hit = false;

            for (int ri = 0; ri < bone.nResources; ri++)
            {
                var model = bone.fRSDResources[ri].Model;
                if (model.Polys == null) continue;

                // Build model transform (pre-multiply order to match OpenGL)
                // OpenGL does: Push bone matrix, Translate, Rotate, Scale, Draw
                // For row-vectors: v * Scale * Rotate * Translate * Bone
                var openTkQuat = new OpenTK.Mathematics.Quaternion(
                    (float)model.rotationQuaternion.x, (float)model.rotationQuaternion.y,
                    (float)model.rotationQuaternion.z, (float)model.rotationQuaternion.w);

                Matrix4 fullTransform = Matrix4.CreateScale(model.resizeX, model.resizeY, model.resizeZ)
                    * Matrix4.CreateFromQuaternion(openTkQuat)
                    * Matrix4.CreateTranslation(model.repositionX, model.repositionY, model.repositionZ)
                    * boneTransform;

                for (int gi = 0; gi < model.Header.numGroups; gi++)
                {
                    if (model.Groups[gi].HiddenQ) continue;

                    int offsetVert = model.Groups[gi].offsetVert;

                    for (int pi = model.Groups[gi].offsetPoly;
                         pi < model.Groups[gi].offsetPoly + model.Groups[gi].numPoly;
                         pi++)
                    {
                        // Transform vertices by the full transform
                        Vector4 v0h = new Vector4(
                            model.Verts[model.Polys[pi].Verts[0] + offsetVert].x,
                            model.Verts[model.Polys[pi].Verts[0] + offsetVert].y,
                            model.Verts[model.Polys[pi].Verts[0] + offsetVert].z, 1.0f) * fullTransform;
                        Vector4 v1h = new Vector4(
                            model.Verts[model.Polys[pi].Verts[1] + offsetVert].x,
                            model.Verts[model.Polys[pi].Verts[1] + offsetVert].y,
                            model.Verts[model.Polys[pi].Verts[1] + offsetVert].z, 1.0f) * fullTransform;
                        Vector4 v2h = new Vector4(
                            model.Verts[model.Polys[pi].Verts[2] + offsetVert].x,
                            model.Verts[model.Polys[pi].Verts[2] + offsetVert].y,
                            model.Verts[model.Polys[pi].Verts[2] + offsetVert].z, 1.0f) * fullTransform;

                        Vector3 v0 = v0h.Xyz / v0h.W;
                        Vector3 v1 = v1h.Xyz / v1h.W;
                        Vector3 v2 = v2h.Xyz / v2h.W;

                        if (RayTriangleIntersect(rayOrigin, rayDir, v0, v1, v2, out float dist))
                        {
                            if (dist > 0 && dist < minDist)
                            {
                                minDist = dist;
                                hit = true;
                            }
                        }
                    }
                }
            }

            return hit;
        }

        public static int GetClosestFieldBone(FieldSkeleton fSkeleton, FieldFrame fFrame,
                                              int px, int py)
        {
            // Get viewport
            int[] vp = new int[4];
            GL.GetInteger(GetPName.Viewport, vp);
            int height = vp[3];

            // Get view and projection matrices from GLRenderer
            Matrix4 view = GLRenderer.ViewMatrix;
            Matrix4 projection = GLRenderer.ProjectionMatrix;

            // Create ray from screen coordinates
            Vector4 viewport = new Vector4(vp[0], vp[1], vp[2], vp[3]);
            float screenY = height - py;

            // Build skeleton root transform
            // Use pre-multiplication to match OpenGL's column-vector convention:
            // OpenGL: M = M * T * R means v' = M * T * R * v (R applied first, then T)
            // Row-vector equivalent: v * R_row * T_row requires building R * T
            double[] rot_mat = new double[16];
            BuildRotationMatrixWithQuaternions(fFrame.rootRotationAlpha, fFrame.rootRotationBeta, fFrame.rootRotationGamma, ref rot_mat);
            Matrix4 rotMatrix = DoubleArrayToMatrix4(rot_mat);

            // Build root transform: rotation first, then translations (in reverse order of OpenGL calls)
            Matrix4 rootTransform = rotMatrix
                * Matrix4.CreateTranslation((float)fFrame.rootTranslationX,
                                           (float)(-fFrame.rootTranslationY),
                                           (float)fFrame.rootTranslationZ);

            // Unproject to create ray (using identity for model since we transform vertices manually)
            Vector3 nearPoint = Unproject(new Vector3(px, screenY, 0.0f), Matrix4.Identity, view, projection, viewport);
            Vector3 farPoint = Unproject(new Vector3(px, screenY, 1.0f), Matrix4.Identity, view, projection, viewport);
            Vector3 rayOrigin = nearPoint;
            Vector3 rayDir = Vector3.Normalize(farPoint - nearPoint);

            // Setup matrix stack for hierarchy traversal
            string[] joint_stack = new string[fSkeleton.bones.Count + 1];
            Matrix4[] matrixStack = new Matrix4[fSkeleton.bones.Count + 2];
            int matrixStackPtr = 0;
            int jsp = 0;

            joint_stack[jsp] = fSkeleton.bones[0].joint_f;
            Matrix4 currentMatrix = rootTransform;
            matrixStack[matrixStackPtr++] = currentMatrix;

            int closestBone = -1;
            float closestDist = float.MaxValue;

            for (int bi = 0; bi < fSkeleton.bones.Count; bi++)
            {
                while (!(fSkeleton.bones[bi].joint_f == joint_stack[jsp]) && jsp > 0)
                {
                    matrixStackPtr--;
                    currentMatrix = matrixStack[matrixStackPtr - 1];
                    jsp--;
                }
                matrixStack[matrixStackPtr++] = currentMatrix;

                BuildRotationMatrixWithQuaternions(fFrame.rotations[bi].alpha, fFrame.rotations[bi].beta, fFrame.rotations[bi].gamma, ref rot_mat);
                rotMatrix = DoubleArrayToMatrix4(rot_mat);
                // Pre-multiply to match OpenGL's transform order
                currentMatrix = rotMatrix * currentMatrix;

                // Test ray intersection with this bone
                if (RayIntersectsBone(rayOrigin, rayDir, fSkeleton.bones[bi], currentMatrix, out float dist))
                {
                    if (dist < closestDist)
                    {
                        closestDist = dist;
                        closestBone = bi;
                    }
                }

                // Pre-multiply translation to match OpenGL's GL.Translate behavior
                currentMatrix = Matrix4.CreateTranslation(0, 0, (float)(-fSkeleton.bones[bi].len)) * currentMatrix;
                jsp++;
                joint_stack[jsp] = fSkeleton.bones[bi].joint_i;
            }

            return closestBone;
        }

        /// <summary>
        /// Tests ray intersection with a specific resource's geometry.
        /// </summary>
        private static bool RayIntersectsResource(Vector3 rayOrigin, Vector3 rayDir,
                                                   FieldRSDResource resource, Matrix4 boneTransform,
                                                   out float minDist)
        {
            minDist = float.MaxValue;
            bool hit = false;

            var model = resource.Model;
            if (model.Polys == null) return false;

            // Build model transform (pre-multiply order to match OpenGL)
            // OpenGL does: Push bone matrix, Translate, Rotate, Scale, Draw
            // For row-vectors: v * Scale * Rotate * Translate * Bone
            var openTkQuat = new OpenTK.Mathematics.Quaternion(
                (float)model.rotationQuaternion.x, (float)model.rotationQuaternion.y,
                (float)model.rotationQuaternion.z, (float)model.rotationQuaternion.w);

            Matrix4 fullTransform = Matrix4.CreateScale(model.resizeX, model.resizeY, model.resizeZ)
                * Matrix4.CreateFromQuaternion(openTkQuat)
                * Matrix4.CreateTranslation(model.repositionX, model.repositionY, model.repositionZ)
                * boneTransform;

            for (int gi = 0; gi < model.Header.numGroups; gi++)
            {
                if (model.Groups[gi].HiddenQ) continue;

                int offsetVert = model.Groups[gi].offsetVert;

                for (int pi = model.Groups[gi].offsetPoly;
                     pi < model.Groups[gi].offsetPoly + model.Groups[gi].numPoly;
                     pi++)
                {
                    // Transform vertices by the full transform
                    Vector4 v0h = new Vector4(
                        model.Verts[model.Polys[pi].Verts[0] + offsetVert].x,
                        model.Verts[model.Polys[pi].Verts[0] + offsetVert].y,
                        model.Verts[model.Polys[pi].Verts[0] + offsetVert].z, 1.0f) * fullTransform;
                    Vector4 v1h = new Vector4(
                        model.Verts[model.Polys[pi].Verts[1] + offsetVert].x,
                        model.Verts[model.Polys[pi].Verts[1] + offsetVert].y,
                        model.Verts[model.Polys[pi].Verts[1] + offsetVert].z, 1.0f) * fullTransform;
                    Vector4 v2h = new Vector4(
                        model.Verts[model.Polys[pi].Verts[2] + offsetVert].x,
                        model.Verts[model.Polys[pi].Verts[2] + offsetVert].y,
                        model.Verts[model.Polys[pi].Verts[2] + offsetVert].z, 1.0f) * fullTransform;

                    Vector3 v0 = v0h.Xyz / v0h.W;
                    Vector3 v1 = v1h.Xyz / v1h.W;
                    Vector3 v2 = v2h.Xyz / v2h.W;

                    if (RayTriangleIntersect(rayOrigin, rayDir, v0, v1, v2, out float dist))
                    {
                        if (dist > 0 && dist < minDist)
                        {
                            minDist = dist;
                            hit = true;
                        }
                    }
                }
            }

            return hit;
        }

        public static int GetClosestFieldBonePiece(FieldSkeleton fSkeleton, FieldFrame fFrame,
                                                   int iBoneSelected, int px, int py)
        {
            // Get viewport
            int[] vp = new int[4];
            GL.GetInteger(GetPName.Viewport, vp);
            int height = vp[3];

            // Get view and projection matrices from GLRenderer
            Matrix4 view = GLRenderer.ViewMatrix;
            Matrix4 projection = GLRenderer.ProjectionMatrix;

            // Create ray from screen coordinates
            Vector4 viewport = new Vector4(vp[0], vp[1], vp[2], vp[3]);
            float screenY = height - py;

            // Build skeleton root transform (pre-multiply to match OpenGL)
            double[] rot_mat = new double[16];
            BuildRotationMatrixWithQuaternions(fFrame.rootRotationAlpha, fFrame.rootRotationBeta, fFrame.rootRotationGamma, ref rot_mat);
            Matrix4 rotMatrix = DoubleArrayToMatrix4(rot_mat);

            Matrix4 rootTransform = rotMatrix
                * Matrix4.CreateTranslation((float)fFrame.rootTranslationX,
                                           (float)(-fFrame.rootTranslationY),
                                           (float)fFrame.rootTranslationZ);

            // Navigate to the selected bone through the hierarchy
            string[] joint_stack = new string[fSkeleton.bones.Count + 1];
            Matrix4[] matrixStack = new Matrix4[fSkeleton.bones.Count + 2];
            int matrixStackPtr = 0;
            int jsp = 0;

            joint_stack[jsp] = fSkeleton.bones[0].joint_f;
            Matrix4 currentMatrix = rootTransform;
            matrixStack[matrixStackPtr++] = currentMatrix;

            for (int iBoneIdx = 0; iBoneIdx < iBoneSelected; iBoneIdx++)
            {
                while (!(fSkeleton.bones[iBoneIdx].joint_f == joint_stack[jsp]) && jsp > 0)
                {
                    matrixStackPtr--;
                    currentMatrix = matrixStack[matrixStackPtr - 1];
                    jsp--;
                }
                matrixStack[matrixStackPtr++] = currentMatrix;

                BuildRotationMatrixWithQuaternions(fFrame.rotations[iBoneIdx].alpha,
                                                   fFrame.rotations[iBoneIdx].beta,
                                                   fFrame.rotations[iBoneIdx].gamma, ref rot_mat);
                rotMatrix = DoubleArrayToMatrix4(rot_mat);
                // Pre-multiply to match OpenGL's transform order
                currentMatrix = rotMatrix * currentMatrix;

                // Pre-multiply translation
                currentMatrix = Matrix4.CreateTranslation(0, 0, (float)(-fSkeleton.bones[iBoneIdx].len)) * currentMatrix;

                jsp++;
                joint_stack[jsp] = fSkeleton.bones[iBoneIdx].joint_i;
            }

            // Navigate to selected bone
            while (!(fSkeleton.bones[iBoneSelected].joint_f == joint_stack[jsp]) && jsp > 0)
            {
                matrixStackPtr--;
                currentMatrix = matrixStack[matrixStackPtr - 1];
                jsp--;
            }

            // Apply selected bone's rotation (pre-multiply)
            BuildRotationMatrixWithQuaternions(fFrame.rotations[iBoneSelected].alpha,
                                               fFrame.rotations[iBoneSelected].beta,
                                               fFrame.rotations[iBoneSelected].gamma, ref rot_mat);
            rotMatrix = DoubleArrayToMatrix4(rot_mat);
            currentMatrix = rotMatrix * currentMatrix;

            // Unproject to create ray
            Vector3 nearPoint = Unproject(new Vector3(px, screenY, 0.0f), Matrix4.Identity, view, projection, viewport);
            Vector3 farPoint = Unproject(new Vector3(px, screenY, 1.0f), Matrix4.Identity, view, projection, viewport);
            Vector3 rayOrigin = nearPoint;
            Vector3 rayDir = Vector3.Normalize(farPoint - nearPoint);

            // Test each resource in the bone
            int closestPiece = -1;
            float closestDist = float.MaxValue;

            for (int iPolyIdx = 0; iPolyIdx < fSkeleton.bones[iBoneSelected].nResources; iPolyIdx++)
            {
                if (RayIntersectsResource(rayOrigin, rayDir, fSkeleton.bones[iBoneSelected].fRSDResources[iPolyIdx],
                                          currentMatrix, out float dist))
                {
                    if (dist < closestDist)
                    {
                        closestDist = dist;
                        closestPiece = iPolyIdx;
                    }
                }
            }

            return closestPiece;
        }

        public static void AddFieldBone(ref FieldBone fBone, ref PModel Model)
        {
            FieldRSDResource tmpRSDResource = new FieldRSDResource()
            {
                ID = "@RSD940102",
                textures = new List<TEX>(),
                Model = Model,
            };

            if (fBone.nResources > 0)
            {
                tmpRSDResource.Model.fileName = fBone.fRSDResources[0].Model.fileName.Substring(0, 4).ToUpper() + 
                                                fBone.nResources.ToString() +
                                                ".P";

                tmpRSDResource.res_file = fBone.fRSDResources[0].res_file.ToUpper() + fBone.nResources.ToString();
            }
            else
            {
                tmpRSDResource.numTextures = 0;

                if (tmpRSDResource.res_file == null)
                {
                    tmpRSDResource.res_file = Model.fileName.Substring(0, 4).ToUpper();
                }
                else
                {
                    tmpRSDResource.res_file = tmpRSDResource.res_file.Substring(0, 4).ToUpper();
                }
                
            }

            fBone.fRSDResources.Add(tmpRSDResource);
            fBone.nResources++;
        }

        public static void RemoveFieldBone(ref FieldBone fBone, ref int b_index)
        {
            FieldRSDResource tmpfRSDResource;

            string tmpRSDResourceName, tmpModelName;
            int iRSDCounter;

            tmpModelName = fBone.fRSDResources[0].Model.fileName.Substring(0, 4).ToUpper();
            tmpRSDResourceName = fBone.fRSDResources[0].res_file.Substring(0, 4).ToUpper();

            if (b_index < fBone.nResources)
            {
                fBone.fRSDResources.RemoveAt(b_index);
                fBone.nResources--;
            }


            // Let's assign the correct names.
            for (iRSDCounter = 0; iRSDCounter < fBone.nResources; iRSDCounter++)
            {
                tmpfRSDResource = fBone.fRSDResources[iRSDCounter];

                tmpfRSDResource.Model.fileName = tmpModelName;
                tmpfRSDResource.res_file = tmpRSDResourceName;

                if (iRSDCounter > 0)
                {
                    tmpfRSDResource.Model.fileName += iRSDCounter.ToString();
                    tmpfRSDResource.res_file += iRSDCounter.ToString();
                }

                tmpfRSDResource.Model.fileName += ".P";

                fBone.fRSDResources[iRSDCounter] = tmpfRSDResource;
            }
        }



        //  ---------------------------------------------------------------------------------------------------
        //  ============================================= SAVING ==============================================
        //  ---------------------------------------------------------------------------------------------------
        public static void MergeResources(ref FieldBone fBone)
        {
            int iResourceIdx;
            FieldRSDResource tmpRSDResource;

            for (iResourceIdx = 1; iResourceIdx < fBone.nResources; iResourceIdx++)
            {
                tmpRSDResource = fBone.fRSDResources[0];
                MergeFieldRSDResources(ref tmpRSDResource, fBone.fRSDResources[iResourceIdx]);
                fBone.fRSDResources[0] = tmpRSDResource;
            }
        }

        public static void ApplyFieldBoneChanges(ref FieldBone fBone, bool merge)
        {
            int ri;
            FieldRSDResource tmpRSDResource;

            for (ri = 0; ri < fBone.nResources; ri++)
            {
                //  Debug.Print "File=", bone.Resources(ri).res_file, bone.Resources(ri).Model.fileName
                if (GL.IsEnabled(EnableCap.Lighting))
                {
                    tmpRSDResource = fBone.fRSDResources[ri];
                    ApplyCurrentVColors(ref tmpRSDResource.Model);
                    fBone.fRSDResources[ri] = tmpRSDResource;
                }

                GL.MatrixMode(MatrixMode.Modelview);
                GL.PushMatrix();

                SetCameraModelViewQuat(fBone.fRSDResources[ri].Model.repositionX, fBone.fRSDResources[ri].Model.repositionY, fBone.fRSDResources[ri].Model.repositionZ,
                                       fBone.fRSDResources[ri].Model.rotationQuaternion,
                                       fBone.fRSDResources[ri].Model.resizeX, fBone.fRSDResources[ri].Model.resizeY, fBone.fRSDResources[ri].Model.resizeZ);

                GL.Scale(fBone.resizeX, fBone.resizeY, fBone.resizeZ);

                tmpRSDResource = fBone.fRSDResources[ri];
                ApplyPChanges(ref tmpRSDResource.Model, false);
                fBone.fRSDResources[ri] = tmpRSDResource;

                GL.MatrixMode(MatrixMode.Modelview);
                GL.PopMatrix();
            }

            if (merge)
            {
                MergeResources(ref fBone);

                if (fBone.nResources > 1)
                {
                    fBone.nResources = 1;
                    while (fBone.fRSDResources.Count > 1) fBone.fRSDResources.RemoveAt(fBone.fRSDResources.Count - 1);

                    tmpRSDResource = fBone.fRSDResources[0];
                    ComputeBoundingBox(ref tmpRSDResource.Model);
                    fBone.fRSDResources[0] = tmpRSDResource;
                }
            }
        }

        public static void ApplyFieldChanges(ref FieldSkeleton fSkeleton, FieldFrame fFrame, bool merge)
        {
            int bi, jsp;
            string[] joint_stack = new string[fSkeleton.bones.Count + 1];

            FieldBone tmpfBone;

            jsp = 0;
            joint_stack[jsp] = fSkeleton.bones[0].joint_f;

            GL.MatrixMode(MatrixMode.Modelview);

            for (bi = 0; bi < fSkeleton.bones.Count; bi++)
            {
                while ((fSkeleton.bones[bi].joint_f != joint_stack[jsp]) && jsp > 0)
                {
                    GL.PopMatrix();
                    jsp--;
                }
                GL.PushMatrix();

                GL.Rotate(fFrame.rotations[bi].beta, 0, 1, 0);
                GL.Rotate(fFrame.rotations[bi].alpha, 1, 0, 0);
                GL.Rotate(fFrame.rotations[bi].gamma, 0, 0, 1);

                tmpfBone = fSkeleton.bones[bi];
                ApplyFieldBoneChanges(ref tmpfBone, merge);
                fSkeleton.bones[bi] = tmpfBone;

                GL.Translate(0, 0, -fSkeleton.bones[bi].len);

                jsp++;
                joint_stack[jsp] = fSkeleton.bones[bi].joint_i;
            }

            while (jsp > 0)
            {
                GL.PopMatrix();
                jsp--;
            }
        }

        public static void WriteFieldBone(ref StringBuilder strHRCContent, ref FieldBone fBone,
                                          string strDirectoryPath)
        {
            int ri;
            string strRSDList;

            FieldRSDResource tmpRSDResource;

            strHRCContent.AppendLine("");
            strHRCContent.AppendLine(fBone.joint_i);
            strHRCContent.AppendLine(fBone.joint_f);
            strHRCContent.AppendLine(fBone.len.ToString("0.0######", CultureInfo.InvariantCulture));

            strRSDList = fBone.nResources.ToString();

            if (fBone.nResources > 0)
            {
                // Write resources (if there is number involved it begins with 1, if we use "0" as first number there are issues).
                //                  first resource has no number.
                for (ri = 0; ri < fBone.nResources; ri++)
                {
                    tmpRSDResource = fBone.fRSDResources[ri];

                    strRSDList = strRSDList + " " + tmpRSDResource.res_file.ToUpper();

                    WriteRSDResource(tmpRSDResource, strDirectoryPath + "\\" + fBone.fRSDResources[ri].res_file.ToUpper() + ".RSD");

                    if (tmpRSDResource.Model.Polys != null)
                        WriteGlobalPModel(ref tmpRSDResource.Model, strDirectoryPath + "\\" + fBone.fRSDResources[ri].Model.fileName.ToUpper());

                    fBone.fRSDResources[ri] = tmpRSDResource;
                }
            }
            else strRSDList += " ";

            strHRCContent.AppendLine(strRSDList);
        }

        public static void WriteFieldSkeleton(ref FieldSkeleton fSkeleton, string fileName)
        {
            int bi;
            StringBuilder strHRCContent = new StringBuilder();
            FieldBone tmpfBone;

            strHRCContent.AppendLine(":HEADER_BLOCK 2");
            strHRCContent.AppendLine(":SKELETON " + fSkeleton.name);
            strHRCContent.AppendLine(":BONES " + fSkeleton.nBones);

            //for (bi = 0; bi < fSkeleton.nBones; bi++)
            for (bi = 0; bi < fSkeleton.bones.Count; bi++)
            {
                tmpfBone = fSkeleton.bones[bi];
                WriteFieldBone(ref strHRCContent, ref tmpfBone, Path.GetDirectoryName(fileName));
                fSkeleton.bones[bi] = tmpfBone;
            }

            File.WriteAllText(fileName.ToUpper(), strHRCContent.ToString());
        }

        public static void CreateDListsFromFieldSkeletonBone(ref FieldBone fBone)
        {
            int ri;
            FieldRSDResource tmpRSDResource;

            for (ri = 0; ri < fBone.nResources; ri++)
            {
                tmpRSDResource = fBone.fRSDResources[ri];
                CreateDListsFromRSDResource(ref tmpRSDResource);
                fBone.fRSDResources[ri] = tmpRSDResource;
            }
        }

        public static void CreateDListsFromFieldSkeleton(ref FieldSkeleton fSkeleton)
        {
            int bi;
            FieldBone tmpfBone;

            //  THIS DOES NOT NEED fAnimation.nBonesCount CHANGE
            for (bi = 0; bi < fSkeleton.bones.Count; bi++)
            {
                tmpfBone = fSkeleton.bones[bi];
                CreateDListsFromFieldSkeletonBone(ref tmpfBone);
                fSkeleton.bones[bi] = tmpfBone;
            }
        }



        //  ---------------------------------------------------------------------------------------------------
        //  ============================================= DESTROY =============================================
        //  ---------------------------------------------------------------------------------------------------
        public static void DestroyFieldBoneResources(ref FieldBone fBone)
        {
            int ri;
            FieldRSDResource tmpRSDResource;

            for (ri = 0; ri < fBone.nResources; ri++)
            {
                tmpRSDResource = fBone.fRSDResources[ri];
                DestroyRSDResources(ref tmpRSDResource);
                fBone.fRSDResources[ri] = tmpRSDResource;
            }

            if (fBone.fRSDResources != null) fBone.fRSDResources.Clear();
        }

        public static void DestroyFieldSkeleton(FieldSkeleton fSkeleton)
        {
            int bi;
            FieldBone tmpfBone;

            //for (bi = 0; bi < fSkeleton.nBones; bi++)
            if (fSkeleton.name != null)
            {
                //  THIS DOES NOT NEED fAnimation.nBonesCount CHANGE
                for (bi = 0; bi < fSkeleton.bones.Count; bi++)
                {
                    tmpfBone = fSkeleton.bones[bi];
                    DestroyFieldBoneResources(ref tmpfBone);
                    fSkeleton.bones[bi] = tmpfBone;
                }

                if (fSkeleton.bones != null) fSkeleton.bones.Clear();
            }
        }



        //  ---------------------------------------------------------------------------------------------------
        //  ========================================== COPY SKELETON ==========================================
        //  ---------------------------------------------------------------------------------------------------
        public static FieldRSDResource CopyfRSDResource(FieldRSDResource fResourceIn)
        {
            FieldRSDResource fResourceOut;

            fResourceOut = new FieldRSDResource()
            {
                ID = fResourceIn.ID,
                numTextures = fResourceIn.numTextures,
                res_file = fResourceIn.res_file,

                textures = new List<TEX>(fResourceIn.textures),

                Model = CopyPModel(fResourceIn.Model),
            };

            foreach (TEX itmTex in fResourceIn.textures) fResourceOut.textures.Add(itmTex);

            return fResourceOut;
        }

        public static FieldBone CopyfBone(FieldBone fBoneIn)
        {
            FieldBone fBoneOut;

            fBoneOut = new FieldBone()
            {
                joint_f = fBoneIn.joint_f,
                joint_i = fBoneIn.joint_i,
                len = fBoneIn.len,
                nResources = fBoneIn.nResources,

                resizeX = fBoneIn.resizeX,
                resizeY = fBoneIn.resizeY,
                resizeZ = fBoneIn.resizeZ,
            };

            fBoneOut.fRSDResources = new List<FieldRSDResource>();

            foreach (FieldRSDResource itmfRSDResource in fBoneIn.fRSDResources) fBoneOut.fRSDResources.Add(CopyfRSDResource(itmfRSDResource));

            return fBoneOut;
        }

        public static FieldSkeleton CopyfSkeleton(FieldSkeleton fSkeletonIn)
        {
            FieldSkeleton fSkeletonOut;

            fSkeletonOut = new FieldSkeleton()
            {
                fileName = fSkeletonIn.fileName,
                name = fSkeletonIn.name,
                nBones = fSkeletonIn.nBones,

                bones = new List<FieldBone>(),
            };

            foreach (FieldBone itmfBone in fSkeletonIn.bones) fSkeletonOut.bones.Add(CopyfBone(itmfBone));

            return fSkeletonOut;
        }




    }
}
