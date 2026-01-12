using System;
using System.Collections.Generic;
using System.Drawing;
using OpenTK.Graphics.OpenGL.Compatibility;
using OpenTK.Mathematics;
using KimeraCS.Rendering;

namespace KimeraCS.Core
{
    using static FF7Skeleton;
    using static FF7FieldSkeleton;
    using static FF7PModel;

    using static FF7BattleSkeleton;

    public static class Utils
    {
        public struct OrderPair
        {
            public float d;
        }

        public struct STIntVector
        {
            public int length;
            public int[] vector;
        }


        //  This is for PEditor
        public struct PairIB
        {
            public int I;
            public float B;
        }

        public const double PIOVER180 = Math.PI / 180;
        public const double QUAT_NORM_TOLERANCE = 0.00001;

        public const double EulRepYes = 1;
        public const double EulParOdd = 1;
        public const double EulFrmR = 1;
        public const float FLT_EPSILON = 1.192092896e-07f;
        public const float MAX_DELTA_SQUARED = 0.001f * 0.001f;

        public const double PI_180 = Math.PI / 180;

        //private int[] Onbits = new int[32];

        public static string strGlobalExceptionMessage;

        // Helper Functions
        //public static bool IsNumeric(string val) => int.TryParse(val, out int _);

        public static void BuildQuaternionFromAxis(ref Vector3 vec, double angle, ref Quaterniond res_quat)
        {
            double sinAngle;
            angle = angle * PIOVER180 / 2;

            sinAngle = Math.Sin(angle);

            res_quat.X = vec.X * sinAngle;
            res_quat.Y = vec.Y * sinAngle;
            res_quat.Z = vec.Z * sinAngle;
            res_quat.W = Math.Cos(angle);
        }

        public static void MultiplyQuaternions(Quaterniond quat_a, Quaterniond quat_b, ref Quaterniond res_quat)
        {
            res_quat.X = quat_a.W * quat_b.X + quat_a.X * quat_b.W + quat_a.Y * quat_b.Z - quat_a.Z * quat_b.Y;
            res_quat.Y = quat_a.W * quat_b.Y + quat_a.Y * quat_b.W + quat_a.Z * quat_b.X - quat_a.X * quat_b.Z;
            res_quat.Z = quat_a.W * quat_b.Z + quat_a.Z * quat_b.W + quat_a.X * quat_b.Y - quat_a.Y * quat_b.X;
            res_quat.W = quat_a.W * quat_b.W - quat_a.X * quat_b.X - quat_a.Y * quat_b.Y - quat_a.Z * quat_b.Z;
        }

        //  Convert Quaternion to Matrix
        public static void BuildMatrixFromQuaternion(Quaterniond quat, ref double[] mat_res)
        {
            double x2, y2, z2;
            double xy, xz, yz;
            double wx, wy, wz;

            x2 = quat.X * quat.X;
            y2 = quat.Y * quat.Y;
            z2 = quat.Z * quat.Z;

            xy = quat.X * quat.Y;
            xz = quat.X * quat.Z;
            yz = quat.Y * quat.Z;

            wx = quat.W * quat.X;
            wy = quat.W * quat.Y;
            wz = quat.W * quat.Z;

            //  This calculation would be a lot more complicated for non-unit length quaternions
            //  Note: The constructor of Matrix4 expects the Matrix in column-major format like expected by
            //  OpenGL
            mat_res[0] = 1 - 2 * (y2 + z2);
            mat_res[4] = 2 * (xy - wz);
            mat_res[8] = 2 * (xz + wy);
            mat_res[12] = 0;
            mat_res[1] = 2 * (xy + wz);
            mat_res[5] = 1 - 2 * (x2 + z2);
            mat_res[9] = 2 * (yz - wx);
            mat_res[13] = 0;
            mat_res[2] = 2 * (xz - wy);
            mat_res[6] = 2 * (yz + wx);
            mat_res[10] = 1 - 2 * (x2 + y2);
            mat_res[14] = 0;
            mat_res[3] = 0;
            mat_res[7] = 0;
            mat_res[11] = 0;
            mat_res[15] = 1;
        }

        public static void BuildRotationMatrixWithQuaternions(double alpha, double beta, double gamma, ref double[] mat_res)
        {
            Quaterniond quat_x = new Quaterniond();
            Quaterniond quat_y = new Quaterniond();
            Quaterniond quat_z = new Quaterniond();
            Quaterniond quat_xy = new Quaterniond();
            Quaterniond quat_xyz = new Quaterniond();

            Vector3 px = new Vector3(1, 0, 0);
            Vector3 py = new Vector3(0, 1, 0);
            Vector3 pz = new Vector3(0, 0, 1);

            BuildQuaternionFromAxis(ref px, alpha, ref quat_x);
            BuildQuaternionFromAxis(ref py, beta, ref quat_y);
            BuildQuaternionFromAxis(ref pz, gamma, ref quat_z);

            MultiplyQuaternions(quat_y, quat_x, ref quat_xy);
            MultiplyQuaternions(quat_xy, quat_z, ref quat_xyz);

            BuildMatrixFromQuaternion(quat_xyz, ref mat_res);
        }


        public static void MultiplyPoint3DByOGLMatrix(double[] matA, Vector3 p_in, ref Vector3 p_out)
        {
            p_out.X = (float)(p_in.X * matA[0] + p_in.Y * matA[4] + p_in.Z * matA[8] + matA[12]);
            p_out.Y = (float)(p_in.X * matA[1] + p_in.Y * matA[5] + p_in.Z * matA[9] + matA[13]);
            p_out.Z = (float)(p_in.X * matA[2] + p_in.Y * matA[6] + p_in.Z * matA[10] + matA[14]);
        }

        public static void ComputeTransformedBoxBoundingBox(double[] MV_matrix, ref Vector3 p_min, ref Vector3 p_max,
                                                            ref Vector3 p_min_trans, ref Vector3 p_max_trans)
        {

            Vector3[] box_pointsV = new Vector3[8];
            Vector3 p_aux_trans = new Vector3();
            int iBoxPoints;

            p_max_trans.X = float.NegativeInfinity;
            p_max_trans.Y = float.NegativeInfinity;
            p_max_trans.Z = float.NegativeInfinity;

            p_min_trans.X = float.PositiveInfinity;
            p_min_trans.Y = float.PositiveInfinity;
            p_min_trans.Z = float.PositiveInfinity;

            box_pointsV[0] = p_min;

            box_pointsV[1].X = p_min.X;
            box_pointsV[1].Y = p_min.Y;
            box_pointsV[1].Z = p_max.Z;

            box_pointsV[2].X = p_min.X;
            box_pointsV[2].Y = p_max.Y;
            box_pointsV[2].Z = p_min.Z;

            box_pointsV[3].X = p_min.X;
            box_pointsV[3].Y = p_max.Y;
            box_pointsV[3].Z = p_max.Z;

            box_pointsV[4] = p_max;

            box_pointsV[5].X = p_max.X;
            box_pointsV[5].Y = p_max.Y;
            box_pointsV[5].Z = p_min.Z;

            box_pointsV[6].X = p_max.X;
            box_pointsV[6].Y = p_min.Y;
            box_pointsV[6].Z = p_max.Z;

            box_pointsV[7].X = p_max.X;
            box_pointsV[7].Y = p_min.Y;
            box_pointsV[7].Z = p_min.Z;

            for (iBoxPoints = 0; iBoxPoints < 8; iBoxPoints++)
            {
                MultiplyPoint3DByOGLMatrix(MV_matrix, box_pointsV[iBoxPoints], ref p_aux_trans);

                if (p_max_trans.X < p_aux_trans.X) p_max_trans.X = p_aux_trans.X;
                if (p_max_trans.Y < p_aux_trans.Y) p_max_trans.Y = p_aux_trans.Y;
                if (p_max_trans.Z < p_aux_trans.Z) p_max_trans.Z = p_aux_trans.Z;

                if (p_min_trans.X > p_aux_trans.X) p_min_trans.X = p_aux_trans.X;
                if (p_min_trans.Y > p_aux_trans.Y) p_min_trans.Y = p_aux_trans.Y;
                if (p_min_trans.Z > p_aux_trans.Z) p_min_trans.Z = p_aux_trans.Z;
            }
        }

        public static void BuildRotationMatrixWithQuaternionsXYZ(double alpha, double beta, double gamma, ref double[] mat_res)
        {
            Quaterniond quat_x = new Quaterniond();
            Quaterniond quat_y = new Quaterniond();
            Quaterniond quat_z = new Quaterniond();
            Quaterniond quat_xy = new Quaterniond();
            Quaterniond quat_xyz = new Quaterniond();
            
            Vector3 px = new Vector3(1, 0, 0);
            Vector3 py = new Vector3(0, 1, 0);
            Vector3 pz = new Vector3(0, 0, 1);

            BuildQuaternionFromAxis(ref px, alpha, ref quat_x);
            BuildQuaternionFromAxis(ref py, beta, ref quat_y);
            BuildQuaternionFromAxis(ref pz, gamma, ref quat_z);

            MultiplyQuaternions(quat_x, quat_y, ref quat_xy);
            MultiplyQuaternions(quat_xy, quat_z, ref quat_xyz);

            BuildMatrixFromQuaternion(quat_xyz, ref mat_res);
        }

        public static Quaterniond GetQuaternionFromEulerUniversal(double y, double x, double z, int i, int j, int k, int n, int s, int f)
        {
            double[] a = new double[3];
            double ti, tj, th, ci, cj, ch, si, sj, sh, cc, cs, sc, ss;

            double t;

            Quaterniond quat_GetQuaternionFromEulerUniversalResult = new Quaterniond();

            if (f == EulFrmR)
            {
                t = x;
                x = z;
                z = t;
            }

            if (n == EulParOdd) y = -y;

            ti = x * 0.5;
            tj = y * 0.5;
            th = z * 0.5;
            ci = Math.Cos(ti);
            cj = Math.Cos(tj);
            ch = Math.Cos(th);
            si = Math.Sin(ti);
            sj = Math.Sin(tj);
            sh = Math.Sin(th);
            cc = ci * ch;
            cs = ci * sh;
            sc = si * ch;
            ss = si * sh;

            if (s == EulRepYes)
            {
                a[i] = cj * (cs + sc); // Could speed up with trig identities.
                a[j] = sj * (cc + ss);
                a[k] = sj * (cs - sc);
                quat_GetQuaternionFromEulerUniversalResult.W = cj * (cc - ss);
            }
            else
            {
                a[i] = cj * sc - sj * cs; // Could speed up with trig identities.
                a[j] = cj * ss + sj * cc;
                a[k] = cj * cs - sj * sc;
                quat_GetQuaternionFromEulerUniversalResult.W = cj * cc + sj * ss;
            }

            if (n == EulParOdd) a[j] = -a[j];

            quat_GetQuaternionFromEulerUniversalResult.X = a[0];
            quat_GetQuaternionFromEulerUniversalResult.Y = a[1];
            quat_GetQuaternionFromEulerUniversalResult.Z = a[2];

            return quat_GetQuaternionFromEulerUniversalResult;
        }

        public static double QuaternionsDot(ref Quaterniond q1, ref Quaterniond q2)
        {
            return q1.X * q2.X + q1.Y * q2.Y + q1.Z * q2.Z + q1.W * q2.W;
        }

        public static void NormalizeQuaternion(ref Quaterniond quat)
        {
            // Don't normalize if we don't have to
            double mag, mag2;

            mag2 = quat.W * quat.W + quat.X * quat.X + quat.Y * quat.Y + quat.Z * quat.Z;

            mag = Math.Sqrt(mag2);

            quat.W /= mag;
            quat.X /= mag;
            quat.Y /= mag;
            quat.Z /= mag;

            //        NEW UDPATE vertex2995 fix for Hojo/Heidegger animations (by L@Zar0)
            //        If Abs(mag2 - 1#) > QUAT_NORM_TOLERANCE Then
            //            mag = Sqr(mag2)
            //            .W = .W / mag
            //            .X = .X / mag
            //            .Y = .Y / mag
            //            .Z = .Z / mag
            //        End If

            //        If .W > 1# Then
            //            .W = 1
            //        End If
        }

        public static Quaterniond QuaternionLerp(ref Quaterniond q1, ref Quaterniond q2, double t)
        {
            double one_minus_t;
            Quaterniond quat_QuaternionLerpResult = new Quaterniond();

            one_minus_t = 1f - t;

            quat_QuaternionLerpResult.X = q1.X * one_minus_t + q2.X * t;
            quat_QuaternionLerpResult.Y = q1.Y * one_minus_t + q2.Y * t;
            quat_QuaternionLerpResult.Z = q1.Z * one_minus_t + q2.Z * t;
            quat_QuaternionLerpResult.W = q1.W * one_minus_t + q2.W * t;

            NormalizeQuaternion(ref quat_QuaternionLerpResult);

            return quat_QuaternionLerpResult;
        }

        public static Quaterniond QuaternionSlerp2(ref Quaterniond q1, ref Quaterniond q2, double t)
        {
            Quaterniond q3 = new Quaterniond();
            Quaterniond quat_QuaternionSlerp2Result = new Quaterniond();

            double dot, angle, one_minus_t, sin_angle, sin_angle_by_t, sin_angle_by_one_t;

            dot = QuaternionsDot(ref q1, ref q2);
            //    dot = cos(theta)
            //    if (dot < 0), q1 and q2 are more than 90 degrees apart,
            //    so we can invert one to reduce spinning
            if (dot < 0)
            {
                dot = -dot;
                q3.X = -q2.X;
                q3.Y = -q2.Y;
                q3.Z = -q2.Z;
                q3.W = -q2.W;
            }
            else
            {
                q3.X = q2.X;
                q3.Y = q2.Y;
                q3.Z = q2.Z;
                q3.W = q2.W;
            }

            if (dot < 0.95)
            {
                angle = Math.Acos(dot);
                one_minus_t = 1f - t;
                sin_angle = Math.Sin(angle);
                sin_angle_by_t = Math.Sin(angle * t);
                sin_angle_by_one_t = Math.Sin(angle * one_minus_t);

                quat_QuaternionSlerp2Result.X = ((q1.X * sin_angle_by_one_t) + q3.X * sin_angle_by_t) / sin_angle;
                quat_QuaternionSlerp2Result.Y = ((q1.Y * sin_angle_by_one_t) + q3.Y * sin_angle_by_t) / sin_angle;
                quat_QuaternionSlerp2Result.Z = ((q1.Z * sin_angle_by_one_t) + q3.Z * sin_angle_by_t) / sin_angle;
                quat_QuaternionSlerp2Result.W = ((q1.W * sin_angle_by_one_t) + q3.W * sin_angle_by_t) / sin_angle;
            }
            else
            {
                quat_QuaternionSlerp2Result = QuaternionLerp(ref q1, ref q3, t);
            }

            return quat_QuaternionSlerp2Result;
        }

        public static double DegToRad(double x)
        {
            return x * Math.PI / 180f;
        }

        public static double RadToDeg(double x)
        {
            return x * 180f / Math.PI;
        }

        public static Quaterniond GetQuaternionFromEulerXYZr(double x, double y, double z)
        {
            return GetQuaternionFromEulerUniversal(DegToRad(x), DegToRad(y), DegToRad(z), 2, 1, 0, 1, 0, 1);
        }

        public static Quaterniond GetQuaternionFromEulerYXZr(double x, double y, double z)
        {
            return GetQuaternionFromEulerUniversal(DegToRad(x), DegToRad(y), DegToRad(z), 2, 0, 1, 0, 0, 1);
        }

        public static Vector3 GetEulerFormMatrixUniversal(double[] mat, int i, int j, int k, int n, int s, int f)
        {
            double sy, cy, t;
            Vector3 up3DGetEulerFormMatrixUniversalResult = new Vector3();

            if (s == EulRepYes)
            {
                sy = Math.Sqrt(mat[i + 4 * j] * mat[i + 4 * j] + mat[i + 4 * k] * mat[i + 4 * k]);
                if (sy > 16f * FLT_EPSILON)
                {
                    up3DGetEulerFormMatrixUniversalResult.X = (float)Math.Atan2(mat[i + 4 * j], mat[i + 4 * k]);
                    up3DGetEulerFormMatrixUniversalResult.Y = (float)Math.Atan2(sy, mat[i + 4 * i]);
                    up3DGetEulerFormMatrixUniversalResult.Z = (float)Math.Atan2(mat[j + 4 * i], -mat[k + 4 * i]);
                }
                else
                {
                    up3DGetEulerFormMatrixUniversalResult.X = (float)Math.Atan2(-mat[j + 4 * k], mat[j + 4 * j]);
                    up3DGetEulerFormMatrixUniversalResult.Y = (float)Math.Atan2(sy, mat[i + 4 * i]);
                    up3DGetEulerFormMatrixUniversalResult.Z = 0;
                }
            }
            else
            {
                cy = Math.Sqrt(mat[i + 4 * i] * mat[i + 4 * i] + mat[j + 4 * i] * mat[j + 4 * i]);
                if (cy >16f * FLT_EPSILON)
                {
                    up3DGetEulerFormMatrixUniversalResult.X = (float)Math.Atan2(mat[k + 4 * j], mat[k + 4 * k]);
                    up3DGetEulerFormMatrixUniversalResult.Y = (float)Math.Atan2(-mat[k + 4 * i], cy);
                    up3DGetEulerFormMatrixUniversalResult.Z = (float)Math.Atan2(mat[j + 4 * i], mat[i + 4 * i]);
                }
                else
                {
                    up3DGetEulerFormMatrixUniversalResult.X = (float)Math.Atan2(-mat[j + 4 * k], mat[j + 4 * j]);
                    up3DGetEulerFormMatrixUniversalResult.Y = (float)Math.Atan2(-mat[k + 4 * i], cy);
                    up3DGetEulerFormMatrixUniversalResult.Z = 0;
                }
            }

            if (n == EulParOdd)
            {
                up3DGetEulerFormMatrixUniversalResult.X = -up3DGetEulerFormMatrixUniversalResult.X;
                up3DGetEulerFormMatrixUniversalResult.Y = -up3DGetEulerFormMatrixUniversalResult.Y;
                up3DGetEulerFormMatrixUniversalResult.Z = -up3DGetEulerFormMatrixUniversalResult.Z;
            }

            if (f == EulFrmR)
            {
                t = up3DGetEulerFormMatrixUniversalResult.X;
                up3DGetEulerFormMatrixUniversalResult.X = up3DGetEulerFormMatrixUniversalResult.Z;
                up3DGetEulerFormMatrixUniversalResult.Z = (float)t;
            }

            up3DGetEulerFormMatrixUniversalResult.X = (float)RadToDeg(up3DGetEulerFormMatrixUniversalResult.X);
            up3DGetEulerFormMatrixUniversalResult.Y = (float)RadToDeg(up3DGetEulerFormMatrixUniversalResult.Y);
            up3DGetEulerFormMatrixUniversalResult.Z = (float)RadToDeg(up3DGetEulerFormMatrixUniversalResult.Z);

            return up3DGetEulerFormMatrixUniversalResult;
        }

        public static Vector3 GetEulerXYZrFromMatrix(double[] mat)
        {
            return GetEulerFormMatrixUniversal(mat, 2, 1, 0, 1, 0, 1);
        }

        public static Vector3 GetEulerYXZrFromMatrix(double[] mat)
        {
            return GetEulerFormMatrixUniversal(mat, 2, 0, 1, 0, 0, 1);
        }

        public static Quaterniond GetQuaternionConjugate(ref Quaterniond quat)
        {
            Quaterniond quat_GetQuaternionConjugateResult =
                new Quaterniond(-quat.X, -quat.Y, -quat.Z, quat.W);

            return quat_GetQuaternionConjugateResult;
        }

        //  Convert from Euler Angles
        public static void BuildQuaternionFromEuler(double alpha, double beta, double gamma, ref Quaterniond res_quat)
        {
            //  Basically we create 3 Quaternions, one for pitch, one for yaw, one for roll
            //  and multiply those together.

            Quaterniond quat_x = new Quaterniond();
            Quaterniond quat_y = new Quaterniond();
            Quaterniond quat_z = new Quaterniond();
            Quaterniond quat_xy = new Quaterniond();

            Vector3 px = new Vector3(1, 0, 0);
            Vector3 py = new Vector3(0, 1, 0);
            Vector3 pz = new Vector3(0, 0, 1);

            BuildQuaternionFromAxis(ref px, alpha, ref quat_x);
            BuildQuaternionFromAxis(ref py, beta, ref quat_y);
            BuildQuaternionFromAxis(ref pz, gamma, ref quat_z);

            MultiplyQuaternions(quat_y, quat_x, ref quat_xy);
            MultiplyQuaternions(quat_xy, quat_z, ref res_quat);

            NormalizeQuaternion(ref res_quat);
        }


        ///////////////////////////////////////////
        // Camera things
        public static void ConcatenateCameraModelView(float cX, float cY, float cZ,
                                                      float alpha, float beta, float gamma,
                                                      float rszX, float rszY, float rszZ)
        {
            double[] rot_mat = new double[16];

            GL.MatrixMode(MatrixMode.Modelview);
            //GL.LoadIdentity();
            GL.Translated(cX, cY, cZ);

            BuildRotationMatrixWithQuaternionsXYZ(alpha, beta, gamma, ref rot_mat);

            GL.MultMatrixd(rot_mat);
            GL.Scaled(rszX, rszY, rszZ);
        }

        public static void ConcatenateCameraModelViewQuat(float cX, float cY, float cZ,
                                                          Quaterniond quat,
                                                          float rszX, float rszY, float rszZ)
        {
            double[] rot_mat = new double[16];

            GL.MatrixMode(MatrixMode.Modelview);
            GL.LoadIdentity();
            GL.Translated(cX, cY, cZ);

            BuildMatrixFromQuaternion(quat, ref rot_mat);

            GL.MultMatrixd(rot_mat);
            GL.Scaled(rszX, rszY, rszZ);
        }

        public static void SetCameraModelView(float cX, float cY, float cZ, 
                                              float alpha, float beta, float gamma,
                                              float rszX, float rszY, float rszZ)
        {
            double[] rot_mat = new double[16];

            GL.MatrixMode(MatrixMode.Modelview);
            GL.LoadIdentity();

            GL.Translated(cX, cY, cZ);

            BuildRotationMatrixWithQuaternionsXYZ(alpha, beta, gamma, ref rot_mat);

            GL.MultMatrixd(rot_mat);

            GL.Scaled(rszX, rszY, rszZ);
        }

        public static void SetCameraModelViewQuat(float cX, float cY, float cZ,
                                                  Quaterniond quat,
                                                  float rszX, float rszY, float rszZ)
        {
            double[] rot_mat = new double[16];

            GL.MatrixMode(MatrixMode.Modelview);
            GL.LoadIdentity();

            GL.Translated(cX, cY, cZ);

            BuildMatrixFromQuaternion(quat, ref rot_mat);

            GL.MultMatrixd(rot_mat);

            GL.Scaled(rszX, rszY, rszZ);
        }

        public static void SetCameraPModel(PModel Model, float cX, float cY, float cZ,
                                           float alpha, float beta, float gamma,
                                           float rszX, float rszY, float rszZ)
        {

            Vector3 p_min = new Vector3();
            Vector3 p_max = new Vector3();
            Vector3 center_model, origin;

            int width, height;

            float model_radius, distance_origin, scene_radius;
            int[] vp = new int[4];

            ComputePModelBoundingBox(Model, ref p_min, ref p_max);

            GL.GetInteger(GetPName.Viewport,vp);
            width = vp[2];
            height = vp[3];

            GL.MatrixMode(MatrixMode.Projection);
            GL.LoadIdentity();

            center_model = new Vector3((p_min.X + p_max.X) / 2,
                                       (p_min.Y + p_max.Y) / 2,
                                       (p_min.Z + p_max.Z) / 2);

            origin = new Vector3();

            model_radius = CalculateDistance(p_min, p_max) / 2;
            distance_origin = CalculateDistance(center_model, origin);
            scene_radius = model_radius + distance_origin;
            gluPerspective(60, (float)width / height, Math.Max(0.1, -cZ - scene_radius), Math.Max(0.1, -cZ + scene_radius));

            SetCameraModelView(cX, cY, cZ, alpha, beta, gamma, rszX, rszY, rszZ);

            // Read matrices directly from OpenGL to ensure picking matches rendering exactly
            float[] projArray = new float[16];
            float[] mvArray = new float[16];
            GL.GetFloat(GetPName.ProjectionMatrix, projArray);
            GL.GetFloat(GetPName.ModelviewMatrix, mvArray);

            GLRenderer.ProjectionMatrix = new Matrix4(
                projArray[0], projArray[1], projArray[2], projArray[3],
                projArray[4], projArray[5], projArray[6], projArray[7],
                projArray[8], projArray[9], projArray[10], projArray[11],
                projArray[12], projArray[13], projArray[14], projArray[15]);

            GLRenderer.ViewMatrix = new Matrix4(
                mvArray[0], mvArray[1], mvArray[2], mvArray[3],
                mvArray[4], mvArray[5], mvArray[6], mvArray[7],
                mvArray[8], mvArray[9], mvArray[10], mvArray[11],
                mvArray[12], mvArray[13], mvArray[14], mvArray[15]);

            GLRenderer.ModelMatrix = Matrix4.Identity;
            GLRenderer.ViewPosition = new OpenTK.Mathematics.Vector3(cX, cY, -cZ);
        }


        public static void SetCameraAroundModel(ref Vector3 p_min, ref Vector3 p_max,
                                                float cX, float cY, float cZ,
                                                float alpha, float beta, float gamma,
                                                float rszX, float rszY, float rszZ)
        {
            float width, height;
            float scene_radius;
            int[] vp = new int[4];

            GL.GetInteger(GetPName.Viewport,vp);
            width = vp[2];
            height = vp[3];

            GL.MatrixMode(MatrixMode.Projection);
            GL.LoadIdentity();

            scene_radius = ComputeSceneRadius(p_min, p_max);

            gluPerspective(60, (float)width / height, Math.Max(0.1, -cZ - scene_radius), Math.Max(0.1, -cZ + scene_radius));

            SetCameraModelView(cX, cY, cZ, alpha, beta, gamma, rszX, rszY, rszZ);

            // Read matrices directly from OpenGL to ensure picking matches rendering exactly
            float[] projArray = new float[16];
            float[] mvArray = new float[16];
            GL.GetFloat(GetPName.ProjectionMatrix, projArray);
            GL.GetFloat(GetPName.ModelviewMatrix, mvArray);

            GLRenderer.ProjectionMatrix = new Matrix4(
                projArray[0], projArray[1], projArray[2], projArray[3],
                projArray[4], projArray[5], projArray[6], projArray[7],
                projArray[8], projArray[9], projArray[10], projArray[11],
                projArray[12], projArray[13], projArray[14], projArray[15]);

            GLRenderer.ViewMatrix = new Matrix4(
                mvArray[0], mvArray[1], mvArray[2], mvArray[3],
                mvArray[4], mvArray[5], mvArray[6], mvArray[7],
                mvArray[8], mvArray[9], mvArray[10], mvArray[11],
                mvArray[12], mvArray[13], mvArray[14], mvArray[15]);

            GLRenderer.ModelMatrix = Matrix4.Identity;
            GLRenderer.ViewPosition = new OpenTK.Mathematics.Vector3(cX, cY, -cZ);
        }

        public static void SetCameraAroundModelQuat(ref Vector3 p_min, ref Vector3 p_max,
                                                    float cX, float cY, float cZ,
                                                    Quaterniond quat,
                                                    float rszX, float rszY, float rszZ)
        {
            float width, height;
            float scene_radius;
            int[] vp = new int[4];

            GL.GetInteger(GetPName.Viewport,vp);
            width = vp[2];
            height = vp[3];

            GL.MatrixMode(MatrixMode.Projection);
            GL.LoadIdentity();

            scene_radius = ComputeSceneRadius(p_min, p_max);

            gluPerspective(60, (float)width / height, Math.Max(0.1, -cZ - scene_radius), Math.Max(0.1, -cZ + scene_radius));

            SetCameraModelViewQuat(cX, cY, cZ, quat, rszX, rszY, rszZ);

            // Read matrices directly from OpenGL to ensure picking matches rendering exactly
            float[] projArray = new float[16];
            float[] mvArray = new float[16];
            GL.GetFloat(GetPName.ProjectionMatrix, projArray);
            GL.GetFloat(GetPName.ModelviewMatrix, mvArray);

            GLRenderer.ProjectionMatrix = new Matrix4(
                projArray[0], projArray[1], projArray[2], projArray[3],
                projArray[4], projArray[5], projArray[6], projArray[7],
                projArray[8], projArray[9], projArray[10], projArray[11],
                projArray[12], projArray[13], projArray[14], projArray[15]);

            GLRenderer.ViewMatrix = new Matrix4(
                mvArray[0], mvArray[1], mvArray[2], mvArray[3],
                mvArray[4], mvArray[5], mvArray[6], mvArray[7],
                mvArray[8], mvArray[9], mvArray[10], mvArray[11],
                mvArray[12], mvArray[13], mvArray[14], mvArray[15]);

            GLRenderer.ModelMatrix = Matrix4.Identity;
            GLRenderer.ViewPosition = new OpenTK.Mathematics.Vector3(cX, cY, -cZ);
        }

        public static bool IsCameraUnderGround()
        {
            Vector3 origin = new Vector3();
            Vector3 originTrans = new Vector3();
            double[] MV_matrix = new double[16];

            GL.GetDouble(GetPName.ModelviewMatrix,MV_matrix);

            InvertMatrix(ref MV_matrix);

            MultiplyPoint3DByOGLMatrix(MV_matrix, origin, ref originTrans);

            return originTrans.Y > -1;
        }

        public static void ResetCamera(ref double alpha, ref double beta, ref double gamma,
                                       ref float panX, ref float panY, ref float panZ,
                                       ref double DIST, int animIndex, int currFrame)
        {
            Vector3 p_min = new Vector3();
            Vector3 p_max = new Vector3();

            //int animIndex;

            if (bLoaded)
            {
                switch (modelType)
                {
                    case ModelType.K_HRC_SKELETON:
                        ComputeFieldBoundingBox(fSkeleton, fAnimation.frames[currFrame],
                                                ref p_min, ref p_max);
                        break;

                    case ModelType.K_AA_SKELETON:
                    case ModelType.K_MAGIC_SKELETON:
                        //if (!bSkeleton.IsBattleLocation)
                        //{
                        ComputeBattleBoundingBox(bSkeleton, bAnimationsPack.SkeletonAnimations[animIndex].frames[currFrame],
                                                 ref p_min, ref p_max);
                        //}
                        break;

                    case ModelType.K_P_FIELD_MODEL:
                    case ModelType.K_P_BATTLE_MODEL:
                    case ModelType.K_P_MAGIC_MODEL:
                    case ModelType.K_3DS_MODEL:
                        ComputePModelBoundingBox(fPModel, ref p_min, ref p_max);
                        break;
                }

                alpha = 200;
                beta = 45;
                gamma = 0;
                panX = 0;
                panY = 0;
                panZ = 0;
                DIST = -2 * ComputeSceneRadius(p_min, p_max);
            }
        }



        ///////////////////////////////////////////
        // Geometric
        public static float CalculateLength3D(Vector3 v)
        {
            return (float)Math.Sqrt(v.X * v.X + v.Y * v.Y + v.Z * v.Z);
        }

        public static Vector3 AddPoint3D(Vector3 v1, Vector3 v2)
        {
            return new Vector3(v1.X + v2.X, v1.Y + v2.Y, v1.Z + v2.Z);
        }

        public static Vector3 SubstractPoint3D(Vector3 v1, Vector3 v2)
        {
            return new Vector3(v1.X - v2.X, v1.Y - v2.Y, v1.Z - v2.Z);
        }

        public static float DotProduct3D(Vector3 v1, Vector3 v2)
        {
            return v1.X * v2.X + v1.Y * v2.Y + v1.Z * v2.Z;
        }

        public static Vector3 CrossProduct3D(Vector3 v1, Vector3 v2)
        {
            return new Vector3(v1.Y * v2.Z - v1.Z * v2.Y,
                               v1.Z * v2.X - v1.X * v2.Z,
                               v1.X * v2.Y - v1.Y * v2.X);
        }

        public static Vector3 DividePoint3D(Vector3 v, float fScalar)
        {
            return new Vector3(v.X / fScalar, v.Y / fScalar, v.Z / fScalar);
        }

        public static float CalculateAngle2Vectors3D(Vector3 v1, Vector3 v2)
        {
            double dAngleRadians;
            float fCalculateAngle2Vectors3DResult;

            dAngleRadians = Math.Acos(DotProduct3D(v1, v2) / (CalculateLength3D(v1) * CalculateLength3D(v2)));

            fCalculateAngle2Vectors3DResult = (float)RadToDeg(dAngleRadians);

            if (float.IsNaN(fCalculateAngle2Vectors3DResult) ||
                float.IsInfinity(fCalculateAngle2Vectors3DResult))
                    fCalculateAngle2Vectors3DResult = 0;

            return fCalculateAngle2Vectors3DResult;
        }

        public static float CalculateAreaPoly3D(Vector3 v0, Vector3 v1, Vector3 v2)
        {
            float a = CalculateDistance(v0, v1);
            float b = CalculateDistance(v1, v2);
            float c = CalculateDistance(v2, v0);
            float s = (a + b + c) / 2;

            return (float)Math.Sqrt(s * (s - a) * (s - b) * (s - c));
        }

        public static Vector3 Normalize(Vector3 v)
        {
            float fLength;

            fLength = CalculateLength3D(v);

            return DividePoint3D(v, fLength);


            //fLength = CalculateLength3D(v);

            //if (fLength > 0)
            //{
            //    fLength = 1 / fLength;

            //    return new Point3D(v.X / fLength, v.Y / fLength, v.Z / fLength);
            //}

            //else return new Point3D(0.0f, 0.0f, 0.0f);
        }

        public static float CalculateDistance(Vector3 v0, Vector3 v1)
        {
            float fDeltaX = v1.X - v0.X;
            float fDeltaY = v1.Y - v0.Y;
            float fDeltaZ = v1.Z - v0.Z;

            return (float)Math.Sqrt(fDeltaX * fDeltaX + fDeltaY * fDeltaY + fDeltaZ * fDeltaZ);
        }

        public static float ComputeSceneRadius(Vector3 p_min, Vector3 p_max)
        {
            float model_radius, distance_origin;

            Vector3 center_model = new Vector3((p_min.X + p_max.X) / 2.0f,
                                               (p_min.Y + p_max.Y) / 2.0f,
                                               (p_min.Z + p_max.Z) / 2.0f);

            Vector3 origin = new Vector3(0, 0, 0);

            model_radius = CalculateDistance(p_min, p_max) / 2;
            distance_origin = CalculateDistance(center_model, origin);

            return model_radius + distance_origin;
        }

        public static Vector3 CalculateNormal(Vector3 p1, Vector3 p2, Vector3 p3)
        {
            Vector3 v1, v2;

            v1 = SubstractPoint3D(p2, p1);
            v2 = SubstractPoint3D(p3, p1);

            return CrossProduct3D(v1, v2);
        }

        public static bool ComparePoints3D(Vector3 a, Vector3 b)
        {
            return (a.X == b.X) && (a.Y == b.Y) && (a.Z == b.Z);
        }

        public static Vector3 CalculateCenteroid(Vector3 p1, Vector3 p2, Vector3 p3)
        {
            return new Vector3((p1.X + p2.X + p3.X) / 3.0f,
                               (p1.Y + p2.Y + p3.Y) / 3.0f,
                               (p1.Z + p2.Z + p3.Z) / 3.0f);
        }

         
        ///////////////////////////////////////////
        // Maths
        public static bool CompareLongs(long val1, long val2)
        {
            if ((val1 ^ val2) < 0) return val1 < 0;
            else return val1 > val2;
        }

        public static void GetSubMatrix(double[] mat, int i, int j, ref double[] mat_out)
        {
            int i2, j2, order, pos;

            order = (int)Math.Sqrt(mat.Length);

            mat_out = new double[(int)Math.Pow(order - 1, 2)];

            for (i2 = 0; i2 < order; i2++)
            {
                if (i2 != i)
                {
                    for (j2 = 0; j2 < order; j2++)
                    {
                        if (j2 != j)
                        {
                            pos = i2 + j2 * (order - 1);
                            if (i2 > i) pos--;
                            if (j2 > j) pos = pos - order + 1;
                            mat_out[pos] = mat[i2 + j2 * order];
                        }
                    }
                }
            }
        }

        public static double GetMatrixDeterminant(ref double[] mat)
        {
            double iGetMatrixDeterminantResult = 0;

            int i, order;
            double det_aux;
            double[] mat_aux = null;

            order = (int)Math.Sqrt(mat.Length);

            if (order > 2)
            {
                for (i = 0; i < order; i++)
                {
                    if (mat[i] != 0)
                    {
                        GetSubMatrix(mat, i, 0, ref mat_aux);
                        det_aux = GetMatrixDeterminant(ref mat_aux) * Math.Pow(-1, i) * mat[i];
                        iGetMatrixDeterminantResult += det_aux;
                    }
                }
            }
            else
            {
                iGetMatrixDeterminantResult = mat[0] * mat[3] - mat[1] * mat[2];
            }

            return iGetMatrixDeterminantResult;
        }

        public static void GetAtachedMatrix(double[] mat, ref double[] mat_out)
        {
            int i, j, order;
            double[] mat_aux = null;

            order = (int)Math.Sqrt(mat.Length);

            mat_out = new double[(int)Math.Pow(order, 2)];

            for (i =  0; i < order; i++)
            {
                for (j = 0; j < order; j++)
                {
                    GetSubMatrix(mat, i, j, ref mat_aux);
                    mat_out[i + j * order] = Math.Pow(-1, i + j) * GetMatrixDeterminant(ref mat_aux);
                }
            }
        }

        public static void TransposeMatrix(ref double[] mat)
        {
            int i, j, order;
            double temp;

            order = (int)Math.Sqrt(mat.Length);

            for (i = 0; i < order; i++)
            {
                for (j = 0; j <= i; j++)
                {
                    temp = mat[i * order + j];
                    mat[i * order + j] = mat[i + j * order];
                    mat[i + j * order] = temp;
                }
            }
        }

        public static void InvertMatrix(ref double[] mat)
        {
            int i, j, order;
            double[] mat_aux = null;
            double det;

            order = (int)Math.Sqrt(mat.Length);
            det = GetMatrixDeterminant(ref mat);

            GetAtachedMatrix(mat, ref mat_aux);

            for (i = 0; i < order; i++)
            {
                for (j = 0; j < order; j++)
                {
                    mat[i + j * order] = mat_aux[i + j * order] / det;
                }
            }

            TransposeMatrix(ref mat);
        }

        //  The value is considered unsigned
        public static int GetBitBlockVUnsigned(byte[] vect, int nBits, ref int FBit)
        {
            int iGetBitBlockVUnsignedResult = 0;
            int baseByte, bi, res, nBytes, unalignedByBits, firstAlignedByte, lastAlignedByte, endBits;
            bool isAligned, cleanEnd;

            if (nBits > 0)
            {
                baseByte = FBit / 8;
                unalignedByBits = FBit % 8;

                if (unalignedByBits + nBits > 8)
                {
                    isAligned = (unalignedByBits == 0);

                    endBits = (FBit + nBits) % 8;
                    cleanEnd = (endBits == 0);

                    nBytes = (nBits - (isAligned ? 0 : 8 - unalignedByBits) - (cleanEnd ? 0 : endBits)) / 8 +
                             (isAligned ? 0 : 1) + (cleanEnd ? 0 : 1);
                    lastAlignedByte = nBytes - (cleanEnd ? 0 : 1) - 1;
                    firstAlignedByte = 0;

                    res = 0;
                    //  Unaligned prefix
                    //  Stored at the begining of the byte
                    if (!isAligned)
                    {
                        res = vect[baseByte];
                        res &= (int)(Math.Pow(2, (8 - unalignedByBits)) - 1);
                        firstAlignedByte = 1;
                    }

                    //  Aligned bytes
                    for (bi = firstAlignedByte; bi <= lastAlignedByte; bi++)
                    {
                        res *= 256;
                        res |= vect[baseByte + bi];
                    }

                    //  Sufix
                    //  Stored at the end of the byte
                    if (!cleanEnd)
                    {
                        res *= (int)Math.Pow(2, endBits);
                        res |= ((vect[baseByte + lastAlignedByte + 1]) / (int)(Math.Pow(2, 8 - endBits)) & (int)(Math.Pow(2, endBits) - 1));
                    }
                }
                else
                {
                    res = vect[baseByte];
                    res /= (int)Math.Pow(2, 8 - (unalignedByBits + nBits));
                    res &= (int)(Math.Pow(2, nBits) - 1);
                }

                iGetBitBlockVUnsignedResult = (short)res;

                FBit += nBits;
            }

            return iGetBitBlockVUnsignedResult;
        }

        public static int ExtendSignInteger(int val, int len)
        {
            int iExtendSignIntegerResult, auxRes;

            //KimeraCS VB6 has this lines but they don't seem to have any effect, right?
            //if (len != 12)
            //{
            //    auxRes = auxRes;
            //}

            if ((val & (int)Math.Pow(2, (len - 1))) != 0)
            {
                auxRes = (int)Math.Pow(2, 16) - 1;
                auxRes ^= (int)(Math.Pow(2, len) - 1);
                auxRes |= val;

                iExtendSignIntegerResult = auxRes;
            }
            else
            {
                iExtendSignIntegerResult = val;
            }

            return iExtendSignIntegerResult;
        }

        public static int GetSignExtendedShort(int src, int valLen)
        {
            int iGetSignExtendedShortResult = 0;

            if (valLen > 0)
            {
                if (valLen < 16)
                {
                    iGetSignExtendedShortResult = ExtendSignInteger(src, valLen);
                }
                else
                {
                    iGetSignExtendedShortResult = src;
                }
            }

            return iGetSignExtendedShortResult;
        }

        public static int GetBitBlockV(byte[] vect, int nBits, ref int FBit)
        {
            int tmpValue;

            tmpValue = GetBitBlockVUnsigned(vect, nBits, ref FBit);

            return GetSignExtendedShort(tmpValue, nBits);
        }

        public static void PutBitBlockV(ref byte[] vect, int nBits, ref int FBit, int iValue)
        {
            int bi, baseByte, nBytes, unalignedByBits;
            int firstAlignedByte, lastAlignedByte, endBits, tmpValue;
            bool isAligned, cleanEnd;

            //  Deal with it as some raw positive value.
            //  Divisions can't be used for bit shifting negative values,
            //  since they round towards 0 instead of minus infinity
            iValue &= (int)(Math.Pow(2, nBits) - 1);

            if (nBits > 0)
            {
                baseByte = FBit / 8;
                unalignedByBits = FBit % 8;

                if (unalignedByBits + nBits > 8)
                {
                    isAligned = (unalignedByBits == 0);

                    endBits = (FBit + nBits) % 8;
                    cleanEnd = (endBits == 0);

                    nBytes = (nBits - (isAligned ? 0 : (8 - unalignedByBits)) - (cleanEnd ? 0 : endBits)) / 8 + 
                             (isAligned ? 0 : 1) + (cleanEnd ? 0 : 1);

                    lastAlignedByte = nBytes - (cleanEnd ? 0 : 1) - 1;
                    firstAlignedByte = 0;

                    Array.Resize(ref vect, baseByte + nBytes);

                    //  Unaligned prefix
                    if (!isAligned)
                    {
                        tmpValue = iValue / (int)(Math.Pow(2, nBits - (8 - unalignedByBits)));
                        tmpValue &= ((int)(Math.Pow(2, (8 - unalignedByBits)) - 1));
                        vect[baseByte] = (byte)(vect[baseByte] | tmpValue);
                        firstAlignedByte = 1;
                    }

                    //  Aligned bytes
                    for (bi = firstAlignedByte; bi <= lastAlignedByte; bi++)
                    {
                        tmpValue = iValue / (int)(Math.Pow(2, ((lastAlignedByte - bi) * 8 + endBits)));
                        vect[baseByte + bi] = (byte)(tmpValue & 0xFF);
                    }

                    // Suffix
                    if (!cleanEnd)
                    {
                        tmpValue = iValue & (int)(Math.Pow(2, endBits) - 1);
                        vect[baseByte + lastAlignedByte + 1] = (byte)(tmpValue * (int)(Math.Pow(2, 8 - endBits)));                            
                    }
                }
                else
                {
                    if (vect.Length - 1 < baseByte)
                    {
                        Array.Resize(ref vect, baseByte + 1);
                        vect[baseByte] = 0;
                    }

                    tmpValue = iValue & (int)Math.Pow(2, nBits) - 1;
                    tmpValue *= (int)Math.Pow(2, 8 - (unalignedByBits + nBits));
                    vect[baseByte] = (byte)(vect[baseByte] | tmpValue);
                }
            }

            FBit += nBits;
        }

        public static float NormalizeAngle180(float fValue)
        {
            float fNormalizeAngle180Result;
            float fDec;

            if (fValue > 0) fDec = 360f;
            else fDec = -360f;

            fNormalizeAngle180Result = fValue;
            while ((fNormalizeAngle180Result > 0 && fValue > 0) || (fNormalizeAngle180Result < 0 && fValue < 0))
                fNormalizeAngle180Result -= fDec;

            if (Math.Abs(fNormalizeAngle180Result) > Math.Abs(fNormalizeAngle180Result + fDec)) fNormalizeAngle180Result += fDec;

            if (fNormalizeAngle180Result >= 180f) fNormalizeAngle180Result -= 360f;

            return fNormalizeAngle180Result;
        }

        public static float GetDegreesFromRaw(int iValue, short key)
        {
            //return (iValue / (float)Math.Pow(2, 12 - key)) * 360;
            float fVal = iValue;
            //fVal = fVal / 4096;
            fVal /= (float)Math.Pow(2, 12 - key);
            fVal *= 360;
            return fVal;
            //return ((float)iValue / 4096) * 360;
        }

        public static int GetRawFromDegrees(float fValue, int key)
        {
            //return (int)((fValue / 360f) * Math.Pow(2, 12 - key));
            float fVal = fValue;
            fVal /= 360;
            //fVal = fVal * 4096;
            fVal *= (float)Math.Pow(2, 12 - key);
            int iVal = (int)Math.Round(fVal);
            return iVal;
            //return (int)(fValue / 360f) * 4096;
        }

        public static int GetBitInteger(int iValue, int iBitIndex)
        {
            return (iValue & (int)Math.Pow(2, iBitIndex)) != 0 ? 1 : 0;
        }

        public static int SetBitInteger(int iValue, int iBitIndex, int iBitValue)
        {
            int iSetBitIntegerResult;

            if (iBitValue == 0) iSetBitIntegerResult = iValue & (~(int)Math.Pow(2, iBitIndex));
            else iSetBitIntegerResult = iValue | ((int)Math.Pow(2, iBitIndex));

            return iSetBitIntegerResult;
        }

        public static int InvertBitInteger(int iValue, int iBitIndex)
        {
            int iInvertBitIntegerResult;

            if (GetBitInteger(iValue, iBitIndex) == 1) iInvertBitIntegerResult = SetBitInteger(iValue, iBitIndex, 0);
            else iInvertBitIntegerResult = SetBitInteger(iValue, iBitIndex, 1);

            return iInvertBitIntegerResult;
        }



        //  -------------------------------------------------------------------------------------------------
        //  ======================================= PEDITOR PROCEDURES ======================================
        //  -------------------------------------------------------------------------------------------------
        public static int GetBrightness(int iR, int iG, int iB)
        {
            return (int)((iR + iG + iB) / 3f);
        }

        public static void FillColorTable(PModel Model, ref List<Color> colorTable,
                                                        ref PairIB[] translationTableVertex, ref PairIB[] translationTablePolys,
                                                        byte iThreshold)
        {
            float v;
            double dv;
            int iC, it, i, iDiff;
            Color cColor;

            colorTable.Clear();

            //colorTable = new Color[Model.Header.numVerts + Model.Header.numPolys];
            translationTablePolys = new PairIB[Model.Header.numPolys];
            translationTableVertex = new PairIB[Model.Header.numVerts];

            for (it = 0; it < Model.Header.numVerts; it++)
            {
                cColor = Model.Vcolors[it];

                v = GetBrightness(cColor.R, cColor.G, cColor.B);
                //  Debug.Print "Brightness(" + Str$(it) + "):" + Str$(v)

                if (v == 0) dv = 255;
                else dv = Math.Round(128 / v, 15);

                iC = -1;
                iDiff = 765;

                for (i = 0; i < colorTable.Count; i++)
                {
                    if (colorTable[i].R <= Math.Min(255, cColor.R + iThreshold) &&
                        colorTable[i].R >= Math.Max(0, cColor.R - iThreshold) &&
                        colorTable[i].G <= Math.Min(255, cColor.G + iThreshold) &&
                        colorTable[i].G >= Math.Max(0, cColor.G - iThreshold) &&
                        colorTable[i].B <= Math.Min(255, cColor.B + iThreshold) &&
                        colorTable[i].B >= Math.Max(0, cColor.B - iThreshold))
                    {
                        if (Math.Abs(cColor.R - colorTable[i].R) +
                            Math.Abs(cColor.G - colorTable[i].G) +
                            Math.Abs(cColor.B - colorTable[i].B) < iDiff)
                        {
                            iDiff = Math.Abs(cColor.R - colorTable[i].R) +
                                    Math.Abs(cColor.G - colorTable[i].G) +
                                    Math.Abs(cColor.B - colorTable[i].B);

                            iC = i;
                        }
                    }
                }

                if (iC == -1)
                {
                    colorTable.Add(Color.FromArgb(255, cColor.R, cColor.G, cColor.B));
                    iC = colorTable.Count - 1;
                }

                translationTableVertex[it].I = iC;
                translationTableVertex[it].B = (float)Math.Round(dv, 7);
            }

            for (it = 0; it < Model.Header.numPolys; it++)
            {
                cColor = Model.Pcolors[it];

                v = GetBrightness(cColor.R, cColor.G, cColor.B);

                if (v == 0) dv = 255;
                else dv = Math.Round(128 / v, 15);

                iC = -1;
                iDiff = 765;

                for (i = 0; i < colorTable.Count; i++)
                {
                    if (colorTable[i].R <= Math.Min(255, cColor.R + iThreshold) &&
                        colorTable[i].R >= Math.Max(0, cColor.R - iThreshold) &&
                        colorTable[i].G <= Math.Min(255, cColor.G + iThreshold) &&
                        colorTable[i].G >= Math.Max(0, cColor.G - iThreshold) &&
                        colorTable[i].B <= Math.Min(255, cColor.B + iThreshold) &&
                        colorTable[i].B >= Math.Max(0, cColor.B - iThreshold))
                    {
                        if (Math.Abs(cColor.R - colorTable[i].R) +
                            Math.Abs(cColor.G - colorTable[i].G) +
                            Math.Abs(cColor.B - colorTable[i].B) < iDiff)
                        {
                            iDiff = Math.Abs(cColor.R - colorTable[i].R) +
                                    Math.Abs(cColor.G - colorTable[i].G) +
                                    Math.Abs(cColor.B - colorTable[i].B);

                            iC = i;
                        }
                    }
                }

                if (iC == -1)
                {
                    colorTable.Add(Color.FromArgb(255, cColor.R, cColor.G, cColor.B));
                    iC = colorTable.Count - 1;
                }

                translationTablePolys[it].I = iC;

                if (dv == 0) translationTablePolys[it].B = 0.001f;
                else translationTablePolys[it].B = (float)Math.Round(dv, 7);
            }
        }

        // -- This function is changed from KimeraVB6. I use direct palette color (no Brightness)
        //public static void FillColorTable(PModel Model, ref List<Color> colorTable, ref int nColors,
        //                                  ref pairIB[] translationTableVertex, ref pairIB[] translationTablePolys,
        //                                  byte iThreshold)
        //{
        //    float v;
        //    double dv;
        //    int tmpR, tmpG, tmpB, iC, it, i, iDiff;
        //    Color cColor;

        //    //colorTable = new Color[Model.Header.numVerts + Model.Header.numPolys];
        //    translationTablePolys = new pairIB[Model.Header.numPolys];
        //    translationTableVertex = new pairIB[Model.Header.numVerts];

        //    for (it = 0; it < Model.Header.numVerts; it++)
        //    {
        //        cColor = Model.Vcolors[it];

        //        v = GetBrightness(cColor.R, cColor.G, cColor.B);
        //        //  Debug.Print "Brightness(" + Str$(it) + "):" + Str$(v)

        //        if (v == 0) dv = 255;
        //        else dv = 128 / v;

        //        tmpR = Math.Min(255, (int)Math.Truncate(cColor.R * dv));
        //        tmpG = Math.Min(255, (int)Math.Truncate(cColor.G * dv));
        //        tmpB = Math.Min(255, (int)Math.Truncate(cColor.B * dv));
        //        iC = -1;
        //        iDiff = 765;

        //        for (i = 0; i < nColors; i++)
        //        {
        //            if ((colorTable[i].R <= Math.Min(255, tmpR + iThreshold) &&
        //                 colorTable[i].R >= Math.Max(0, tmpR - iThreshold)) &&
        //                (colorTable[i].G <= Math.Min(255, tmpG + iThreshold) &&
        //                 colorTable[i].G >= Math.Max(0, tmpG - iThreshold)) &&
        //                (colorTable[i].B <= Math.Min(255, tmpB + iThreshold) &&
        //                 colorTable[i].B >= Math.Max(0, tmpB - iThreshold)))
        //            {
        //                if (Math.Abs(tmpR - colorTable[i].R) +
        //                    Math.Abs(tmpG - colorTable[i].G) +
        //                    Math.Abs(tmpB - colorTable[i].B) < iDiff)
        //                {
        //                    iDiff = Math.Abs(tmpR - colorTable[i].R) +
        //                            Math.Abs(tmpG - colorTable[i].G) +
        //                            Math.Abs(tmpB - colorTable[i].B);

        //                    iC = i;
        //                }
        //            }
        //        }

        //        if (iC == -1)
        //        {
        //            colorTable.Add(Color.FromArgb(255, tmpR, tmpG, tmpB));

        //            nColors++;
        //        }

        //        if (iC == -1) iC = nColors - 1;

        //        translationTableVertex[it].I = iC;
        //        translationTableVertex[it].B = (float)dv;
        //    }

        //    for (it = 0; it < Model.Header.numPolys; it++)
        //    {
        //        cColor = Model.Pcolors[it];

        //        v = GetBrightness(cColor.R, cColor.G, cColor.B);

        //        if (v == 0) dv = 128;
        //        else dv = 128 / v;

        //        tmpR = Math.Min(255, (int)Math.Truncate(cColor.R * dv));
        //        tmpG = Math.Min(255, (int)Math.Truncate(cColor.G * dv));
        //        tmpB = Math.Min(255, (int)Math.Truncate(cColor.B * dv));

        //        iC = -1;
        //        iDiff = 765;

        //        for (i = 0; i < nColors; i++)
        //        {
        //            if ((colorTable[i].R <= Math.Min(255, tmpR + iThreshold) &&
        //                 colorTable[i].R >= Math.Max(0, tmpR - iThreshold)) &&
        //                (colorTable[i].G <= Math.Min(255, tmpG + iThreshold) &&
        //                 colorTable[i].G >= Math.Max(0, tmpG - iThreshold)) &&
        //                (colorTable[i].B <= Math.Min(255, tmpB + iThreshold) &&
        //                 colorTable[i].B >= Math.Max(0, tmpB - iThreshold)))
        //            {
        //                if (Math.Abs(tmpR - colorTable[i].R) +
        //                    Math.Abs(tmpG - colorTable[i].G) +
        //                    Math.Abs(tmpB - colorTable[i].B) < iDiff)
        //                {
        //                    iDiff = Math.Abs(tmpR - colorTable[i].R) +
        //                            Math.Abs(tmpG - colorTable[i].G) +
        //                            Math.Abs(tmpB - colorTable[i].B);

        //                    iC = i;
        //                }
        //            }
        //        }

        //        if (iC == -1)
        //        {
        //            colorTable.Add(Color.FromArgb(255, tmpR, tmpG, tmpB));

        //            iC = nColors - 1;
        //        }

        //        translationTablePolys[it].I = iC;

        //        if (dv == 0) translationTablePolys[it].B = 0.001f;
        //        else translationTablePolys[it].B = (float)dv;
        //    }
        //}

        public static void ApplyColorTable(ref PModel Model, List<Color> colorTable, PairIB[] translationTableVertex,
                                                                                     PairIB[] translationTablePolys)
        {
            int iPolyIdx, iVertIdx;
            Color cColor;

            for (iVertIdx = 0; iVertIdx < Model.Header.numVerts; iVertIdx++)
            {
                cColor = colorTable[translationTableVertex[iVertIdx].I];

                if (!Model.Groups[GetVertexGroup(Model, iVertIdx)].HiddenQ)
                {
                    Model.Vcolors[iVertIdx] = Color.FromArgb(255, cColor.R, cColor.G, cColor.B);
                }
            }

            // -- This is not in KimeraVB6
            for (iPolyIdx = 0; iPolyIdx < Model.Header.numPolys; iPolyIdx++)
            {
                cColor = colorTable[translationTablePolys[iPolyIdx].I];

                if (!Model.Groups[GetPolygonGroup(Model, iPolyIdx)].HiddenQ)
                {
                    Model.Pcolors[iPolyIdx] = Color.FromArgb(255, cColor.R, cColor.G, cColor.B);
                }
            }
        }

        // -- This function is changed from KimeraVB6. I use direct palette color (no Brightness)
        //public static void ApplyColorTable(ref PModel Model, List<Color> colorTable, pairIB[] translationTableVertex,
        //                                                                             pairIB[] translationTablePolys)
        //{
        //    int vi, pi;
        //    Color cColor;
        //    float dv;

        //    for (vi = 0; vi < Model.Header.numVerts; vi++)
        //    {
        //        cColor = colorTable[translationTableVertex[iVertIdx].I];
        //        dv = translationTableVertex[iVertIdx].B;

        //        Model.Vcolors[iVertIdx] = Color.FromArgb(255,
        //                                           (byte)Math.Max(0, Math.Min(255, Math.Ceiling(cColor.R / dv))),
        //                                           (byte)Math.Max(0, Math.Min(255, Math.Ceiling(cColor.G / dv))),
        //                                           (byte)Math.Max(0, Math.Min(255, Math.Ceiling(cColor.B / dv))));
        //    }

        //}

        //public static void ChangeBrightness(ref PModel Model, int iFactor, Color[] vcolorsOriginal, Color[] pcolorsOriginal)
        public static void ChangeBrightness(ref PModel Model, int iFactor, Color[] vcolorsOriginal)
        {
            int iVertIdx;

            for (iVertIdx = 0; iVertIdx < Model.Header.numVerts; iVertIdx++)
            {

                Model.Vcolors[iVertIdx] = Color.FromArgb(vcolorsOriginal[iVertIdx].A,       // 255
                                Math.Max(0, Math.Min(255, vcolorsOriginal[iVertIdx].R + iFactor)),
                                Math.Max(0, Math.Min(255, vcolorsOriginal[iVertIdx].G + iFactor)),
                                Math.Max(0, Math.Min(255, vcolorsOriginal[iVertIdx].B + iFactor)));
            }

            ComputePColors(ref Model);

            //for (iPolyIdx = 0; iPolyIdx < Model.Header.numPolys; iPolyIdx++)
            //{
            //    Model.Pcolors[iPolyIdx] = Color.FromArgb(255,
            //                         Math.Max(0, Math.Min(255, pcolorsOriginal[iPolyIdx].R + iFactor)),
            //                         Math.Max(0, Math.Min(255, pcolorsOriginal[iPolyIdx].G + iFactor)),
            //                         Math.Max(0, Math.Min(255, pcolorsOriginal[iPolyIdx].B + iFactor)));
            //}
        }

        public static void UpdateTranslationTable(ref PairIB[] translationTableVertex, 
                                                  PModel Model, int pIndex, int cIndex)
        {
            int iVertIdx, iGroupIdx, iDiff, baseVert;

            iDiff = Model.Header.numVerts - 1 - (translationTableVertex.Length - 1);

            iGroupIdx = GetPolygonGroup(Model, pIndex);
            baseVert = Model.Groups[iGroupIdx].offsetVert + Model.Groups[iGroupIdx].numVert - 1 - iDiff;

            Array.Resize(ref translationTableVertex, Model.Header.numVerts);

            for (iVertIdx = Model.Header.numVerts - 1; iVertIdx >= baseVert + 1; iVertIdx--)
            {
                translationTableVertex[iVertIdx].I = translationTableVertex[iVertIdx - iDiff].I;
                translationTableVertex[iVertIdx].B = translationTableVertex[iVertIdx - iDiff].B;
            }

            for (iVertIdx = baseVert + 1; iVertIdx <= baseVert + iDiff; iVertIdx++)
            {
                translationTableVertex[iVertIdx].I = cIndex;
                translationTableVertex[iVertIdx].B = 1;
            }
        }

        public static float CalculatePoint2LineProjectionPosition(Vector3 q, Vector3 p1, Vector3 p2)
        {
            float alpha;
            Vector3 vdUP3D = new Vector3(p2.X - p1.X, p2.Y - p1.Y, p2.Z - p1.Z);

            alpha = (float)((vdUP3D.X * (q.X - p1.X) + vdUP3D.Y * (q.Y - p1.Y) + vdUP3D.Z * (q.Z - p1.Z)) /
                            (Math.Pow(vdUP3D.X, 2) + Math.Pow(vdUP3D.Y, 2) + Math.Pow(vdUP3D.Z, 2)));

            if (alpha > 1) alpha = 1;
            if (alpha < -1) alpha = -1;

            return alpha;
        }

        public static Vector3 CalculatePoint2LineProjection(Vector3 q, Vector3 p1, Vector3 p2)
        {
            float alpha;

            alpha = CalculatePoint2LineProjectionPosition(q, p1, p2);

            return GetPointInLine(p1, p2, alpha);
        }

        // This function will check if there are any duplicated vertices or duplicated polys indexes
        // in Add Polygon feature of PEditor.
        // iArrayVNP = VertexNewPoly        iVCNP = VertexCountNewPoly
        public static bool ValidateAddPolygonVertices(PModel Model, ushort[] iArrayVNP, int iVCNP)
        {
            bool bValidateVerts = true;
            int iGrpv0, iGrpv1, iGrpv2;
            Vector3 p3Dv0, p3Dv1, p3Dv2;

            if (iVCNP > 1)
            {
                switch (iVCNP)
                {
                    case 2:
                        p3Dv0 = Model.Verts[iArrayVNP[0]];
                        p3Dv1 = Model.Verts[iArrayVNP[1]];

                        iGrpv0 = GetVertexGroup(Model, iArrayVNP[0]);
                        iGrpv1 = GetVertexGroup(Model, iArrayVNP[1]);

                        if (iArrayVNP[0] == iArrayVNP[1] || ComparePoints3D(p3Dv0, p3Dv1) || iGrpv0 != iGrpv1)
                            bValidateVerts = false;

                        break;

                    case 3:
                        p3Dv0 = Model.Verts[iArrayVNP[0]];
                        p3Dv1 = Model.Verts[iArrayVNP[1]];
                        p3Dv2 = Model.Verts[iArrayVNP[2]];

                        iGrpv0 = GetVertexGroup(Model, iArrayVNP[0]);
                        iGrpv1 = GetVertexGroup(Model, iArrayVNP[1]);
                        iGrpv2 = GetVertexGroup(Model, iArrayVNP[2]);

                        if (iArrayVNP[0] == iArrayVNP[1] || iArrayVNP[0] == iArrayVNP[2] || iArrayVNP[1] == iArrayVNP[2] ||
                            ComparePoints3D(p3Dv0, p3Dv1) || ComparePoints3D(p3Dv0, p3Dv2) || ComparePoints3D(p3Dv1, p3Dv2) ||
                            iGrpv0 != iGrpv1 || iGrpv0 != iGrpv2)
                            bValidateVerts = false;

                        break;
                }
            }

            return bValidateVerts;
        }



        //  -------------------------------------------------------------------------------------------------
        //  ================================= OPENTK MATRIX UTILITIES =======================================
        //  -------------------------------------------------------------------------------------------------

        /// <summary>
        /// Creates a perspective projection matrix (replaces gluPerspective)
        /// </summary>
        /// <param name="fovDegrees">Field of view in degrees</param>
        /// <param name="aspect">Aspect ratio (width/height)</param>
        /// <param name="near">Near clipping plane</param>
        /// <param name="far">Far clipping plane</param>
        public static Matrix4 CreatePerspectiveMatrix(float fovDegrees, float aspect, float near, float far)
        {
            return Matrix4.CreatePerspectiveFieldOfView(
                MathHelper.DegreesToRadians(fovDegrees),
                aspect,
                near,
                far);
        }

        /// <summary>
        /// Creates a view matrix looking at a target (replaces gluLookAt)
        /// </summary>
        public static Matrix4 CreateLookAtMatrix(OpenTK.Mathematics.Vector3 eye, OpenTK.Mathematics.Vector3 target, OpenTK.Mathematics.Vector3 up)
        {
            return Matrix4.LookAt(eye, target, up);
        }

        /// <summary>
        /// Projects a 3D world coordinate to 2D screen coordinates (replaces gluProject)
        /// </summary>
        /// <param name="worldPos">World position to project</param>
        /// <param name="model">Model matrix</param>
        /// <param name="view">View matrix</param>
        /// <param name="projection">Projection matrix</param>
        /// <param name="viewport">Viewport (x, y, width, height)</param>
        /// <returns>Screen coordinates (x, y, depth)</returns>
        public static OpenTK.Mathematics.Vector3 Project(OpenTK.Mathematics.Vector3 worldPos, Matrix4 model, Matrix4 view, Matrix4 projection, Vector4 viewport)
        {
            Vector4 clipPos = new Vector4(worldPos, 1.0f) * model * view * projection;

            if (Math.Abs(clipPos.W) < float.Epsilon)
                return OpenTK.Mathematics.Vector3.Zero;

            OpenTK.Mathematics.Vector3 ndc = clipPos.Xyz / clipPos.W;

            float winX = viewport.Z * (ndc.X + 1.0f) / 2.0f + viewport.X;
            float winY = viewport.W * (ndc.Y + 1.0f) / 2.0f + viewport.Y;
            float winZ = (ndc.Z + 1.0f) / 2.0f;

            return new OpenTK.Mathematics.Vector3(winX, winY, winZ);
        }

        /// <summary>
        /// Unprojects 2D screen coordinates to 3D world coordinates (replaces gluUnProject)
        /// </summary>
        /// <param name="screenPos">Screen position (x, y, depth)</param>
        /// <param name="model">Model matrix</param>
        /// <param name="view">View matrix</param>
        /// <param name="projection">Projection matrix</param>
        /// <param name="viewport">Viewport (x, y, width, height)</param>
        /// <returns>World coordinates</returns>
        public static OpenTK.Mathematics.Vector3 Unproject(OpenTK.Mathematics.Vector3 screenPos, Matrix4 model, Matrix4 view, Matrix4 projection, Vector4 viewport)
        {
            // Use row-vector convention to match Project function: pos * M * V * P
            Matrix4 mvp = model * view * projection;
            Matrix4 invMvp = mvp.Inverted();

            Vector4 ndc = new Vector4(
                2.0f * (screenPos.X - viewport.X) / viewport.Z - 1.0f,
                2.0f * (screenPos.Y - viewport.Y) / viewport.W - 1.0f,
                2.0f * screenPos.Z - 1.0f,
                1.0f);

            // Row vector multiplication: ndc * invMvp
            Vector4 worldPos = ndc * invMvp;

            if (Math.Abs(worldPos.W) < float.Epsilon)
                return OpenTK.Mathematics.Vector3.Zero;

            return worldPos.Xyz / worldPos.W;
        }

        /// <summary>
        /// Creates a model-view matrix from camera parameters (replaces SetCameraModelView pattern)
        /// </summary>
        public static Matrix4 CreateModelViewMatrix(float cX, float cY, float cZ,
                                                     float alpha, float beta, float gamma,
                                                     float scaleX, float scaleY, float scaleZ)
        {
            // Build rotation from quaternions (matching existing BuildRotationMatrixWithQuaternionsXYZ)
            var quatX = OpenTK.Mathematics.Quaternion.FromAxisAngle(OpenTK.Mathematics.Vector3.UnitX, MathHelper.DegreesToRadians(alpha));
            var quatY = OpenTK.Mathematics.Quaternion.FromAxisAngle(OpenTK.Mathematics.Vector3.UnitY, MathHelper.DegreesToRadians(beta));
            var quatZ = OpenTK.Mathematics.Quaternion.FromAxisAngle(OpenTK.Mathematics.Vector3.UnitZ, MathHelper.DegreesToRadians(gamma));
            var rotation = quatX * quatY * quatZ;

            Matrix4 rotationMatrix = Matrix4.CreateFromQuaternion(rotation);
            Matrix4 translationMatrix = Matrix4.CreateTranslation(cX, cY, cZ);
            Matrix4 scaleMatrix = Matrix4.CreateScale(scaleX, scaleY, scaleZ);

            return translationMatrix * rotationMatrix * scaleMatrix;
        }

        /// <summary>
        /// Creates a model-view matrix from camera parameters using a quaternion for rotation
        /// </summary>
        public static Matrix4 CreateModelViewMatrixQuat(float cX, float cY, float cZ,
                                                         Quaterniond quat,
                                                         float scaleX, float scaleY, float scaleZ)
        {
            var openTkQuat = new OpenTK.Mathematics.Quaternion((float)quat.X, (float)quat.Y, (float)quat.Z, (float)quat.W);

            Matrix4 rotationMatrix = Matrix4.CreateFromQuaternion(openTkQuat);
            Matrix4 translationMatrix = Matrix4.CreateTranslation(cX, cY, cZ);
            Matrix4 scaleMatrix = Matrix4.CreateScale(scaleX, scaleY, scaleZ);

            return translationMatrix * rotationMatrix * scaleMatrix;
        }

        /// <summary>
        /// Converts the double[] matrix (16 elements, column-major) to OpenTK Matrix4
        /// </summary>
        public static Matrix4 ToMatrix4(double[] mat)
        {
            return new Matrix4(
                (float)mat[0], (float)mat[1], (float)mat[2], (float)mat[3],
                (float)mat[4], (float)mat[5], (float)mat[6], (float)mat[7],
                (float)mat[8], (float)mat[9], (float)mat[10], (float)mat[11],
                (float)mat[12], (float)mat[13], (float)mat[14], (float)mat[15]);
        }

        /// <summary>
        /// Converts OpenTK Matrix4 to double[] (16 elements, column-major)
        /// </summary>
        public static double[] FromMatrix4(Matrix4 mat)
        {
            return new double[]
            {
                mat.M11, mat.M12, mat.M13, mat.M14,
                mat.M21, mat.M22, mat.M23, mat.M24,
                mat.M31, mat.M32, mat.M33, mat.M34,
                mat.M41, mat.M42, mat.M43, mat.M44
            };
        }

        /// <summary>
        /// Alias for ToMatrix4 - converts double[] matrix to OpenTK Matrix4
        /// </summary>
        public static Matrix4 DoubleArrayToMatrix4(double[] mat)
        {
            return ToMatrix4(mat);
        }


        //  -------------------------------------------------------------------------------------------------
        //  ================================= OPENGL LEGACY HELPER FUNCTIONS ================================
        //  -------------------------------------------------------------------------------------------------

        /// <summary>
        /// Set blend mode (legacy helper) - uses Defines.BLEND_MODE enum
        /// </summary>
        public static void SetBlendMode(BlendMode bmMode)
        {
            if (bmMode == BlendMode.Disabled)
            {
                GL.Disable(EnableCap.Blend);
            }
            else
            {
                GL.Enable(EnableCap.Blend);
                GL.BlendEquation(BlendEquationMode.FuncAdd);

                switch (bmMode)
                {
                    case BlendMode.Average:
                        GL.BlendFunc(BlendingFactor.SrcAlpha, BlendingFactor.OneMinusSrcAlpha);
                        break;

                    case BlendMode.Add:
                        GL.BlendFunc(BlendingFactor.One, BlendingFactor.One);
                        break;

                    case BlendMode.Subtract:
                        GL.BlendFunc(BlendingFactor.One, BlendingFactor.One);
                        GL.BlendEquation(BlendEquationMode.FuncReverseSubtract);
                        break;

                    case BlendMode._25P:
                        GL.BlendFunc(BlendingFactor.SrcAlpha, BlendingFactor.One);
                        break;

                    case BlendMode.None:
                        GL.BlendFunc(BlendingFactor.One, BlendingFactor.Zero);
                        break;
                }
            }
        }

        /// <summary>
        /// Clear the OpenGL panel with default background color
        /// </summary>
        public static void ClearPanel()
        {
            GL.ClearColor(0.4f, 0.4f, 0.65f, 0);
            GL.Clear(ClearBufferMask.ColorBufferBit | ClearBufferMask.DepthBufferBit);
        }

        /// <summary>
        /// Set default OpenGL render state
        /// </summary>
        public static void SetDefaultOGLRenderState()
        {
            GL.PolygonMode(TriangleFace.FrontAndBack, PolygonMode.Fill);
            GL.CullFace(TriangleFace.Front);
            GL.Enable(EnableCap.CullFace);

            GL.DepthFunc(DepthFunction.Lequal);
            GL.Enable(EnableCap.DepthTest);
            GL.DepthMask(true);

            GL.TexParameteri(TextureTarget.Texture2d, TextureParameterName.TextureMinFilter, (int)TextureMinFilter.Linear);
            GL.TexParameteri(TextureTarget.Texture2d, TextureParameterName.TextureMagFilter, (int)TextureMagFilter.Linear);

            SetBlendMode(BlendMode.Disabled);
        }

        /// <summary>
        /// GLU perspective wrapper - also loads the matrix into the current matrix mode
        /// </summary>
        public static void gluPerspective(double fov, double aspect, double zNear, double zFar)
        {
            Matrix4 perspectiveMatrix = CreatePerspectiveMatrix((float)fov, (float)aspect, (float)zNear, (float)zFar);
            double[] matArray = Matrix4ToDoubleArray(perspectiveMatrix);
            GL.MultMatrixd(matArray);
        }

        /// <summary>
        /// GLU ortho 2D wrapper - sets up 2D orthographic projection
        /// </summary>
        public static void gluOrtho2D(double left, double right, double bottom, double top)
        {
            Matrix4 orthoMatrix = Matrix4.CreateOrthographicOffCenter((float)left, (float)right, (float)bottom, (float)top, -1, 1);
            double[] matArray = Matrix4ToDoubleArray(orthoMatrix);
            GL.MultMatrixd(matArray);
        }

        /// <summary>
        /// Get projected coordinates (world to screen)
        /// </summary>
        public static Vector3 GetProjectedCoords(Vector3 p)
        {
            float[] mm = new float[16];
            float[] pm = new float[16];
            int[] vp = new int[4];

            GL.GetFloat(GetPName.ModelviewMatrix, mm);
            GL.GetFloat(GetPName.ProjectionMatrix, pm);
            GL.GetInteger(GetPName.Viewport, vp);

            Matrix4 modelView = new Matrix4(
                mm[0], mm[1], mm[2], mm[3],
                mm[4], mm[5], mm[6], mm[7],
                mm[8], mm[9], mm[10], mm[11],
                mm[12], mm[13], mm[14], mm[15]);

            Matrix4 projection = new Matrix4(
                pm[0], pm[1], pm[2], pm[3],
                pm[4], pm[5], pm[6], pm[7],
                pm[8], pm[9], pm[10], pm[11],
                pm[12], pm[13], pm[14], pm[15]);

            Vector4 viewport = new Vector4(vp[0], vp[1], vp[2], vp[3]);
            Vector3 result = Project(p, Matrix4.Identity, modelView, projection, viewport);

            return new Vector3(result.X, result.Y, result.Z);
        }

        /// <summary>
        /// Get unprojected coordinates (screen to world)
        /// </summary>
        public static Vector3 GetUnProjectedCoords(Vector3 p)
        {
            float[] mm = new float[16];
            float[] pm = new float[16];
            int[] vp = new int[4];

            GL.GetFloat(GetPName.ModelviewMatrix, mm);
            GL.GetFloat(GetPName.ProjectionMatrix, pm);
            GL.GetInteger(GetPName.Viewport, vp);

            Matrix4 modelView = new Matrix4(
                mm[0], mm[1], mm[2], mm[3],
                mm[4], mm[5], mm[6], mm[7],
                mm[8], mm[9], mm[10], mm[11],
                mm[12], mm[13], mm[14], mm[15]);

            Matrix4 projection = new Matrix4(
                pm[0], pm[1], pm[2], pm[3],
                pm[4], pm[5], pm[6], pm[7],
                pm[8], pm[9], pm[10], pm[11],
                pm[12], pm[13], pm[14], pm[15]);

            Vector4 viewport = new Vector4(vp[0], vp[1], vp[2], vp[3]);
            // Note: Y coordinate is flipped in screen space
            OpenTK.Mathematics.Vector3 screenPos = new OpenTK.Mathematics.Vector3(p.X, vp[3] - p.Y, p.Z);
            OpenTK.Mathematics.Vector3 result = Unproject(screenPos, Matrix4.Identity, modelView, projection, viewport);

            return new Vector3(result.X, result.Y, result.Z);
        }

        /// <summary>
        /// Get projected vertex coordinates
        /// </summary>
        public static Vector3 GetVertexProjectedCoords(Vector3[] lstVerts, int iVertIdx)
        {
            GL.Clear(ClearBufferMask.DepthBufferBit);
            return GetProjectedCoords(lstVerts[iVertIdx]);
        }

        /// <summary>
        /// Get depth Z of a point
        /// </summary>
        public static float GetDepthZ(Vector3 pUP3D)
        {
            return GetProjectedCoords(pUP3D).Z;
        }

        /// <summary>
        /// Get eye space coordinates
        /// </summary>
        public static Vector3 GetEyeSpaceCoords(Vector3 p)
        {
            float[] mm = new float[16];
            GL.GetFloat(GetPName.ModelviewMatrix, mm);

            return new Vector3(
                p.X * mm[0] + p.Y * mm[4] + p.Z * mm[8] + mm[12],
                p.X * mm[1] + p.Y * mm[5] + p.Z * mm[9] + mm[13],
                p.X * mm[2] + p.Y * mm[6] + p.Z * mm[10] + mm[14]);
        }

        /// <summary>
        /// Get vertex color with lighting applied
        /// </summary>
        public static Color GetVertColor(Vector3 p, Vector3 n, Color c)
        {
            byte[] pcolor = new byte[4];
            int[] vp0 = new int[4];
            int[] vp = new int[4];

            GL.GetInteger(GetPName.Viewport, vp0);
            GL.Viewport(0, 0, 3, 3);

            GL.GetInteger(GetPName.Viewport, vp);
            GL.MatrixMode(MatrixMode.Projection);
            GL.PushMatrix();

            GL.Clear(ClearBufferMask.ColorBufferBit | ClearBufferMask.DepthBufferBit);

            GL.PointSize(100.0f);

            GL.Begin(PrimitiveType.Points);
            GL.Color4f(c.R / 255.0f, c.G / 255.0f, c.B / 255.0f, 1.0f);
            GL.ColorMaterial(TriangleFace.FrontAndBack, ColorMaterialParameter.AmbientAndDiffuse);

            GL.Normal3f(n.X, n.Y, n.Z);
            GL.Vertex3f(p.X, p.Y, p.Z);
            GL.End();

            GL.Flush();
            GL.ReadBuffer(ReadBufferMode.Back);

            GL.ReadPixels(1, 1, 1, 1, OpenTK.Graphics.OpenGL.Compatibility.PixelFormat.Rgb, PixelType.UnsignedByte, pcolor);

            Color result = Color.FromArgb(255, pcolor[0] * 2, pcolor[1] * 2, pcolor[2] * 2);

            GL.PopMatrix();
            GL.Viewport(vp0[0], vp0[1], vp0[2], vp0[3]);

            return result;
        }

    }
}
