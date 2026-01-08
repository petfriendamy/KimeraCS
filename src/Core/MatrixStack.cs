using System;
using System.Collections.Generic;
using OpenTK.Mathematics;

namespace KimeraCS.Rendering
{
    /// <summary>
    /// CPU-side matrix stack to replace legacy OpenGL matrix stack operations.
    /// Provides equivalent functionality to GL.PushMatrix, GL.PopMatrix, GL.LoadIdentity,
    /// GL.Translate, GL.Rotate, GL.Scale, and GL.MultMatrix.
    /// </summary>
    public class MatrixStack
    {
        private readonly Stack<Matrix4> _stack = new Stack<Matrix4>();

        /// <summary>
        /// The current matrix at the top of the stack.
        /// </summary>
        public Matrix4 Current { get; private set; } = Matrix4.Identity;

        /// <summary>
        /// Push the current matrix onto the stack (like GL.PushMatrix).
        /// </summary>
        public void Push()
        {
            _stack.Push(Current);
        }

        /// <summary>
        /// Pop a matrix from the stack and make it current (like GL.PopMatrix).
        /// </summary>
        public void Pop()
        {
            if (_stack.Count > 0)
            {
                Current = _stack.Pop();
            }
            else
            {
                Current = Matrix4.Identity;
            }
        }

        /// <summary>
        /// Reset the current matrix to identity (like GL.LoadIdentity).
        /// </summary>
        public void LoadIdentity()
        {
            Current = Matrix4.Identity;
        }

        /// <summary>
        /// Load a specific matrix as current.
        /// </summary>
        public void LoadMatrix(Matrix4 matrix)
        {
            Current = matrix;
        }

        /// <summary>
        /// Multiply the current matrix by another matrix (like GL.MultMatrix).
        /// </summary>
        public void MultMatrix(Matrix4 matrix)
        {
            Current = Current * matrix;
        }

        /// <summary>
        /// Multiply the current matrix by a double[16] matrix in OpenGL column-major order.
        /// </summary>
        public void MultMatrix(double[] matrix)
        {
            if (matrix == null || matrix.Length < 16)
                return;

            Matrix4 m = new Matrix4(
                (float)matrix[0], (float)matrix[1], (float)matrix[2], (float)matrix[3],
                (float)matrix[4], (float)matrix[5], (float)matrix[6], (float)matrix[7],
                (float)matrix[8], (float)matrix[9], (float)matrix[10], (float)matrix[11],
                (float)matrix[12], (float)matrix[13], (float)matrix[14], (float)matrix[15]
            );
            Current = Current * m;
        }

        /// <summary>
        /// Apply a translation to the current matrix (like GL.Translate).
        /// </summary>
        public void Translate(float x, float y, float z)
        {
            Current = Current * Matrix4.CreateTranslation(x, y, z);
        }

        /// <summary>
        /// Apply a translation to the current matrix (double overload).
        /// </summary>
        public void Translate(double x, double y, double z)
        {
            Translate((float)x, (float)y, (float)z);
        }

        /// <summary>
        /// Apply a rotation around an axis to the current matrix (like GL.Rotate).
        /// Angle is in degrees.
        /// </summary>
        public void Rotate(float angle, float x, float y, float z)
        {
            float radians = MathHelper.DegreesToRadians(angle);
            Vector3 axis = new Vector3(x, y, z);
            if (axis.LengthSquared > 0)
            {
                axis.Normalize();
                Current = Current * Matrix4.CreateFromAxisAngle(axis, radians);
            }
        }

        /// <summary>
        /// Apply a rotation around an axis (double overload).
        /// </summary>
        public void Rotate(double angle, double x, double y, double z)
        {
            Rotate((float)angle, (float)x, (float)y, (float)z);
        }

        /// <summary>
        /// Apply scaling to the current matrix (like GL.Scale).
        /// </summary>
        public void Scale(float x, float y, float z)
        {
            Current = Current * Matrix4.CreateScale(x, y, z);
        }

        /// <summary>
        /// Apply scaling (double overload).
        /// </summary>
        public void Scale(double x, double y, double z)
        {
            Scale((float)x, (float)y, (float)z);
        }

        /// <summary>
        /// Clear the stack and reset to identity.
        /// </summary>
        public void Clear()
        {
            _stack.Clear();
            Current = Matrix4.Identity;
        }

        /// <summary>
        /// Get the number of matrices on the stack (not including current).
        /// </summary>
        public int Count => _stack.Count;

        /// <summary>
        /// Pop all matrices from the stack until a specific depth is reached.
        /// </summary>
        public void PopToDepth(int depth)
        {
            while (_stack.Count > depth && _stack.Count > 0)
            {
                Pop();
            }
        }
    }

    /// <summary>
    /// Manages separate projection and modelview matrix stacks,
    /// similar to the legacy OpenGL GL.MatrixMode system.
    /// </summary>
    public static class MatrixManager
    {
        /// <summary>
        /// The projection matrix stack.
        /// </summary>
        public static MatrixStack Projection { get; } = new MatrixStack();

        /// <summary>
        /// The modelview matrix stack.
        /// </summary>
        public static MatrixStack ModelView { get; } = new MatrixStack();

        /// <summary>
        /// The currently active matrix stack (for compatibility with GL.MatrixMode pattern).
        /// </summary>
        public static MatrixStack ActiveStack { get; private set; } = ModelView;

        /// <summary>
        /// Set the active matrix stack (like GL.MatrixMode).
        /// </summary>
        public static void SetMatrixMode(CpuMatrixMode mode)
        {
            ActiveStack = mode == CpuMatrixMode.Projection ? Projection : ModelView;
        }

        /// <summary>
        /// Push the active stack.
        /// </summary>
        public static void PushMatrix() => ActiveStack.Push();

        /// <summary>
        /// Pop the active stack.
        /// </summary>
        public static void PopMatrix() => ActiveStack.Pop();

        /// <summary>
        /// Load identity on the active stack.
        /// </summary>
        public static void LoadIdentity() => ActiveStack.LoadIdentity();

        /// <summary>
        /// Translate the active stack.
        /// </summary>
        public static void Translate(float x, float y, float z) => ActiveStack.Translate(x, y, z);
        public static void Translate(double x, double y, double z) => ActiveStack.Translate(x, y, z);

        /// <summary>
        /// Rotate the active stack.
        /// </summary>
        public static void Rotate(float angle, float x, float y, float z) => ActiveStack.Rotate(angle, x, y, z);
        public static void Rotate(double angle, double x, double y, double z) => ActiveStack.Rotate(angle, x, y, z);

        /// <summary>
        /// Scale the active stack.
        /// </summary>
        public static void Scale(float x, float y, float z) => ActiveStack.Scale(x, y, z);
        public static void Scale(double x, double y, double z) => ActiveStack.Scale(x, y, z);

        /// <summary>
        /// Multiply the active stack by a matrix.
        /// </summary>
        public static void MultMatrix(Matrix4 matrix) => ActiveStack.MultMatrix(matrix);
        public static void MultMatrix(double[] matrix) => ActiveStack.MultMatrix(matrix);

        /// <summary>
        /// Sync the current matrices to GLRenderer for shader use.
        /// </summary>
        public static void SyncToRenderer()
        {
            GLRenderer.ProjectionMatrix = Projection.Current;
            GLRenderer.ModelMatrix = ModelView.Current;
        }

        /// <summary>
        /// Reset both stacks to identity.
        /// </summary>
        public static void Reset()
        {
            Projection.Clear();
            ModelView.Clear();
        }
    }

    /// <summary>
    /// CPU-side matrix mode enum for MatrixManager (not related to GL.MatrixMode).
    /// </summary>
    public enum CpuMatrixMode
    {
        Projection,
        Modelview
    }
}
