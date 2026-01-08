using System;
using System.IO;
using OpenTK.Graphics.OpenGL;
using OpenTK.Mathematics;

namespace KimeraCS.Rendering
{
    public class ShaderProgram : IDisposable
    {
        public int Handle { get; private set; }
        private bool _disposed;

        public ShaderProgram(string vertexSource, string fragmentSource)
        {
            int vertexShader = CompileShader(ShaderType.VertexShader, vertexSource);
            int fragmentShader = CompileShader(ShaderType.FragmentShader, fragmentSource);

            Handle = GL.CreateProgram();
            GL.AttachShader(Handle, vertexShader);
            GL.AttachShader(Handle, fragmentShader);
            GL.LinkProgram(Handle);

            int success;
            GL.GetProgram(Handle, GetProgramParameterName.LinkStatus, out success);
            if (success == 0)
            {
                string infoLog = GL.GetProgramInfoLog(Handle);
                throw new Exception($"Shader program linking failed: {infoLog}");
            }

            GL.DetachShader(Handle, vertexShader);
            GL.DetachShader(Handle, fragmentShader);
            GL.DeleteShader(vertexShader);
            GL.DeleteShader(fragmentShader);
        }

        public static ShaderProgram FromFiles(string vertexPath, string fragmentPath)
        {
            string vertexSource = File.ReadAllText(vertexPath);
            string fragmentSource = File.ReadAllText(fragmentPath);
            return new ShaderProgram(vertexSource, fragmentSource);
        }

        private static int CompileShader(ShaderType type, string source)
        {
            int shader = GL.CreateShader(type);
            GL.ShaderSource(shader, source);
            GL.CompileShader(shader);

            int success;
            GL.GetShader(shader, ShaderParameter.CompileStatus, out success);
            if (success == 0)
            {
                string infoLog = GL.GetShaderInfoLog(shader);
                throw new Exception($"{type} compilation failed: {infoLog}");
            }

            return shader;
        }

        public void Use()
        {
            GL.UseProgram(Handle);
        }

        public int GetUniformLocation(string name)
        {
            return GL.GetUniformLocation(Handle, name);
        }

        public void SetBool(string name, bool value)
        {
            GL.Uniform1(GetUniformLocation(name), value ? 1 : 0);
        }

        public void SetInt(string name, int value)
        {
            GL.Uniform1(GetUniformLocation(name), value);
        }

        public void SetFloat(string name, float value)
        {
            GL.Uniform1(GetUniformLocation(name), value);
        }

        public void SetVector2(string name, Vector2 value)
        {
            GL.Uniform2(GetUniformLocation(name), value.X, value.Y);
        }

        public void SetVector3(string name, Vector3 value)
        {
            GL.Uniform3(GetUniformLocation(name), value.X, value.Y, value.Z);
        }

        public void SetVector4(string name, Vector4 value)
        {
            GL.Uniform4(GetUniformLocation(name), value.X, value.Y, value.Z, value.W);
        }

        public void SetMatrix4(string name, Matrix4 value)
        {
            GL.UniformMatrix4(GetUniformLocation(name), false, ref value);
        }

        public void SetMatrix4(string name, ref Matrix4 value)
        {
            GL.UniformMatrix4(GetUniformLocation(name), false, ref value);
        }

        protected virtual void Dispose(bool disposing)
        {
            if (!_disposed)
            {
                if (Handle != 0)
                {
                    GL.DeleteProgram(Handle);
                    Handle = 0;
                }
                _disposed = true;
            }
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        ~ShaderProgram()
        {
            Dispose(false);
        }
    }
}
