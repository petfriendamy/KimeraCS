using System;
using System.Collections.Generic;
using OpenTK.Mathematics;

namespace KimeraCS.Rendering
{
    using static FF7PModel;
    using static Utils;

    /// <summary>
    /// Helper functions to create mesh objects for visualization (axes, bounding boxes, normals, etc.)
    /// These replace legacy immediate-mode OpenGL drawing with modern VAO/VBO meshes.
    /// </summary>
    public static class VisualizationHelpers
    {
        /// <summary>
        /// Create a LineMesh for displaying vertex or face normals.
        /// </summary>
        public static LineMesh CreateNormalsMesh(PGroup group, PPolygon[] polys, Point3D[] verts,
                                                   Point3D[] normals, int[] normalsIndex,
                                                   bool showVertexNormals, bool showFaceNormals,
                                                   float normalsScale, int normalsColor)
        {
            if (group.HiddenQ || normals == null || normals.Length == 0)
                return null;

            var vertices = new List<LineVertex>();

            float r = (normalsColor & 0x1) == 0x1 ? 1.0f : 0.0f;
            float g = (normalsColor & 0x2) == 0x2 ? 1.0f : 0.0f;
            float b = (normalsColor & 0x4) == 0x4 ? 1.0f : 0.0f;
            var color = new Vector4(r, g, b, 1.0f);

            for (int iPolyIdx = group.offsetPoly; iPolyIdx < group.offsetPoly + group.numPoly; iPolyIdx++)
            {
                if (showVertexNormals)
                {
                    for (int iVertIdx = 0; iVertIdx < 3; iVertIdx++)
                    {
                        int vertIdx = polys[iPolyIdx].Verts[iVertIdx] + group.offsetVert;
                        float x = verts[vertIdx].x;
                        float y = verts[vertIdx].y;
                        float z = verts[vertIdx].z;

                        int normIdx = normalsIndex[vertIdx];
                        if (normIdx >= 0 && normIdx < normals.Length)
                        {
                            float xn = x + normals[normIdx].x * normalsScale;
                            float yn = y + normals[normIdx].y * normalsScale;
                            float zn = z + normals[normIdx].z * normalsScale;

                            vertices.Add(new LineVertex { Position = new Vector3(x, y, z), Color = color });
                            vertices.Add(new LineVertex { Position = new Vector3(xn, yn, zn), Color = color });
                        }
                    }
                }
                else if (showFaceNormals)
                {
                    Point3D v0 = verts[polys[iPolyIdx].Verts[0] + group.offsetVert];
                    Point3D v1 = verts[polys[iPolyIdx].Verts[1] + group.offsetVert];
                    Point3D v2 = verts[polys[iPolyIdx].Verts[2] + group.offsetVert];

                    Point3D centroid = CalculateCenteroid(v0, v1, v2);
                    Point3D normal = CalculateNormal(v0, v1, v2);
                    normal = Normalize(normal);

                    vertices.Add(new LineVertex
                    {
                        Position = new Vector3(centroid.x, centroid.y, centroid.z),
                        Color = color
                    });
                    vertices.Add(new LineVertex
                    {
                        Position = new Vector3(
                            centroid.x + (-normal.x * normalsScale),
                            centroid.y + (-normal.y * normalsScale),
                            centroid.z + (-normal.z * normalsScale)),
                        Color = color
                    });
                }
            }

            if (vertices.Count == 0)
                return null;

            var mesh = new LineMesh();
            mesh.Upload(vertices.ToArray());
            return mesh;
        }

        /// <summary>
        /// Create a LineMesh for a wireframe bounding box.
        /// </summary>
        public static LineMesh CreateBoundingBoxMesh(float maxX, float maxY, float maxZ,
                                                       float minX, float minY, float minZ,
                                                       float r, float g, float b)
        {
            var color = new Vector4(r, g, b, 1.0f);
            var vertices = new LineVertex[]
            {
                // From max corner
                new LineVertex { Position = new Vector3(maxX, maxY, maxZ), Color = color },
                new LineVertex { Position = new Vector3(maxX, maxY, minZ), Color = color },
                new LineVertex { Position = new Vector3(maxX, maxY, maxZ), Color = color },
                new LineVertex { Position = new Vector3(maxX, minY, maxZ), Color = color },
                new LineVertex { Position = new Vector3(maxX, maxY, maxZ), Color = color },
                new LineVertex { Position = new Vector3(minX, maxY, maxZ), Color = color },

                // From min corner
                new LineVertex { Position = new Vector3(minX, minY, minZ), Color = color },
                new LineVertex { Position = new Vector3(minX, minY, maxZ), Color = color },
                new LineVertex { Position = new Vector3(minX, minY, minZ), Color = color },
                new LineVertex { Position = new Vector3(minX, maxY, minZ), Color = color },
                new LineVertex { Position = new Vector3(minX, minY, minZ), Color = color },
                new LineVertex { Position = new Vector3(maxX, minY, minZ), Color = color },

                // Remaining edges
                new LineVertex { Position = new Vector3(maxX, minY, minZ), Color = color },
                new LineVertex { Position = new Vector3(maxX, maxY, minZ), Color = color },
                new LineVertex { Position = new Vector3(maxX, minY, minZ), Color = color },
                new LineVertex { Position = new Vector3(maxX, minY, maxZ), Color = color },

                new LineVertex { Position = new Vector3(minX, maxY, minZ), Color = color },
                new LineVertex { Position = new Vector3(minX, maxY, maxZ), Color = color },
                new LineVertex { Position = new Vector3(minX, maxY, minZ), Color = color },
                new LineVertex { Position = new Vector3(maxX, maxY, minZ), Color = color },

                new LineVertex { Position = new Vector3(minX, minY, maxZ), Color = color },
                new LineVertex { Position = new Vector3(minX, maxY, maxZ), Color = color },
                new LineVertex { Position = new Vector3(minX, minY, maxZ), Color = color },
                new LineVertex { Position = new Vector3(maxX, minY, maxZ), Color = color },
            };

            var mesh = new LineMesh();
            mesh.Upload(vertices);
            return mesh;
        }

        /// <summary>
        /// Create a LineMesh for a wireframe bounding box from PBoundingBox struct.
        /// </summary>
        public static LineMesh CreateBoundingBoxMesh(PBoundingBox bbox, float r, float g, float b)
        {
            return CreateBoundingBoxMesh(bbox.max_x, bbox.max_y, bbox.max_z,
                                         bbox.min_x, bbox.min_y, bbox.min_z,
                                         r, g, b);
        }

        /// <summary>
        /// Create a LineMesh for 3D coordinate axes.
        /// </summary>
        public static LineMesh CreateAxesMesh(float axisLength, bool isBattleLocation)
        {
            var vertices = new List<LineVertex>();
            var red = new Vector4(1, 0, 0, 1);
            var green = new Vector4(0, 1, 0, 1);
            var blue = new Vector4(0, 0, 1, 1);

            // X axis (always red)
            vertices.Add(new LineVertex { Position = Vector3.Zero, Color = red });
            vertices.Add(new LineVertex { Position = new Vector3(axisLength, 0, 0), Color = red });

            if (isBattleLocation)
            {
                // Y axis (green in battle)
                vertices.Add(new LineVertex { Position = Vector3.Zero, Color = green });
                vertices.Add(new LineVertex { Position = new Vector3(0, -axisLength, 0), Color = green });

                // Z axis (blue in battle)
                vertices.Add(new LineVertex { Position = Vector3.Zero, Color = blue });
                vertices.Add(new LineVertex { Position = new Vector3(0, 0, axisLength), Color = blue });
            }
            else
            {
                // Y axis (blue in field - swapped)
                vertices.Add(new LineVertex { Position = Vector3.Zero, Color = blue });
                vertices.Add(new LineVertex { Position = new Vector3(0, -axisLength, 0), Color = blue });

                // Z axis (green in field - swapped)
                vertices.Add(new LineVertex { Position = Vector3.Zero, Color = green });
                vertices.Add(new LineVertex { Position = new Vector3(0, 0, axisLength), Color = green });
            }

            var mesh = new LineMesh();
            mesh.Upload(vertices.ToArray());
            return mesh;
        }

        /// <summary>
        /// Create a LineMesh for P editor coordinate axes (simpler version).
        /// </summary>
        public static LineMesh CreateAxesMeshPE(float axisLength)
        {
            var red = new Vector4(1, 0, 0, 1);
            var green = new Vector4(0, 1, 0, 1);
            var blue = new Vector4(0, 0, 1, 1);

            var vertices = new LineVertex[]
            {
                // X axis
                new LineVertex { Position = Vector3.Zero, Color = red },
                new LineVertex { Position = new Vector3(axisLength, 0, 0), Color = red },
                // Y axis
                new LineVertex { Position = Vector3.Zero, Color = green },
                new LineVertex { Position = new Vector3(0, axisLength, 0), Color = green },
                // Z axis
                new LineVertex { Position = Vector3.Zero, Color = blue },
                new LineVertex { Position = new Vector3(0, 0, axisLength), Color = blue },
            };

            var mesh = new LineMesh();
            mesh.Upload(vertices);
            return mesh;
        }

        /// <summary>
        /// Create a LineMesh for the ground grid lines.
        /// </summary>
        public static LineMesh CreateGroundGridMesh()
        {
            var vertices = new List<LineVertex>();
            var color = new Vector4(0.0f, 0.0f, 0.0f, 1.0f);

            // Main cross
            vertices.Add(new LineVertex { Position = new Vector3(0, 0, 50), Color = color });
            vertices.Add(new LineVertex { Position = new Vector3(0, 0, -50), Color = color });
            vertices.Add(new LineVertex { Position = new Vector3(-50, 0, 0), Color = color });
            vertices.Add(new LineVertex { Position = new Vector3(50, 0, 0), Color = color });

            var mesh = new LineMesh();
            mesh.Upload(vertices.ToArray());
            return mesh;
        }

        /// <summary>
        /// Create a mesh for the ground plane (as two triangles).
        /// Returns vertices for a simple quad.
        /// </summary>
        public static GroupMesh CreateGroundPlaneMesh()
        {
            var vertices = new Vertex[]
            {
                new Vertex { Position = new Vector3(300, 0, 300), Normal = Vector3.UnitY, Color = new Vector4(0.9f, 0.9f, 1f, 1f) },
                new Vertex { Position = new Vector3(300, 0, -300), Normal = Vector3.UnitY, Color = new Vector4(0.9f, 0.9f, 1f, 1f) },
                new Vertex { Position = new Vector3(-300, 0, -300), Normal = Vector3.UnitY, Color = new Vector4(0.9f, 0.9f, 1f, 1f) },
                new Vertex { Position = new Vector3(-300, 0, 300), Normal = Vector3.UnitY, Color = new Vector4(0.9f, 0.9f, 1f, 1f) },
            };

            var indices = new uint[] { 0, 1, 2, 0, 2, 3 };

            var mesh = new GroupMesh();
            mesh.Upload(vertices, indices);
            return mesh;
        }

        /// <summary>
        /// Create a mesh for the shadow (triangle fan as triangles).
        /// </summary>
        public static GroupMesh CreateShadowMesh(Point3D pMin, Point3D pMax, int numSegments = 20)
        {
            float cx = (pMin.x + pMax.x) / 2;
            float cz = (pMin.z + pMax.z) / 2;

            Point3D pMinAux = pMin;
            Point3D pMaxAux = pMax;
            pMinAux.y = 0;
            pMaxAux.y = 0;
            float radius = CalculateDistance(pMinAux, pMaxAux) / 2;

            var vertices = new List<Vertex>();
            var indices = new List<uint>();

            // Center vertex (dark, semi-transparent)
            vertices.Add(new Vertex
            {
                Position = new Vector3(cx, 0, cz),
                Normal = Vector3.UnitY,
                Color = new Vector4(0.1f, 0.1f, 0.1f, 0.5f)
            });

            // Edge vertices (dark, transparent)
            for (int i = 0; i <= numSegments; i++)
            {
                float angle = i * 2 * (float)Math.PI / numSegments;
                vertices.Add(new Vertex
                {
                    Position = new Vector3(
                        radius * (float)Math.Sin(angle) + cx,
                        0,
                        radius * (float)Math.Cos(angle) + cz),
                    Normal = Vector3.UnitY,
                    Color = new Vector4(0.1f, 0.1f, 0.1f, 0.0f)
                });
            }

            // Create triangle fan as individual triangles
            for (int i = 1; i <= numSegments; i++)
            {
                indices.Add(0);  // Center
                indices.Add((uint)i);
                indices.Add((uint)(i + 1));
            }

            var mesh = new GroupMesh();
            mesh.Upload(vertices.ToArray(), indices.ToArray());
            return mesh;
        }

        /// <summary>
        /// Create a LineMesh for wireframe rendering of a PModel.
        /// </summary>
        public static LineMesh CreateWireframeMesh(PModel model)
        {
            var vertices = new List<LineVertex>();
            var color = new Vector4(0, 0, 0, 1);

            for (int iGroupIdx = 0; iGroupIdx < model.Header.numGroups; iGroupIdx++)
            {
                if (model.Groups[iGroupIdx].HiddenQ)
                    continue;

                for (int iPolyIdx = model.Groups[iGroupIdx].offsetPoly;
                     iPolyIdx < model.Groups[iGroupIdx].offsetPoly + model.Groups[iGroupIdx].numPoly;
                     iPolyIdx++)
                {
                    int v0 = model.Polys[iPolyIdx].Verts[0] + model.Groups[iGroupIdx].offsetVert;
                    int v1 = model.Polys[iPolyIdx].Verts[1] + model.Groups[iGroupIdx].offsetVert;
                    int v2 = model.Polys[iPolyIdx].Verts[2] + model.Groups[iGroupIdx].offsetVert;

                    Vector3 p0 = new Vector3(model.Verts[v0].x, model.Verts[v0].y, model.Verts[v0].z);
                    Vector3 p1 = new Vector3(model.Verts[v1].x, model.Verts[v1].y, model.Verts[v1].z);
                    Vector3 p2 = new Vector3(model.Verts[v2].x, model.Verts[v2].y, model.Verts[v2].z);

                    // Three edges per triangle
                    vertices.Add(new LineVertex { Position = p0, Color = color });
                    vertices.Add(new LineVertex { Position = p1, Color = color });

                    vertices.Add(new LineVertex { Position = p1, Color = color });
                    vertices.Add(new LineVertex { Position = p2, Color = color });

                    vertices.Add(new LineVertex { Position = p2, Color = color });
                    vertices.Add(new LineVertex { Position = p0, Color = color });
                }
            }

            if (vertices.Count == 0)
                return null;

            var mesh = new LineMesh();
            mesh.Upload(vertices.ToArray());
            return mesh;
        }

        /// <summary>
        /// Create a PointMesh for skeleton joint visualization.
        /// </summary>
        public static PointMesh CreateJointPointsMesh(Vector3[] jointPositions, Vector4 color, float pointSize = 5.0f)
        {
            var vertices = new LineVertex[jointPositions.Length];
            for (int i = 0; i < jointPositions.Length; i++)
            {
                vertices[i] = new LineVertex { Position = jointPositions[i], Color = color };
            }

            var mesh = new PointMesh();
            mesh.PointSize = pointSize;
            mesh.Upload(vertices);
            return mesh;
        }

        /// <summary>
        /// Create a LineMesh for skeleton bone connections.
        /// </summary>
        public static LineMesh CreateBoneLinesMesh(Vector3[] startPoints, Vector3[] endPoints, Vector4 color)
        {
            if (startPoints.Length != endPoints.Length)
                throw new ArgumentException("Start and end point arrays must have the same length");

            var vertices = new LineVertex[startPoints.Length * 2];
            for (int i = 0; i < startPoints.Length; i++)
            {
                vertices[i * 2] = new LineVertex { Position = startPoints[i], Color = color };
                vertices[i * 2 + 1] = new LineVertex { Position = endPoints[i], Color = color };
            }

            var mesh = new LineMesh();
            mesh.Upload(vertices);
            return mesh;
        }

        /// <summary>
        /// Create a 2D overlay LineMesh for axis labels (X, Y, Z).
        /// Coordinates should be in screen space (0,0 = bottom-left).
        /// </summary>
        public static LineMesh CreateAxisLabels2D(float xPosX, float xPosY, float yPosX, float yPosY,
                                                    float zPosX, float zPosY, float letterSize,
                                                    bool isBattleLocation)
        {
            var vertices = new List<LineVertex>();
            var black = new Vector4(0, 0, 0, 1);
            float w = letterSize;
            float h = letterSize * 1.5f;

            // Draw X at xPos
            vertices.Add(new LineVertex { Position = new Vector3(xPosX - w, xPosY - h, 0), Color = black });
            vertices.Add(new LineVertex { Position = new Vector3(xPosX + w, xPosY + h, 0), Color = black });
            vertices.Add(new LineVertex { Position = new Vector3(xPosX - w, xPosY + h, 0), Color = black });
            vertices.Add(new LineVertex { Position = new Vector3(xPosX + w, xPosY - h, 0), Color = black });

            if (isBattleLocation)
            {
                // Draw Y at yPos
                vertices.Add(new LineVertex { Position = new Vector3(yPosX - w, yPosY - h, 0), Color = black });
                vertices.Add(new LineVertex { Position = new Vector3(yPosX + w, yPosY + h, 0), Color = black });
                vertices.Add(new LineVertex { Position = new Vector3(yPosX - w, yPosY + h, 0), Color = black });
                vertices.Add(new LineVertex { Position = new Vector3(yPosX, yPosY, 0), Color = black });

                // Draw Z at zPos
                vertices.Add(new LineVertex { Position = new Vector3(zPosX + w, zPosY + h, 0), Color = black });
                vertices.Add(new LineVertex { Position = new Vector3(zPosX - w, zPosY + h, 0), Color = black });
                vertices.Add(new LineVertex { Position = new Vector3(zPosX + w, zPosY + h, 0), Color = black });
                vertices.Add(new LineVertex { Position = new Vector3(zPosX - w, zPosY - h, 0), Color = black });
                vertices.Add(new LineVertex { Position = new Vector3(zPosX + w, zPosY - h, 0), Color = black });
                vertices.Add(new LineVertex { Position = new Vector3(zPosX - w, zPosY - h, 0), Color = black });
            }
            else
            {
                // Draw Y at zPos (swapped)
                vertices.Add(new LineVertex { Position = new Vector3(zPosX - w, zPosY - h, 0), Color = black });
                vertices.Add(new LineVertex { Position = new Vector3(zPosX + w, zPosY + h, 0), Color = black });
                vertices.Add(new LineVertex { Position = new Vector3(zPosX - w, zPosY + h, 0), Color = black });
                vertices.Add(new LineVertex { Position = new Vector3(zPosX, zPosY, 0), Color = black });

                // Draw Z at yPos (swapped)
                vertices.Add(new LineVertex { Position = new Vector3(yPosX + w, yPosY + h, 0), Color = black });
                vertices.Add(new LineVertex { Position = new Vector3(yPosX - w, yPosY + h, 0), Color = black });
                vertices.Add(new LineVertex { Position = new Vector3(yPosX + w, yPosY + h, 0), Color = black });
                vertices.Add(new LineVertex { Position = new Vector3(yPosX - w, yPosY - h, 0), Color = black });
                vertices.Add(new LineVertex { Position = new Vector3(yPosX + w, yPosY - h, 0), Color = black });
                vertices.Add(new LineVertex { Position = new Vector3(yPosX - w, yPosY - h, 0), Color = black });
            }

            var mesh = new LineMesh();
            mesh.Upload(vertices.ToArray());
            return mesh;
        }

        /// <summary>
        /// Create a simple line mesh from two points.
        /// </summary>
        public static LineMesh CreateSimpleLineMesh(Vector3 start, Vector3 end, Vector4 color)
        {
            var vertices = new LineVertex[]
            {
                new LineVertex { Position = start, Color = color },
                new LineVertex { Position = end, Color = color }
            };

            var mesh = new LineMesh();
            mesh.Upload(vertices);
            return mesh;
        }

        /// <summary>
        /// Create a LineMesh for multiple colored lines.
        /// </summary>
        public static LineMesh CreateMultiColorLineMesh(List<(Vector3 Start, Vector3 End, Vector4 Color)> lines)
        {
            var vertices = new LineVertex[lines.Count * 2];
            for (int i = 0; i < lines.Count; i++)
            {
                vertices[i * 2] = new LineVertex { Position = lines[i].Start, Color = lines[i].Color };
                vertices[i * 2 + 1] = new LineVertex { Position = lines[i].End, Color = lines[i].Color };
            }

            var mesh = new LineMesh();
            mesh.Upload(vertices);
            return mesh;
        }
    }
}
