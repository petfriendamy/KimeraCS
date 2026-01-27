using Assimp;
using OpenTK.Mathematics;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using System.Linq;
//using System.Numerics;
using System.Text;


namespace KimeraCS.Core
{
    using static FF7BattleAnimation;
    using static FF7BattleAnimationsPack;
    using static FF7BattleSkeleton;
    using static FF7FieldAnimation;
    using static FF7FieldRSDResource;
    using static FF7FieldSkeleton;
    using static FF7PModel;
    using static FF7TEXTexture;
    using static Utils;

    public static class ModelConverter
    {
        private static AssimpContext _context = new();

        #region IMPORT FUNCTIONS

        public static bool IsValidImport(string path)
        {
            string fileExt = Path.GetExtension(path).ToLower();
            return _context.IsImportFormatSupported(fileExt);
        }

        public static string GetFileFilter()
        {
            var extensions = _context.GetSupportedImportFormats();
            var builder = new StringBuilder("External Models|");
            foreach (var filter in extensions)
            {
                builder.Append('*' + filter + ';');
            }
            builder.Remove(builder.Length - 1, 1); //remove last semicolon
            return builder.ToString();
        }

        /// <summary>
        /// Extracts texture file paths from an Assimp scene's materials.
        /// </summary>
        /// <param name="scene">The Assimp scene to extract textures from</param>
        /// <returns>Dictionary mapping material index to texture file path</returns>
        private static Dictionary<int, string> ExtractTexturePathsFromScene(Scene scene)
        {
            var texturePaths = new Dictionary<int, string>();

            if (!scene.HasMaterials)
                return texturePaths;

            for (int i = 0; i < scene.MaterialCount; i++)
            {
                var material = scene.Materials[i];

                // Check for diffuse texture
                if (material.HasTextureDiffuse)
                {
                    var textureSlot = material.TextureDiffuse;
                    if (!string.IsNullOrEmpty(textureSlot.FilePath))
                    {
                        texturePaths[i] = textureSlot.FilePath;
                    }
                }
            }

            return texturePaths;
        }

        /// <summary>
        /// Loads textures from an Assimp scene, searching in the specified base directory.
        /// </summary>
        /// <param name="scene">The Assimp scene containing texture references</param>
        /// <param name="baseDirectory">Directory to search for texture files (usually the model's directory)</param>
        /// <returns>List of loaded TEX textures, with indices matching material texture references</returns>
        private static List<TEX> LoadTexturesFromScene(Scene scene, string baseDirectory)
        {
            var textures = new List<TEX>();
            var texturePaths = ExtractTexturePathsFromScene(scene);
            var loadedTextures = new List<string>(); // Map path to index to avoid duplicates

            foreach (var kvp in texturePaths)
            {
                string texturePath = kvp.Value;

                // Normalize the path
                string normalizedPath = texturePath.Replace('/', Path.DirectorySeparatorChar)
                                                   .Replace('\\', Path.DirectorySeparatorChar);

                // Try to find the texture file
                string fullPath = "";

                // Try as relative path from base directory
                string relativePath = Path.Combine(baseDirectory, normalizedPath);
                if (File.Exists(relativePath))
                {
                    fullPath = relativePath;
                }
                // Try just the filename in base directory
                else if (File.Exists(Path.Combine(baseDirectory, Path.GetFileName(normalizedPath))))
                {
                    fullPath = Path.Combine(baseDirectory, Path.GetFileName(normalizedPath));
                }
                // Try as absolute path
                else if (File.Exists(normalizedPath))
                {
                    fullPath = normalizedPath;
                }

                if (!string.IsNullOrEmpty(fullPath))
                {
                    // Check if already loaded
                    if (!loadedTextures.Contains(fullPath))
                    {
                        TEX tex = new();
                        LoadImageAsTEXTexture(fullPath, ref tex);
                        textures.Add(tex);
                        loadedTextures.Add(fullPath);
                    }
                }
            }

            return textures;
        }

        public static Scene? LoadSceneFromFile(string filePath)
        {
            if (IsValidImport(filePath))
                try
                {
                    return _context.ImportFile(filePath,
                    PostProcessSteps.Triangulate |
                    PostProcessSteps.GenerateNormals |
                    PostProcessSteps.FlipUVs);
                }
                catch (Exception ex)
                {
                    int temp = 1;
                }
            return null;
        }

        private static void GetVerts(Mesh mesh, out Vector3[] vertsV, bool bAdjust)
        {
            vertsV = new Vector3[mesh.VertexCount];

            for (int i = 0; i < mesh.VertexCount; i++)
            {
                float x = mesh.Vertices[i].X;
                float y = mesh.Vertices[i].Y;
                float z = mesh.Vertices[i].Z;

                if (bAdjust)
                {
                    // Export negates X and Z when bAdjust is true
                    // Import must also negate X and Z to reverse this
                    x = -x;
                    z = -z;
                }

                vertsV[i].X = x;
                vertsV[i].Y = y;
                vertsV[i].Z = z;
            }
        }

        private static void GetFaces(Mesh mesh, out PPolygon[] facesV, bool bAdjust)
        {
            // Count triangulated faces
            int triCount = 0;
            foreach (var face in mesh.Faces)
            {
                if (face.IndexCount >= 3)
                    triCount += face.IndexCount - 2; // Fan triangulation
            }

            facesV = new PPolygon[triCount];
            int fi = 0;

            foreach (var face in mesh.Faces)
            {
                if (face.IndexCount < 3) continue;

                // Triangulate using fan method (works for convex polygons)
                for (int t = 0; t < face.IndexCount - 2; t++)
                {
                    facesV[fi].tag1 = 0;
                    facesV[fi].tag2 = PPOLY_TAG2;

                    facesV[fi].Verts = new ushort[3];
                    // Export reverses winding order {0,2,1} only when bAdjust is true
                    // So import should reverse only when bAdjust is true to match
                    if (bAdjust)
                    {
                        // Reverse winding order (0, 2, 1) to undo export's reversal
                        facesV[fi].Verts[0] = (ushort)face.Indices[0];
                        facesV[fi].Verts[1] = (ushort)face.Indices[t + 2];
                        facesV[fi].Verts[2] = (ushort)face.Indices[t + 1];
                    }
                    else
                    {
                        // Keep original winding order (0, 1, 2)
                        facesV[fi].Verts[0] = (ushort)face.Indices[0];
                        facesV[fi].Verts[1] = (ushort)face.Indices[t + 1];
                        facesV[fi].Verts[2] = (ushort)face.Indices[t + 2];
                    }

                    facesV[fi].Normals = new ushort[3];
                    facesV[fi].Edges = new ushort[3];
                    fi++;
                }
            }
        }

        private static void GetTexCoords(Mesh mesh, out Vector2[] texCoordsV, bool processCoords)
        {
            texCoordsV = Array.Empty<Vector2>();

            // Check if mesh has texture coordinates (channel 0)
            if (processCoords && mesh.HasTextureCoords(0))
            {
                texCoordsV = new Vector2[mesh.VertexCount];

                for (int i = 0; i < mesh.VertexCount; i++)
                {
                    texCoordsV[i].X = mesh.TextureCoordinateChannels[0][i].X;
                    texCoordsV[i].Y = mesh.TextureCoordinateChannels[0][i].Y;
                }
            }
        }

        private static void GetVColors(Mesh mesh, Material material, out Color[] vcolorsV)
        {
            vcolorsV = new Color[mesh.VertexCount];

            // Check if mesh has vertex colors (channel 0)
            if (mesh.HasVertexColors(0))
            {
                for (int i = 0; i < mesh.VertexCount; i++)
                {
                    var c = mesh.VertexColorChannels[0][i];
                    // Color.FromArgb expects (alpha, red, green, blue)
                    // Assimp Color4D indexer returns [R, G, B, A]
                    vcolorsV[i] = Color.FromArgb(
                        (int)(c[3] * 255),  // A
                        (int)(c[0] * 255),  // R
                        (int)(c[1] * 255),  // G
                        (int)(c[2] * 255)); // B
                }
            }
            else if (material != null && material.HasColorDiffuse)
            {
                // Use material diffuse color for all vertices
                var diffuse = material.ColorDiffuse;
                var color = Color.FromArgb(
                    255,
                    (int)(diffuse[0] * 255),
                    (int)(diffuse[1] * 255),
                    (int)(diffuse[2] * 255));

                for (int i = 0; i < mesh.VertexCount; i++)
                {
                    vcolorsV[i] = color;
                }
            }
            else
            {
                // Default to white
                for (int i = 0; i < mesh.VertexCount; i++)
                {
                    vcolorsV[i] = Color.FromArgb(255, 255, 255, 255);
                }
            }
        }

        private static void GetPColors(Mesh mesh, Material material, int faceCount, out Color[] pcolorsV)
        {
            pcolorsV = new Color[faceCount];

            // If mesh has vertex colors, sample from first vertex of each face
            // This handles per-face colors that were exported as per-vertex colors
            if (mesh.HasVertexColors(0))
            {
                int fi = 0;
                foreach (var face in mesh.Faces)
                {
                    if (face.IndexCount < 3) continue;

                    // For triangulated faces, sample from first vertex
                    for (int t = 0; t < face.IndexCount - 2 && fi < faceCount; t++)
                    {
                        int firstVertIdx = face.Indices[0];
                        if (firstVertIdx < mesh.VertexColorChannels[0].Count)
                        {
                            var c = mesh.VertexColorChannels[0][firstVertIdx];
                            // Color.FromArgb expects (alpha, red, green, blue)
                            // Assimp Color4D indexer returns [R, G, B, A]
                            pcolorsV[fi] = Color.FromArgb(
                                (int)(c[3] * 255),  // A
                                (int)(c[0] * 255),  // R
                                (int)(c[1] * 255),  // G
                                (int)(c[2] * 255)); // B
                        }
                        else
                        {
                            pcolorsV[fi] = Color.FromArgb(255, 255, 255, 255);
                        }
                        fi++;
                    }
                }
                return;
            }

            // Fall back to material diffuse color
            Color polyColor;
            if (material != null && material.HasColorDiffuse)
            {
                var diffuse = material.ColorDiffuse;
                polyColor = Color.FromArgb(
                    255,
                    (int)(diffuse[0] * 255),
                    (int)(diffuse[1] * 255),
                    (int)(diffuse[2] * 255));
            }
            else
            {
                polyColor = Color.FromArgb(255, 255, 255, 255);
            }

            for (int i = 0; i < faceCount; i++)
            {
                pcolorsV[i] = polyColor;
            }
        }
        
        /// <summary>
         /// Calculates bone length for a node, using metadata if available or inferring from children.
         /// </summary>
         /// <param name="node">The bone node to calculate length for</param>
         /// <param name="bAdjust">Whether coordinate adjustments are applied</param>
         /// <returns>The bone length (positive value)</returns>
        private static double CalculateBoneLengthFromNode(Node node, bool bAdjust)
        {
            // 1. First try to get FF7_BoneLength metadata (Kimera exports)
            /*if (node.Metadata.TryGetValue("FF7_BoneLength", out var entry))
            {
                if (entry.Data is double d)
                    return d;
                if (entry.Data is float f)
                    return f;
            }*/

            // 2. Calculate from first child bone's translation
            // The child bone's transform contains this bone's length (as Z translation)
            foreach (var child in node.Children)
            {
                // Skip model nodes - they're not bones
                if (child.Name.StartsWith("Model_"))
                    continue;

                var childTransform = ToMatrix4(child.Transform);
                float z = childTransform.M43;

                // The child's Z translation represents the parent's bone length
                // Account for coordinate system direction
                if (bAdjust)
                {
                    // When bAdjust, positive Z was used in export
                    return z;
                }
                else
                {
                    // When !bAdjust, negative Z was used in export
                    return -z;
                }
            }

            // 3. No children (leaf bone) - default to 0
            return 0;
        }

        /// <summary>
        /// Calculates bone length for a battle skeleton node, using metadata if available or inferring from children.
        /// Battle skeletons use opposite Z direction from field skeletons.
        /// </summary>
        /// <param name="node">The bone node to calculate length for</param>
        /// <param name="bAdjust">Whether coordinate adjustments are applied</param>
        /// <returns>The bone length (positive value)</returns>
        private static float CalculateBattleBoneLengthFromNode(Node node, bool bAdjust)
        {
            // 1. First try to get FF7_BoneLength metadata (Kimera exports)
            if (node.Metadata.TryGetValue("FF7_BoneLength", out var entry))
            {
                if (entry.Data is double d)
                    return (float)d;
                if (entry.Data is float f)
                    return f;
            }

            // 2. Calculate from first child bone's translation
            foreach (var child in node.Children)
            {
                // Skip model nodes - they're not bones
                if (child.Name.StartsWith("Model_"))
                    continue;

                var childTransform = ToMatrix4(child.Transform);
                float z = childTransform.M43;

                // Battle skeletons use opposite Z direction from field skeletons
                if (bAdjust)
                {
                    // When bAdjust, negative Z was used in export
                    return -z;
                }
                else
                {
                    // When !bAdjust, positive Z was used in export
                    return z;
                }
            }

            // 3. No children (leaf bone) - default to 0
            return 0;
        }

        /// <summary>
        /// Flattens an Assimp node hierarchy into an ordered list for import.
        /// </summary>
        private static List<(Node node, Node? parent, int depth)> FlattenHierarchy(Node root)
        {
            var result = new List<(Node, Node?, int)>();
            FlattenHierarchyRecursive(root, null, 0, result);
            return result;
        }

        private static void FlattenHierarchyRecursive(Node node, Node? parent, int depth,
            List<(Node, Node?, int)> result)
        {
            result.Add((node, parent, depth));
            foreach (var child in node.Children)
            {
                FlattenHierarchyRecursive(child, node, depth + 1, result);
            }
        }

        private static void AddMeshAsGroup(Mesh mesh, Material material, ref PModel Model, bool bAdjust, int textureID = -1)
        {
            // bAdjust is passed to GetVerts and GetFaces to handle coordinate system conversion
            // Export negates X and Z vertices and reverses winding when bAdjust=true
            // Import does the same to reverse the transformation
            GetVerts(mesh, out Vector3[] vertsV, bAdjust);
            GetFaces(mesh, out PPolygon[] facesV, bAdjust);
            GetTexCoords(mesh, out Vector2[] texcoordsV, (textureID >= 0));
            GetVColors(mesh, material, out Color[] vcolorsV);
            GetPColors(mesh, material, facesV.Length, out Color[] pcolorsV);

            AddGroup(ref Model, vertsV, facesV, texcoordsV, vcolorsV, pcolorsV, Math.Max(textureID, 0));
        }

        /// <summary>
        /// Converts an Assimp Scene to a PModel.
        /// Supports any format that Assimp can load (3DS, OBJ, FBX, GLTF, etc.)
        /// </summary>
        /// <param name="scene">The Assimp scene to convert</param>
        /// <param name="outModel">The output PModel (should be initialized)</param>
        /// <param name="bAdjustModel">Whether to adjust coordinates for FF7 (mirror Z, rotate 180° Y)</param>
        public static void ConvertSceneToPModel(Scene scene, ref PModel outModel, bool bAdjustModel)
        {
            // Ensure model arrays are initialized (not null) even if scene is empty
            if (outModel.Polys == null) outModel.Polys = Array.Empty<PPolygon>();
            if (outModel.Verts == null) outModel.Verts = Array.Empty<Vector3>();
            if (outModel.Normals == null) outModel.Normals = Array.Empty<Vector3>();

            if (!scene.HasMeshes) return;
            foreach (var mesh in scene.Meshes)
            {
                // Validate mesh limits
                if (mesh.VertexCount > 0xFFFF)
                {
                    throw new TooManyVerticesException(
                        $"Mesh '{mesh.Name}' has {mesh.VertexCount} vertices. " +
                        "The max number of vertices allowed for a FF7 .P model is 65535.");
                }

                if (mesh.FaceCount > 0xFFFF)
                {
                    throw new TooManyFacesException(
                        $"Mesh '{mesh.Name}' has {mesh.FaceCount} faces. " +
                        "The max number of faces allowed for a FF7 .P model is 65535.");
                }

                // Get the material for this mesh
                Material material = new() { Name = "Default" };
                if (scene.HasMaterials && mesh.MaterialIndex >= 0 && mesh.MaterialIndex < scene.MaterialCount)
                {
                    material = scene.Materials[mesh.MaterialIndex];
                }

                AddMeshAsGroup(mesh, material, ref outModel, bAdjustModel);
            }

            // Initialize PModel header and properties
            outModel.Header.version = 1;
            outModel.Header.off04 = 1;
            outModel.Header.unknown = new int[16];

            outModel.resizeX = 1;
            outModel.resizeY = 1;
            outModel.resizeZ = 1;
            outModel.repositionX = 0;
            outModel.repositionY = 0;
            outModel.repositionZ = 0;
            outModel.rotateAlpha = 0;
            outModel.rotateBeta = 0;
            outModel.rotateGamma = 0;
            outModel.rotationQuaternion.X = 0;
            outModel.rotationQuaternion.Y = 0;
            outModel.rotationQuaternion.Z = 0;
            outModel.rotationQuaternion.W = 1;
        }

        //
        // FIELD SKELETON CONVERSION - IMPORT
        //

        /// <summary>
        /// Converts an Assimp Scene to a FieldSkeleton.
        /// </summary>
        /// <param name="scene">The Assimp Scene to convert</param>
        /// <param name="filePath">File path for the file being loaded</param>
        /// <param name="bAdjust">Whether to apply FF7 coordinate adjustments</param>
        /// <returns>A FieldSkeleton constructed from the scene</returns>
        public static FieldSkeleton ConvertSceneToFieldSkeleton(Scene scene, string filePath, bool bAdjust = true)
        {
            //get texture files (if they exist)
            var textures = new List<TEX>();
            var tex_names = new List<string>();
            var directory = Path.GetDirectoryName(filePath);
            if (Directory.Exists(directory))
            {
                textures = LoadTexturesFromScene(scene, directory);
                tex_names =
                    (from t in textures
                     select Path.GetFileNameWithoutExtension(t.TEXfileName)).ToList();
            }

            //get skeleton data
            var skeleton = new FieldSkeleton
            {
                fileName = filePath,
                name = Path.GetFileName(filePath),
                bones = new List<FieldBone>(),
                textures_pool = textures
            };

            if (scene.RootNode == null)
            {
                skeleton.nBones = 0;
                return skeleton;
            }

            // Flatten hierarchy
            var flatNodes = FlattenHierarchy(scene.RootNode);

            // Find the armature node if it exists (created during export as container for bones)
            var armatureNode = scene.RootNode.Children.FirstOrDefault(c => c.Name == "Armature");

            // Skip the root node and armature node itself, process children as bones
            var boneNodes = flatNodes.Where(n => n.node != scene.RootNode &&
                n.node != armatureNode &&
                !n.node.Name.StartsWith("Model_") && // Skip model child nodes
                (n.node.HasMeshes || n.node.ChildCount > 0 || n.node.Metadata.ContainsKey("FF7_BoneIndex"))).ToList();

            foreach (var (node, parent, depth) in boneNodes)
            {
                var bone = new FieldBone
                {
                    fRSDResources = new List<FieldRSDResource>(),
                    nResources = 0
                };

                // Get joint names from metadata or node names
                bone.joint_i = node.Name;//GetFF7MetadataString(node, "FF7_JointI") ?? node.Name;

                // Determine parent joint name - if parent is root/armature, use "root"
                string? parentJointF = null; //GetFF7MetadataString(node, "FF7_JointF");
                if (parentJointF == null)
                {
                    if (parent == null || parent == scene.RootNode || parent == armatureNode)
                        parentJointF = "root";
                    else
                        parentJointF = parent.Name;
                }
                bone.joint_f = parentJointF;

                // Calculate bone length - tries metadata first, then infers from child transforms
                bone.len = CalculateBoneLengthFromNode(node, bAdjust);

                // Get resize values
                bone.resizeX = 1.0f; //GetFF7MetadataFloat(node, "FF7_ResizeX") ?? 1.0f;
                bone.resizeY = 1.0f; // GetFF7MetadataFloat(node, "FF7_ResizeY") ?? 1.0f;
                bone.resizeZ = 1.0f; // GetFF7MetadataFloat(node, "FF7_ResizeZ") ?? 1.0f;

                // Process meshes attached to this node or its model children
                var meshNodes = new List<Node> { node };
                meshNodes.AddRange(node.Children.Where(c => c.HasMeshes && c.Name.StartsWith("Model_")));

                foreach (var meshNode in meshNodes)
                {
                    if (!meshNode.HasMeshes) continue;

                    // Create RSD resource
                    var resource = new FieldRSDResource
                    {
                        ID = "@RSD940102",
                        numTextures = 0,
                        textures = new List<TEX>()
                    };

                    // Create a PModel from the meshes
                    resource.Model = new PModel();
                    InitializePModel(ref resource.Model);

                    foreach (var meshIndex in meshNode.MeshIndices)
                    {
                        if (meshIndex < scene.MeshCount)
                        {
                            var mesh = scene.Meshes[meshIndex];
                            Material material = new() { Name = "Default" };
                            int texId = -1;
                            if (scene.HasMaterials && mesh.MaterialIndex >= 0 && mesh.MaterialIndex < scene.MaterialCount)
                            {
                                material = scene.Materials[mesh.MaterialIndex];
                                foreach (var tex in material.GetAllMaterialTextures())
                                {
                                    if (!string.IsNullOrEmpty(tex.FilePath))
                                    {
                                        string name = Path.GetFileNameWithoutExtension(tex.FilePath);
                                        int index = tex_names.IndexOf(name);
                                        if (index >= 0)
                                        {
                                            resource.textures.Add(textures[index]);
                                            resource.numTextures++;
                                            resource.res_file = tex.FilePath;
                                            texId = index;
                                            break;
                                        }
                                    }
                                }
                            }
                            AddMeshAsGroup(mesh, material, ref resource.Model, bAdjust, texId);
                        }
                    }

                    if (resource.Model.Header.numGroups > 0)
                    {
                        // Extract model transform from the mesh node
                        if (meshNode != node)
                        {
                            var modelTransform = ToMatrix4(meshNode.Transform);
                            resource.Model.repositionX = modelTransform.M41;
                            resource.Model.repositionY = modelTransform.M42;
                            resource.Model.repositionZ = modelTransform.M43;

                            // Extract scale
                            resource.Model.resizeX = modelTransform.Row0.Length;
                            resource.Model.resizeY = modelTransform.Row1.Length;
                            resource.Model.resizeZ = modelTransform.Row2.Length;
                        }

                        string resFile = $"{Path.GetFileNameWithoutExtension(filePath)}{skeleton.bones.Count}";
                        resource.Model.fileName = resFile + ".P";
                        resource.res_file = resFile;

                        AssignRealGID(ref resource.Model);
                        ComputeBoundingBox(ref resource.Model);
                        CreateDListsFromPModel(ref resource.Model);

                        bone.fRSDResources.Add(resource);
                        bone.nResources++;
                    }
                }

                skeleton.bones.Add(bone);
            }

            skeleton.nBones = skeleton.bones.Count;
            return skeleton;
        }

        /// <summary>
        /// Converts a System.Numerics.Quaternion to FF7 Euler angles (alpha, beta, gamma).
        /// This reverses the EulerToQuaternionFF7 conversion.
        /// </summary>
        /// <param name="quat">The quaternion to convert</param>
        /// <param name="alpha">Output X rotation in degrees</param>
        /// <param name="beta">Output Y rotation in degrees</param>
        /// <param name="gamma">Output Z rotation in degrees</param>
        private static void QuaternionToEulerFF7(System.Numerics.Quaternion quat, out float alpha, out float beta, out float gamma)
        {
            // Convert System.Numerics.Quaternion to OpenTK.Mathematics.Quaterniond
            var quatD = new Quaterniond(quat.X, quat.Y, quat.Z, quat.W);

            // Use existing utility function to convert to matrix, then extract Euler
            double[] mat = new double[16];
            BuildMatrixFromQuaternion(quatD, ref mat);

            // GetEulerYXZrFromMatrix returns Vector3 where:
            // X = beta (Y rotation), Y = alpha (X rotation), Z = gamma (Z rotation)
            var euler = GetEulerYXZrFromMatrix(mat);

            alpha = (float)euler.Y;  // X rotation
            beta = (float)euler.X;   // Y rotation
            gamma = (float)euler.Z;  // Z rotation
        }

        /// <summary>
        /// Extracts a FieldAnimation from an Assimp Scene.
        /// This reverses the AddFieldAnimationToScene export transformation.
        /// </summary>
        /// <param name="scene">The Assimp scene containing animation data</param>
        /// <param name="skeleton">The FieldSkeleton to use for bone name matching</param>
        /// <param name="filePath">File path for naming the animation</param>
        /// <param name="bAdjust">Whether coordinate adjustments were applied during export</param>
        /// <returns>A FieldAnimation extracted from the scene, or an empty animation if none found</returns>
        public static FieldAnimation ExtractFieldAnimationFromScene(Scene scene, FieldSkeleton skeleton, string filePath, bool bAdjust = true)
        {
            var animation = new FieldAnimation
            {
                version = 1,
                nBones = skeleton.nBones,
                nFrames = 0,
                rotationOrder = new byte[] { 1, 0, 2 },  // YXZ order
                unused = 0,
                runtime_data = new int[5],
                frames = new List<FieldFrame>(),
                strFieldAnimationFile = Path.GetFileName(filePath).ToUpper()
            };

            // Check if scene has animations
            if (!scene.HasAnimations || scene.AnimationCount == 0)
                return animation;

            // Get the first animation (typically there's only one)
            var sceneAnim = scene.Animations[0];

            // Find the Armature channel (contains root transforms)
            NodeAnimationChannel? armatureChannel = null;
            foreach (var channel in sceneAnim.NodeAnimationChannels)
            {
                if (channel.NodeName == "Armature")
                {
                    armatureChannel = channel;
                    break;
                }
            }

            if (armatureChannel == null)
                return animation;

            // Build a map of bone names to channel indices for quick lookup
            var boneChannelMap = new Dictionary<string, NodeAnimationChannel>();
            foreach (var channel in sceneAnim.NodeAnimationChannels)
            {
                if (channel.NodeName != "Armature")
                {
                    boneChannelMap[channel.NodeName] = channel;
                }
            }

            // Determine frame count from armature channel
            int frameCount = armatureChannel.PositionKeyCount;
            animation.nFrames = frameCount;

            // Direction multipliers for reversing export transformations
            // Export: xDir = bAdjust ? -1 : 1, yDir = bAdjust ? 1 : -1, zDir = bAdjust ? 1 : -1
            float xDir = bAdjust ? -1 : 1;
            float yDir = bAdjust ? 1 : -1;
            float zDir = bAdjust ? 1 : -1;

            // 180° Z rotation for removing from quaternion if bAdjust
            var rot180ZInverse = System.Numerics.Quaternion.CreateFromAxisAngle(
                System.Numerics.Vector3.UnitZ, -(float)Math.PI);

            // Process each frame
            for (int fi = 0; fi < frameCount; fi++)
            {
                var frame = new FieldFrame
                {
                    rotations = new List<FieldRotation>()
                };

                // Extract root translation (reverse the direction multipliers)
                if (fi < armatureChannel.PositionKeyCount)
                {
                    var posKey = armatureChannel.PositionKeys[fi];
                    // Divide by direction multipliers to reverse (equivalent to multiply by same value)
                    frame.rootTranslationX = posKey.Value.X * xDir;
                    frame.rootTranslationY = posKey.Value.Y * yDir;
                    frame.rootTranslationZ = posKey.Value.Z * zDir;
                }

                // Extract root rotation
                if (fi < armatureChannel.RotationKeyCount)
                {
                    var rotKey = armatureChannel.RotationKeys[fi];
                    var quat = rotKey.Value;

                    // If bAdjust, remove the 180° Z rotation that was applied during export
                    if (bAdjust)
                    {
                        // Export did: rootQuat = Concatenate(rootQuat, rot180Z)
                        // To reverse: multiply by inverse of rot180Z
                        quat = System.Numerics.Quaternion.Concatenate(quat, rot180ZInverse);
                    }

                    // Convert quaternion back to Euler
                    QuaternionToEulerFF7(quat, out float alpha, out float beta, out float gamma);

                    // If bAdjust, reverse the sign adjustments on alpha and gamma
                    if (bAdjust)
                    {
                        alpha = -alpha;  // Was negated during export
                        gamma = -gamma;  // Was negated during export
                    }

                    frame.rootRotationAlpha = alpha;
                    frame.rootRotationBeta = beta;
                    frame.rootRotationGamma = gamma;
                }

                // Extract bone rotations
                for (int bi = 0; bi < skeleton.nBones; bi++)
                {
                    var bone = skeleton.bones[bi];
                    var rotation = new FieldRotation(0, 0, 0);

                    // Find the channel for this bone by joint_i name
                    if (boneChannelMap.TryGetValue(bone.joint_i, out var boneChannel))
                    {
                        if (fi < boneChannel.RotationKeyCount)
                        {
                            var rotKey = boneChannel.RotationKeys[fi];
                            var quat = rotKey.Value;

                            // Convert quaternion back to Euler
                            QuaternionToEulerFF7(quat, out float alpha, out float beta, out float gamma);

                            // If bAdjust, reverse the sign adjustments on alpha and gamma
                            if (bAdjust)
                            {
                                alpha = -alpha;  // Was negated during export
                                gamma = -gamma;  // Was negated during export
                            }

                            rotation = new FieldRotation(alpha, beta, gamma);
                        }
                    }

                    frame.rotations.Add(rotation);
                }

                animation.frames.Add(frame);
            }

            return animation;
        }

        /// <summary>
        /// Detects if an Assimp Scene represents a battle location (flat hierarchy with Piece_ nodes)
        /// versus a battle skeleton (hierarchical with Bone_ nodes).
        /// </summary>
        /// <param name="scene">The Assimp Scene to check</param>
        /// <returns>True if the scene contains battle location structure</returns>
        private static bool IsBattleLocationScene(Scene scene)
        {
            if (scene.RootNode == null)
                return false;

            // Find the armature node
            var armatureNode = scene.RootNode.Children.FirstOrDefault(c => c.Name == "Armature");
            if (armatureNode == null)
                return false;

            // Check for Piece_ nodes (battle location) vs Bone_ nodes (battle skeleton)
            bool hasPieceNodes = armatureNode.Children.Any(c => c.Name.StartsWith("Piece_"));
            bool hasBoneNodes = armatureNode.Children.Any(c => c.Name.StartsWith("Bone_"));

            // If we have Piece_ nodes and no Bone_ nodes, it's a battle location
            return hasPieceNodes && !hasBoneNodes;
        }

        /// <summary>
        /// Converts a battle location scene to a BattleSkeleton.
        /// Battle locations have flat hierarchy where each piece is independent.
        /// </summary>
        private static BattleSkeleton ConvertSceneToBattleLocation(Scene scene, BattleSkeleton skeleton,
            Node? armatureNode, List<string> tex_names, bool bAdjust)
        {
            if (armatureNode == null)
            {
                skeleton.nBones = 0;
                skeleton.nJoints = 0;
                return skeleton;
            }

            // Sort piece nodes by their index (Piece_0, Piece_1, etc.)
            var pieceNodes = armatureNode.Children
                .Where(c => c.Name.StartsWith("Piece_"))
                .Select(c => {
                    int index = 0;
                    if (c.Name.Length > 6)
                        int.TryParse(c.Name.Substring(6), out index);
                    return (node: c, index);
                })
                .OrderBy(x => x.index)
                .ToList();

            int pieceIndex = 0;
            foreach (var (pieceNode, _) in pieceNodes)
            {
                var bone = new BattleBone
                {
                    Models = new List<PModel>(),
                    nModels = 0,
                    hasModel = 0,
                    // For battle locations, parentBone is just the piece index (not a hierarchy reference)
                    parentBone = pieceIndex,
                    len = 0,
                    resizeX = 1.0f,
                    resizeY = 1.0f,
                    resizeZ = 1.0f
                };

                // Find Model_ children under this piece
                var modelNodes = pieceNode.Children
                    .Where(c => c.Name.StartsWith("Model_") && c.HasMeshes)
                    .ToList();

                // If piece node itself has meshes, include it
                if (pieceNode.HasMeshes)
                    modelNodes.Insert(0, pieceNode);

                foreach (var modelNode in modelNodes)
                {
                    var model = new PModel();
                    InitializePModel(ref model);

                    // Extract model transform from the node
                    // Export order was: BoneScale * Translate * RotX * RotY * RotZ * ModelScale
                    // We need to decompose and reverse the bAdjust transformations
                    var modelTransform = ToMatrix4(modelNode.Transform);

                    // Decompose the transform matrix to extract position, rotation, scale
                    // Note: This is a simplified extraction - for complex transforms,
                    // we extract what we can from the matrix
                    Vector3 translation = modelTransform.ExtractTranslation();
                    Vector3 scale = modelTransform.ExtractScale();

                    // Reverse bAdjust coordinate transformation for position
                    model.repositionX = bAdjust ? -translation.X : translation.X;
                    model.repositionY = translation.Y;
                    model.repositionZ = bAdjust ? -translation.Z : translation.Z;

                    // Extract rotation (approximate - assumes standard rotation order)
                    // For more accurate results, we'd need to properly decompose the rotation
                    var rotationMatrix = new Matrix3(modelTransform);
                    // Normalize by scale to get pure rotation
                    if (scale.X != 0) rotationMatrix.Row0 /= scale.X;
                    if (scale.Y != 0) rotationMatrix.Row1 /= scale.Y;
                    if (scale.Z != 0) rotationMatrix.Row2 /= scale.Z;

                    // Extract Euler angles (X-Y-Z order)
                    // This is an approximation and may have gimbal lock issues
                    float rotAlpha = (float)MathHelper.RadiansToDegrees(Math.Atan2(-rotationMatrix.M23, rotationMatrix.M33));
                    float rotBeta = (float)MathHelper.RadiansToDegrees(Math.Asin(rotationMatrix.M13));
                    float rotGamma = (float)MathHelper.RadiansToDegrees(Math.Atan2(-rotationMatrix.M12, rotationMatrix.M11));

                    // Reverse bAdjust for rotations
                    model.rotateAlpha = bAdjust ? -rotAlpha : rotAlpha;
                    model.rotateBeta = rotBeta;
                    model.rotateGamma = bAdjust ? -rotGamma : rotGamma;

                    model.resizeX = scale.X;
                    model.resizeY = scale.Y;
                    model.resizeZ = scale.Z;

                    // Process meshes
                    foreach (var meshIndex in modelNode.MeshIndices)
                    {
                        if (meshIndex < scene.MeshCount)
                        {
                            var mesh = scene.Meshes[meshIndex];
                            Material material = new() { Name = "Default" };
                            int texId = -1;
                            if (scene.HasMaterials && mesh.MaterialIndex >= 0 && mesh.MaterialIndex < scene.MaterialCount)
                            {
                                material = scene.Materials[mesh.MaterialIndex];
                                foreach (var tex in material.GetAllMaterialTextures())
                                {
                                    if (!string.IsNullOrEmpty(tex.FilePath))
                                    {
                                        string name = Path.GetFileNameWithoutExtension(tex.FilePath);
                                        int index = tex_names.IndexOf(name);
                                        if (index >= 0)
                                        {
                                            texId = index;
                                            break;
                                        }
                                    }
                                }
                            }
                            AddMeshAsGroup(mesh, material, ref model, bAdjust, texId);
                        }
                    }

                    if (model.Header.numGroups > 0)
                    {
                        model.fileName = $"Piece{pieceIndex}_Model{bone.nModels}.P";

                        AssignRealGID(ref model);
                        ComputeBoundingBox(ref model);
                        CreateDListsFromPModel(ref model);

                        bone.Models.Add(model);
                        bone.nModels++;
                        bone.hasModel = 1;

                        // Set bone length from model diameter (like LoadBattleLocationPiece does)
                        if (bone.len == 0)
                        {
                            bone.len = ComputeDiameter(model.BoundingBox) / 2;
                        }
                    }
                }

                skeleton.bones.Add(bone);
                pieceIndex++;
            }

            skeleton.nBones = skeleton.bones.Count;
            skeleton.nJoints = skeleton.bones.Count;
            skeleton.nWeapons = 0; // Battle locations don't have weapons

            return skeleton;
        }

        /// <summary>
        /// Converts an Assimp Scene to a BattleSkeleton.
        /// </summary>
        /// <param name="scene">The Assimp Scene to convert</param>
        /// <param name="skeletonType">The skeleton type to create</param>
        /// <param name="filePath">The path of the file being loaded</param>
        /// <param name="bAdjust">Whether to apply FF7 coordinate adjustments</param>
        /// <returns>A BattleSkeleton constructed from the scene</returns>
        public static BattleSkeleton ConvertSceneToBattleSkeleton(Scene scene, string filePath, bool bAdjust = true)
        {
            //get texture files (if they exist)
            var textures = new List<TEX>();
            var tex_names = new List<string>();
            var directory = Path.GetDirectoryName(filePath);
            if (Directory.Exists(directory))
            {
                textures = LoadTexturesFromScene(scene, directory);
                tex_names =
                    (from t in textures
                     select Path.GetFileNameWithoutExtension(t.TEXfileName)).ToList();
            }

            //get skeleton
            var skeleton = new BattleSkeleton
            {
                fileName = Path.GetFileName(filePath).ToUpper(),
                skeletonType = SkeletonType.PC,
                bones = new List<BattleBone>(),
                wpModels = new List<PModel>(),
                textures = textures,
                unk1 = 1,
                unk2 = 1,
                unk3 = 0,
                unk4 = 0,
                unk5 = 0,
                unk6 = 0,
                IsBattleLocation = false,
                CanHaveLimitBreak = false
            };

            if (scene.RootNode == null)
            {
                skeleton.nBones = 0;
                skeleton.TexIDS = Array.Empty<uint>();
                return skeleton;
            }

            skeleton.nTextures = skeleton.textures.Count;
            skeleton.TexIDS =
                (from tex in textures
                 select tex.texID).ToArray();

            // Find the armature node if it exists (created during export as container for bones)
            var armatureNode = scene.RootNode.Children.FirstOrDefault(c => c.Name == "Armature");

            // Detect if this is a battle location (flat hierarchy with Piece_ nodes)
            bool isBattleLocation = IsBattleLocationScene(scene);
            skeleton.IsBattleLocation = isBattleLocation;

            if (isBattleLocation)
            {
                // Battle location import: flat hierarchy with Piece_ nodes
                return ConvertSceneToBattleLocation(scene, skeleton, armatureNode, tex_names, bAdjust);
            }

            // Regular battle skeleton import: hierarchical with Bone_ nodes
            // Flatten hierarchy and assign indices
            var flatNodes = FlattenHierarchy(scene.RootNode);

            // Build node-to-index map (excluding root, armature, and weapons)
            var nodeToIndex = new Dictionary<Node, int>();
            var boneNodesList = new List<(Node node, Node? parent)>();

            foreach (var (node, parent, depth) in flatNodes)
            {
                // Skip root node, armature node, weapon armature node, and weapons node
                if (node == scene.RootNode || node == armatureNode ||
                    node.Name == "WeaponArmature" || node.Name == "Weapons")
                    continue;

                // Skip weapon children (under Weapons or WeaponArmature)
                if (parent != null && (parent.Name == "Weapons" || parent.Name == "WeaponArmature"))
                    continue;

                // Skip model child nodes (they contain mesh data, not bone data)
                if (node.Name.StartsWith("Model_"))
                    continue;

                // Check if this is a bone node (has bone metadata or has meshes/children)
                if (node.Name.StartsWith("Bone_") ||
                    node.Metadata.ContainsKey("FF7_ParentIndex") ||
                    node.HasMeshes ||
                    node.ChildCount > 0)
                {
                    nodeToIndex[node] = boneNodesList.Count;
                    boneNodesList.Add((node, parent));
                }
            }

            // Process bone nodes
            foreach (var (node, parent) in boneNodesList)
            {
                var bone = new BattleBone
                {
                    Models = new List<PModel>(),
                    nModels = 0,
                    hasModel = 0,
                    resizeX = 1.0f,
                    resizeY = 1.0f,
                    resizeZ = 1.0f
                };

                // Get parent index
                //var storedParent = GetFF7MetadataInt(node, "FF7_ParentIndex");
                //if (storedParent.HasValue)
                //{
                //    bone.parentBone = storedParent.Value;
                //}
                //else
                if (parent != null && nodeToIndex.TryGetValue(parent, out int parentIndex))
                {
                    bone.parentBone = parentIndex;
                }
                else
                {
                    bone.parentBone = -1; // Root
                }

                // Calculate bone length - tries metadata first, then infers from child transforms
                bone.len = CalculateBattleBoneLengthFromNode(node, bAdjust);

                // Get resize values
                // bone.resizeX = GetFF7MetadataFloat(node, "FF7_ResizeX") ?? 1.0f;
                // bone.resizeY = GetFF7MetadataFloat(node, "FF7_ResizeY") ?? 1.0f;
                // bone.resizeZ = GetFF7MetadataFloat(node, "FF7_ResizeZ") ?? 1.0f;

                // Process meshes
                var meshNodes = new List<Node> { node };
                meshNodes.AddRange(node.Children.Where(c => c.HasMeshes && c.Name.StartsWith("Model_")));

                foreach (var meshNode in meshNodes)
                {
                    if (!meshNode.HasMeshes) continue;

                    var model = new PModel();
                    InitializePModel(ref model);

                    foreach (var meshIndex in meshNode.MeshIndices)
                    {
                        if (meshIndex < scene.MeshCount)
                        {
                            var mesh = scene.Meshes[meshIndex];
                            Material material = new() { Name = "Default" };
                            int texId = -1;
                            if (scene.HasMaterials && mesh.MaterialIndex >= 0 && mesh.MaterialIndex < scene.MaterialCount)
                            {
                                material = scene.Materials[mesh.MaterialIndex];
                                foreach (var tex in material.GetAllMaterialTextures())
                                {
                                    if (!string.IsNullOrEmpty(tex.FilePath))
                                    {
                                        string name = Path.GetFileNameWithoutExtension(tex.FilePath);
                                        int index = tex_names.IndexOf(name);
                                        if (index >= 0)
                                        {
                                            texId = index;
                                            break;
                                        }
                                    }
                                }
                            }
                            AddMeshAsGroup(mesh, material, ref model, bAdjust, texId);
                        }
                    }

                    if (model.Header.numGroups > 0)
                    {
                        // Extract model transform
                        if (meshNode != node)
                        {
                            var modelTransform = ToMatrix4(meshNode.Transform);
                            model.repositionX = modelTransform.M41;
                            model.repositionY = modelTransform.M42;
                            model.repositionZ = modelTransform.M43;
                        }

                        // Set fileName BEFORE adding to list (PModel is a struct, so it's copied)
                        model.fileName = $"Bone{skeleton.bones.Count}_Model{bone.nModels}.P";

                        // Process model before adding (struct is copied when added to list)
                        AssignRealGID(ref model);
                        ComputeBoundingBox(ref model);
                        CreateDListsFromPModel(ref model);

                        bone.Models.Add(model);
                        bone.nModels++;
                        bone.hasModel = 1;
                    }
                }

                skeleton.bones.Add(bone);
            }

            // Process weapons node
            // Current hierarchy: Armature -> WeaponArmature -> Weapons -> Weapon_N
            Node? weaponsNode = null;
            Node? weaponArmatureNode = null;

            // First, try to find WeaponArmature as child of Armature
            if (armatureNode != null)
            {
                weaponArmatureNode = armatureNode.Children.FirstOrDefault(c => c.Name == "WeaponArmature");
                if (weaponArmatureNode != null)
                {
                    weaponsNode = weaponArmatureNode.Children.FirstOrDefault(c => c.Name == "Weapons");

                    if (weaponsNode != null)
                    {
                        // Sort weapon children by their index (Weapon_0, Weapon_1, etc.)
                        // This ensures weapons are imported in the correct order regardless of
                        // how Assimp orders the children in the loaded scene.
                        var sortedWeaponChildren = weaponsNode.Children
                            .Where(c => c.HasMeshes && c.Name.StartsWith("Weapon_"))
                            .Select(c => {
                                int index = 0;
                                if (c.Name.Length > 7)
                                    int.TryParse(c.Name.Substring(7), out index);
                                return (node: c, index);
                            })
                            .OrderBy(x => x.index)
                            .ToList();

                        foreach (var (weaponChild, weaponIndex) in sortedWeaponChildren)
                        {
                            var wpModel = new PModel();
                            InitializePModel(ref wpModel);

                            foreach (var meshIndex in weaponChild.MeshIndices)
                            {
                                if (meshIndex < scene.MeshCount)
                                {
                                    var mesh = scene.Meshes[meshIndex];
                                    Material material = new() { Name = "Default" };
                                    int texId = -1;
                                    if (scene.HasMaterials && mesh.MaterialIndex >= 0 && mesh.MaterialIndex < scene.MaterialCount)
                                    {
                                        material = scene.Materials[mesh.MaterialIndex];
                                        foreach (var tex in material.GetAllMaterialTextures())
                                        {
                                            if (!string.IsNullOrEmpty(tex.FilePath))
                                            {
                                                string name = Path.GetFileNameWithoutExtension(tex.FilePath);
                                                int index = tex_names.IndexOf(name);
                                                if (index >= 0)
                                                {
                                                    texId = index;
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                    AddMeshAsGroup(mesh, material, ref wpModel, bAdjust, texId);
                                }
                            }

                            // Set unique fileName for mesh caching (required for GLRenderer)
                            wpModel.fileName = $"Weapon_{weaponIndex}.P";

                            // Process model before adding (struct is copied when added to list)
                            AssignRealGID(ref wpModel);
                            ComputeBoundingBox(ref wpModel);
                            CreateDListsFromPModel(ref wpModel);

                            skeleton.wpModels.Add(wpModel);
                        }
                    }
                }
            }

            skeleton.nBones = skeleton.bones.Count;
            skeleton.nWeapons = skeleton.wpModels.Count;
            skeleton.nJoints = skeleton.nBones; // Simplified

            return skeleton;
        }

        /// <summary>
        /// Extracts a BattleAnimation from an Assimp Scene.
        /// This reverses the AddBattleAnimationToScene export transformation.
        /// </summary>
        /// <param name="scene">The Assimp scene containing animation data</param>
        /// <param name="skeleton">The BattleSkeleton to use for bone matching</param>
        /// <param name="animationsPack">The BattleAnimationsPack to put data into</param>
        /// <param name="animIndex">Where to put the extracted animation</param>
        /// <param name="bAdjust">Whether coordinate adjustments were applied during export</param>
        public static void ExtractBattleAnimationsFromScene(Scene scene,
           ref BattleSkeleton skeleton, ref BattleAnimationsPack animationsPack, int animIndex,
           bool insertNew, bool bAdjust = true)
        {
            /*if (animIndex < 0 || animIndex > animationsPack.nAnimations)
                throw new ArgumentOutOfRangeException(nameof(animIndex));

            if (insertNew) //add a new animation
            {
                skeleton.nsSkeletonAnims++;
                animationsPack.nbSkeletonAnims++;
                animationsPack.nAnimations++;
            }
            else if (animIndex == animationsPack.nbSkeletonAnims)
                throw new ArgumentOutOfRangeException(nameof(animIndex));*/

            var battleAnimation = new BattleAnimation
                {
                    nBones = skeleton.nBones > 1 ? skeleton.nBones + 1 : 1,  // +1 for root transformation if multi-bone
                    numFrames = 0,
                    numFramesShort = 0,
                    blockSize = 0,
                    blockSizeShort = 0,
                    key = 0,
                    framesRawData = null,
                    padding4bytes = null,
                    frames = new List<BattleFrame>()
                };

            // Check if scene has animations
            if (scene.HasAnimations && scene.AnimationCount > 0)
            {
                // Get the first animation
                var sceneAnim = scene.Animations[0];

                // Find the Armature channel (contains root position and rotation)
                NodeAnimationChannel? armatureChannel = null;
                foreach (var channel in sceneAnim.NodeAnimationChannels)
                {
                    if (channel.NodeName == "Armature")
                    {
                        armatureChannel = channel;
                        break;
                    }
                }

                if (armatureChannel != null)
                {
                    // Build a map of bone names to channels for quick lookup
                    var boneChannelMap = new Dictionary<string, NodeAnimationChannel>();
                    foreach (var channel in sceneAnim.NodeAnimationChannels)
                    {
                        if (channel.NodeName.StartsWith("Bone_"))
                        {
                            boneChannelMap[channel.NodeName] = channel;
                        }
                    }

                    // Determine frame count from armature channel
                    int frameCount = armatureChannel.PositionKeyCount;
                    battleAnimation.numFrames = frameCount;
                    battleAnimation.numFramesShort = (ushort)frameCount;

                    // Direction multipliers for reversing export transformations
                    // Export used: xDir = bAdjust ? -1 : 1, zDir = bAdjust ? -1 : 1
                    float xDir = bAdjust ? -1 : 1;
                    float zDir = bAdjust ? -1 : 1;

                    // 180° Z rotation inverse for removing from quaternion if bAdjust
                    var rot180ZInverse = System.Numerics.Quaternion.CreateFromAxisAngle(
                        System.Numerics.Vector3.UnitZ, -(float)Math.PI);

                    // Bone index offset: if skeleton has >1 bone, animation bone[0] is root rotation
                    // and bone[1+] are actual bones; otherwise bone[0] is both root and only bone
                    int itmpbones = skeleton.nBones > 1 ? 1 : 0;

                    // Process each frame
                    for (int fi = 0; fi < frameCount; fi++)
                    {
                        var frame = new BattleFrame
                        {
                            startX = 0,
                            startY = 0,
                            startZ = 0,
                            bones = new List<BattleFrameBone>()
                        };

                        // Extract root translation from armature position
                        if (fi < armatureChannel.PositionKeyCount)
                        {
                            var posKey = armatureChannel.PositionKeys[fi];
                            // Reverse direction multipliers
                            frame.startX = (int)(posKey.Value.X * xDir);
                            frame.startY = (int)posKey.Value.Y;  // Y was not negated in export
                            frame.startZ = (int)(posKey.Value.Z * zDir);
                        }

                        // Extract root rotation from armature rotation (goes into bone[0])
                        BattleFrameBone rootBone = new BattleFrameBone(0, 0, 0);
                        if (fi < armatureChannel.RotationKeyCount)
                        {
                            var rotKey = armatureChannel.RotationKeys[fi];
                            var quat = rotKey.Value;

                            // If bAdjust, remove the 180° Z rotation that was applied during export
                            if (bAdjust)
                            {
                                quat = System.Numerics.Quaternion.Concatenate(quat, rot180ZInverse);
                            }

                            // Convert quaternion back to Euler
                            QuaternionToEulerFF7(quat, out float alpha, out float beta, out float gamma);

                            // If bAdjust, reverse the sign adjustments on alpha and gamma
                            if (bAdjust)
                            {
                                alpha = -alpha;
                                gamma = -gamma;
                            }

                            // BattleFrameBone stores degrees (QuaternionToEulerFF7 returns degrees)
                            rootBone = new BattleFrameBone(alpha, beta, gamma);
                        }
                        frame.bones.Add(rootBone);

                        // Extract rotations for each skeleton bone
                        for (int bi = 0; bi < skeleton.nBones; bi++)
                        {
                            var boneName = $"Bone_{bi}";
                            var boneRotation = new BattleFrameBone(0, 0, 0);

                            if (boneChannelMap.TryGetValue(boneName, out var boneChannel))
                            {
                                if (fi < boneChannel.RotationKeyCount)
                                {
                                    var rotKey = boneChannel.RotationKeys[fi];
                                    var quat = rotKey.Value;

                                    // Convert quaternion back to Euler
                                    QuaternionToEulerFF7(quat, out float alpha, out float beta, out float gamma);

                                    // If bAdjust, reverse the sign adjustments on alpha and gamma
                                    if (bAdjust)
                                    {
                                        alpha = -alpha;
                                        gamma = -gamma;
                                    }

                                    // BattleFrameBone stores degrees (QuaternionToEulerFF7 returns degrees)
                                    boneRotation = new BattleFrameBone(alpha, beta, gamma);
                                }
                            }

                            frame.bones.Add(boneRotation);
                        }

                        battleAnimation.frames.Add(frame);
                    }

                    // Insert the extracted animation into the battle animations pack
                    if (insertNew)
                        animationsPack.SkeletonAnimations.Insert(animIndex, battleAnimation);
                    else
                        animationsPack.SkeletonAnimations[animIndex] = battleAnimation;

                    // Extract the weapon animation if it exists
                    var weaponAnimation = ExtractWeaponAnimation(sceneAnim, battleAnimation, bAdjust);
                    if (weaponAnimation != null)
                    {
                        var a = (BattleAnimation)weaponAnimation;
                        if (insertNew)
                            animationsPack.WeaponAnimations.Insert(animIndex, a);
                        else
                            animationsPack.WeaponAnimations[animIndex] = a;
                    }
                }
            }
        }

        /// <summary>
        /// Extracts a weapon BattleAnimation from an Assimp Scene.
        /// The weapon animation is stored as part of the main animation (WeaponArmature channel).
        /// Position is stored relative to skeleton, so we need skeleton animation to convert back to absolute.
        /// </summary>
        /// <param name="sceneAnim">The Assimp animation containing the weapon animation</param>
        /// <param name="skeletonAnimation">The skeleton animation (needed to convert relative positions to absolute)</param>
        /// <param name="bAdjust">Whether coordinate adjustments were applied during export</param>
        /// <returns>A BattleAnimation for the weapon, or null if no weapon animation found</returns>
        private static BattleAnimation? ExtractWeaponAnimation(Animation sceneAnim, BattleAnimation skeletonAnimation, bool bAdjust = true)
        {
            NodeAnimationChannel? weaponArmatureChannel = null;
            foreach (var channel in sceneAnim.NodeAnimationChannels)
            {
                if (channel.NodeName == "WeaponArmature")
                {
                    weaponArmatureChannel = channel;
                    break;
                }
            }

            if (weaponArmatureChannel == null)
                return null;

            var weaponAnimation = new BattleAnimation
            {
                nBones = 1,  // Weapon animation only has one bone (the weapon rotation)
                numFrames = 0,
                numFramesShort = 0,
                blockSize = 0,
                blockSizeShort = 0,
                key = 0,
                framesRawData = null,
                padding4bytes = null,
                frames = new List<BattleFrame>()
            };

            // Determine frame count from weapon armature channel
            int frameCount = weaponArmatureChannel.PositionKeyCount;
            weaponAnimation.numFrames = frameCount;
            weaponAnimation.numFramesShort = (ushort)frameCount;

            // Direction multipliers for reversing export transformations
            // Export used: relX * xDir, relY, relZ * zDir (X and Z negated when bAdjust)
            float xDir = bAdjust ? -1 : 1;
            float zDir = bAdjust ? -1 : 1;

            // Process each frame
            for (int fi = 0; fi < frameCount; fi++)
            {
                var frame = new BattleFrame
                {
                    startX = 0,
                    startY = 0,
                    startZ = 0,
                    bones = new List<BattleFrameBone>()
                };

                // Extract weapon translation (relative to skeleton)
                // Then convert back to absolute by adding skeleton position
                if (fi < weaponArmatureChannel.PositionKeyCount)
                {
                    var posKey = weaponArmatureChannel.PositionKeys[fi];

                    // Reverse the coordinate transforms to get relative position in FF7 space
                    float relX = posKey.Value.X * xDir;
                    float relY = posKey.Value.Y;
                    float relZ = posKey.Value.Z * zDir;

                    // Add skeleton position to get absolute weapon position
                    if (fi < skeletonAnimation.frames.Count)
                    {
                        var skelFrame = skeletonAnimation.frames[fi];
                        frame.startX = (int)relX + skelFrame.startX;
                        frame.startY = (int)relY + skelFrame.startY;
                        frame.startZ = (int)relZ + skelFrame.startZ;
                    }
                    else
                    {
                        frame.startX = (int)relX;
                        frame.startY = (int)relY;
                        frame.startZ = (int)relZ;
                    }
                }

                // Extract weapon rotation (goes into bones[0])
                // Note: Weapon rotation does NOT have 180° Z rotation applied (inherited from parent Armature)
                var weaponBone = new BattleFrameBone(0, 0, 0);
                if (fi < weaponArmatureChannel.RotationKeyCount)
                {
                    var rotKey = weaponArmatureChannel.RotationKeys[fi];
                    var quat = rotKey.Value;

                    // Convert quaternion back to Euler
                    QuaternionToEulerFF7(quat, out float alpha, out float beta, out float gamma);

                    // If bAdjust, reverse the sign adjustments on alpha and gamma
                    if (bAdjust)
                    {
                        alpha = -alpha;
                        gamma = -gamma;
                    }

                    weaponBone = new BattleFrameBone(alpha, beta, gamma);
                }
                frame.bones.Add(weaponBone);

                weaponAnimation.frames.Add(frame);
            }

            return weaponAnimation;
        }

        /// <summary>
        /// Initializes a new empty PModel with default values.
        /// </summary>
        private static void InitializePModel(ref PModel model)
        {
            model.Header = new PHeader
            {
                version = 1,
                off04 = 1,
                numGroups = 0,
                unknown = new int[16]
            };
            model.Groups = Array.Empty<PGroup>();
            model.Verts = Array.Empty<Vector3>();
            model.Normals = Array.Empty<Vector3>();
            model.Polys = null;
            model.TexCoords = Array.Empty<Vector2>();
            model.Vcolors = Array.Empty<Color>();
            model.Pcolors = Array.Empty<Color>();
            model.Hundrets = Array.Empty<PHundret>();

            model.resizeX = 1;
            model.resizeY = 1;
            model.resizeZ = 1;
            model.repositionX = 0;
            model.repositionY = 0;
            model.repositionZ = 0;
            model.rotateAlpha = 0;
            model.rotateBeta = 0;
            model.rotateGamma = 0;
            model.rotationQuaternion = new Quaterniond(0, 0, 0, 1);
        }

        #endregion

        #region EXPORT FUNCTIONS

        /// <summary>
        /// Gets the Assimp format ID from a file extension.
        /// </summary>
        private static string GetExportFormatFromExtension(string extension)
        {
            extension = extension.ToLowerInvariant();

            // Map common extensions to Assimp format IDs
            return extension switch
            {
                ".obj" => "obj",
                ".objnomtl" => "objnomtl",
                ".stl" => "stl",
                ".stlb" => "stlb",
                ".ply" => "ply",
                ".plyb" => "plyb",
                ".dae" => "collada",
                ".gltf" => "gltf2",
                ".glb" => "glb2",
                ".fbx" => "fbx",
                ".3ds" => "3ds",
                ".x" => "x",
                _ => ""
            };
        }

        public static bool IsValidExport(string path)
        {
            string fileExt = Path.GetExtension(path).ToLower();
            return _context.IsExportFormatSupported(fileExt);
        }

        /// <summary>
        /// Gets a file filter string for export dialogs.
        /// </summary>
        public static string GetExportFileFilter()
        {
            var builder = new StringBuilder();

            // Common export formats
            builder.Append("Wavefront OBJ (*.obj)|*.obj|");
            builder.Append("Collada DAE (*.dae)|*.dae|");
            builder.Append("GLTF 2.0 (*.gltf)|*.gltf|");
            builder.Append("GLTF Binary (*.glb)|*.glb|");
            builder.Append("STL (*.stl)|*.stl|");
            builder.Append("PLY (*.ply)|*.ply|");
            builder.Append("3D Studio (*.3ds)|*.3ds|");
            builder.Append("FBX (*.fbx)|*.fbx");

            return builder.ToString();
        }

        /// <summary>
        /// Gets a file filter string for skeleton export dialogs.
        /// Only includes formats that support node hierarchies.
        /// </summary>
        public static string GetSkeletonExportFileFilter()
        {
            var builder = new StringBuilder();

            // Only formats that support hierarchies/skeletons
            builder.Append("GLTF 2.0 (*.gltf)|*.gltf|");
            builder.Append("GLTF Binary (*.glb)|*.glb|");
            builder.Append("Collada DAE (*.dae)|*.dae|");
            builder.Append("FBX (*.fbx)|*.fbx");

            return builder.ToString();
        }

        /// <summary>
        /// Stores an FF7-specific metadata value in an Assimp node.
        /// </summary>
        /*private static void SetFF7Metadata(Node node, string key, object value)
        {
            if (value is double d)
                node.Metadata[key] = new Metadata.Entry(MetaDataType.Double, d);
            else if (value is float f)
                node.Metadata[key] = new Metadata.Entry(MetaDataType.Float, f);
            else if (value is int i)
                node.Metadata[key] = new Metadata.Entry(MetaDataType.Int32, i);
            else if (value is string s)
                node.Metadata[key] = new Metadata.Entry(MetaDataType.String, s);
        }

        /// <summary>
        /// Retrieves an FF7-specific float metadata value from an Assimp node.
        /// </summary>
        private static float? GetFF7MetadataFloat(Node node, string key)
        {
            if (node.Metadata.TryGetValue(key, out var entry))
            {
                if (entry.Data is float f)
                    return f;
                if (entry.Data is double d)
                    return (float)d;
            }
            return null;
        }

        /// <summary>
        /// Retrieves an FF7-specific int metadata value from an Assimp node.
        /// </summary>
        private static int? GetFF7MetadataInt(Node node, string key)
        {
            if (node.Metadata.TryGetValue(key, out var entry))
            {
                if (entry.Data is int i)
                    return i;
            }
            return null;
        }

        /// <summary>
        /// Retrieves an FF7-specific string metadata value from an Assimp node.
        /// </summary>
        private static string? GetFF7MetadataString(Node node, string key)
        {
            if (node.Metadata.TryGetValue(key, out var entry))
            {
                if (entry.Data is string s)
                    return s;
            }
            return null;
        }*/

        /// <summary>
        /// Converts a PModel to an Assimp Scene for export to various formats.
        /// </summary>
        /// <param name="model">The PModel to convert</param>
        /// <param name="bAdjust">Whether to reverse FF7 coordinate adjustments (un-mirror Z, un-rotate 180° Y)</param>
        /// <param name="usePolygonColors">If true, uses polygon colors (field models). If false, uses vertex colors (battle models). Default is true.</param>
        /// <returns>An Assimp Scene ready for export</returns>
        public static Scene ConvertPModelToScene(PModel model, bool bAdjust, bool usePolygonColors = true)
        {
            var scene = new Scene();
            scene.RootNode = new Node("Root");

            // Create a mesh for each group in the PModel
            int meshIndex = 0;
            for (int g = 0; g < model.Header.numGroups; g++)
            {
                var group = model.Groups[g];

                // Skip hidden or empty groups
                if (group.HiddenQ || group.numVert == 0 || group.numPoly == 0)
                    continue;

                var mesh = ConvertGroupToMesh(model, g, bAdjust, usePolygonColors);
                mesh.Name = $"Group_{g}";

                // Assign material index
                mesh.MaterialIndex = meshIndex;

                scene.Meshes.Add(mesh);

                // Add mesh index directly to RootNode (simplest structure)
                scene.RootNode.MeshIndices.Add(meshIndex);

                // Create a material for this group
                var material = CreateMaterial(g);
                scene.Materials.Add(material);
                meshIndex++;
            }

            // If no meshes were added (all groups hidden/empty), add a default material
            if (scene.MeshCount == 0)
            {
                scene.Materials.Add(new Material { Name = "Default" });
            }

            return scene;
        }

        /// <summary>
        /// Converts a single PModel group to an Assimp Mesh.
        /// </summary>
        /// <param name="model">The PModel containing the group</param>
        /// <param name="groupIndex">Index of the group to convert</param>
        /// <param name="bAdjust">Whether to reverse FF7 coordinate adjustments</param>
        /// <param name="usePolygonColors">If true, duplicates vertices per polygon to preserve polygon colors (for field models).
        /// If false, shares vertices and uses vertex colors (for battle models).</param>
        private static Mesh ConvertGroupToMesh(PModel model, int groupIndex, bool bAdjust, bool usePolygonColors = true)
        {
            var group = model.Groups[groupIndex];
            var mesh = new Mesh(PrimitiveType.Triangle);

            int vertStart = group.offsetVert;
            int polyStart = group.offsetPoly;
            int texStart = group.offsetTex;

            bool hasTexCoords = model.TexCoords != null &&
                               model.TexCoords.Length > 0 &&
                               texStart >= 0;

            if (usePolygonColors)
            {
                // Polygon colors mode: duplicate vertices per polygon to preserve per-face colors
                // This is used for field models where each polygon has a single flat color
                ConvertGroupWithPolygonColors(model, group, mesh, vertStart, polyStart, texStart, hasTexCoords, bAdjust);
            }
            else
            {
                // Vertex colors mode: share vertices with original indexing
                // This is used for battle models where colors are interpolated across vertices
                ConvertGroupWithVertexColors(model, group, mesh, vertStart, polyStart, texStart, hasTexCoords, bAdjust);
            }

            // Set UV component count (2 = UV, not UVW) - required for Assimp export
            if (mesh.TextureCoordinateChannels[0].Count > 0)
            {
                mesh.UVComponentCount[0] = 2;
            }

            return mesh;
        }

        /// <summary>
        /// Converts a group using polygon colors (duplicates vertices per polygon).
        /// Used for field models where each polygon has a single flat color.
        /// </summary>
        private static void ConvertGroupWithPolygonColors(PModel model, PGroup group, Mesh mesh,
            int vertStart, int polyStart, int texStart, bool hasTexCoords, bool bAdjust)
        {
            bool hasPcolors = model.Pcolors != null && model.Pcolors.Length > 0;

            // Create vertices per-polygon to preserve per-polygon colors
            // Each polygon gets 3 unique vertices with that polygon's color
            for (int p = 0; p < group.numPoly; p++)
            {
                var poly = model.Polys[polyStart + p];

                // Get polygon color (or default to white)
                System.Numerics.Vector4 polyColor;
                if (hasPcolors && (polyStart + p) < model.Pcolors?.Length)
                {
                    var pc = model.Pcolors[polyStart + p];
                    polyColor = new System.Numerics.Vector4(pc.R / 255f, pc.G / 255f, pc.B / 255f, pc.A / 255f);
                }
                else
                {
                    polyColor = new System.Numerics.Vector4(1f, 1f, 1f, 1f); // Default white
                }

                // Process each vertex of the polygon
                // When bAdjust flips X and Z, reverse winding to maintain correct face orientation
                int[] windingOrder = bAdjust ? new[] { 0, 2, 1 } : new[] { 0, 1, 2 };

                for (int vi = 0; vi < 3; vi++)
                {
                    int polyVertIdx = windingOrder[vi];
                    int localVertIdx = poly.Verts[polyVertIdx];
                    int globalVertIdx = vertStart + localVertIdx;

                    // Add vertex position
                    AddVertexPosition(mesh, model, globalVertIdx, bAdjust);

                    // Add normal
                    AddVertexNormal(mesh, model, poly.Normals[polyVertIdx]);

                    // Add texture coordinate
                    AddVertexTexCoord(mesh, model, texStart + localVertIdx, hasTexCoords);

                    // Add vertex color - use polygon color directly to preserve per-face colors
                    mesh.VertexColorChannels[0].Add(polyColor);
                }

                // Add face with sequential indices (since we duplicated vertices)
                var face = new Face();
                int baseIdx = p * 3;
                face.Indices.Add(baseIdx);
                face.Indices.Add(baseIdx + 1);
                face.Indices.Add(baseIdx + 2);
                mesh.Faces.Add(face);
            }
        }

        /// <summary>
        /// Converts a group using vertex colors (shares vertices with original indexing).
        /// Used for battle models where colors are interpolated across vertices.
        /// </summary>
        private static void ConvertGroupWithVertexColors(PModel model, PGroup group, Mesh mesh,
            int vertStart, int polyStart, int texStart, bool hasTexCoords, bool bAdjust)
        {
            bool hasVcolors = model.Vcolors != null && model.Vcolors.Length > 0;

            // First pass: add all vertices for this group
            for (int v = 0; v < group.numVert; v++)
            {
                int globalVertIdx = vertStart + v;

                // Add vertex position
                AddVertexPosition(mesh, model, globalVertIdx, bAdjust);

                // Add normal (use the normal index for this vertex)
                int normalIdx = globalVertIdx < model.NormalIndex?.Length ? model.NormalIndex[globalVertIdx] : 0;
                AddVertexNormal(mesh, model, normalIdx);

                // Add texture coordinate
                AddVertexTexCoord(mesh, model, texStart + v, hasTexCoords);

                // Add vertex color
                System.Numerics.Vector4 vertColor;
                if (hasVcolors && globalVertIdx < model.Vcolors?.Length)
                {
                    var vc = model.Vcolors[globalVertIdx];
                    vertColor = new System.Numerics.Vector4(vc.R / 255f, vc.G / 255f, vc.B / 255f, vc.A / 255f);
                }
                else
                {
                    vertColor = new System.Numerics.Vector4(1f, 1f, 1f, 1f); // Default white
                }
                mesh.VertexColorChannels[0].Add(vertColor);
            }

            // Second pass: add faces with indices relative to this group's vertices
            for (int p = 0; p < group.numPoly; p++)
            {
                var poly = model.Polys[polyStart + p];

                var face = new Face();
                // When bAdjust flips X and Z, reverse winding to maintain correct face orientation
                if (bAdjust)
                {
                    face.Indices.Add(poly.Verts[0]);
                    face.Indices.Add(poly.Verts[2]);
                    face.Indices.Add(poly.Verts[1]);
                }
                else
                {
                    face.Indices.Add(poly.Verts[0]);
                    face.Indices.Add(poly.Verts[1]);
                    face.Indices.Add(poly.Verts[2]);
                }
                mesh.Faces.Add(face);
            }
        }

        /// <summary>
        /// Adds a vertex position to the mesh.
        /// </summary>
        private static void AddVertexPosition(Mesh mesh, PModel model, int globalVertIdx, bool bAdjust)
        {
            if (globalVertIdx < model.Verts.Length)
            {
                var vert = model.Verts[globalVertIdx];
                float x = vert.X;
                float y = vert.Y;
                float z = vert.Z;

                if (bAdjust)
                {
                    // Reverse FF7 adjustments: un-mirror X and Z planes
                    x = -x;
                    z = -z;
                }

                mesh.Vertices.Add(new System.Numerics.Vector3(x, y, z));
            }
            else
            {
                mesh.Vertices.Add(new System.Numerics.Vector3(0, 0, 0));
            }
        }

        /// <summary>
        /// Adds a vertex normal to the mesh.
        /// </summary>
        private static void AddVertexNormal(Mesh mesh, PModel model, int normalIdx)
        {
            if (normalIdx < model.Normals.Length)
            {
                var normal = model.Normals[normalIdx];
                mesh.Normals.Add(new System.Numerics.Vector3(normal.X, normal.Y, normal.Z));
            }
            else
            {
                mesh.Normals.Add(new System.Numerics.Vector3(0, 1, 0));
            }
        }

        /// <summary>
        /// Adds a texture coordinate to the mesh.
        /// </summary>
        private static void AddVertexTexCoord(Mesh mesh, PModel model, int texIdx, bool hasTexCoords)
        {
            if (hasTexCoords && model.TexCoords != null && texIdx < model.TexCoords.Length)
            {
                var tc = model.TexCoords[texIdx];
                mesh.TextureCoordinateChannels[0].Add(new System.Numerics.Vector3(tc.X, tc.Y, 0));
            }
            else
            {
                // Add default UV even when no texcoords - some exporters require this
                mesh.TextureCoordinateChannels[0].Add(new System.Numerics.Vector3(0, 0, 0));
            }
        }

        /// <summary>
        /// Exports a PModel to a file using Assimp.
        /// </summary>
        /// <param name="model">The PModel to export</param>
        /// <param name="filePath">Output file path (extension determines format)</param>
        /// <param name="bAdjust">Whether to reverse FF7 coordinate adjustments</param>
        /// <returns>True if export succeeded</returns>
        /// <summary>
        /// Exports a PModel to a file using Assimp.
        /// </summary>
        /// <param name="model">The PModel to export</param>
        /// <param name="filePath">Output file path (extension determines format)</param>
        /// <param name="bAdjust">Whether to reverse FF7 coordinate adjustments</param>
        /// <param name="usePolygonColors">If true, uses polygon colors (field models). If false, uses vertex colors (battle models). Default is true.</param>
        /// <returns>True if export succeeded</returns>
        public static bool ExportPModel(PModel model, string filePath, bool bAdjust, bool usePolygonColors = true)
        {
            string outputDir = Path.GetDirectoryName(filePath) ?? ".";
            string baseName = Path.GetFileNameWithoutExtension(filePath);

            var scene = ConvertPModelToScene(model, bAdjust, usePolygonColors);

            string formatId = GetExportFormatFromExtension(Path.GetExtension(filePath));
            if (string.IsNullOrEmpty(formatId))
            {
                throw new NotSupportedException(
                    $"Export format not supported for extension: {Path.GetExtension(filePath)}");
            }

            return _context.ExportFile(scene, filePath, formatId);
        }

        /// <summary>
        /// Exports textures from a list of TEX and returns a mapping of texture index to filename.
        /// </summary>
        /// <param name="textures">The list of textures to export</param>
        /// <param name="outputDirectory">Directory to save PNG files</param>
        /// <param name="baseFileName">Base name for texture files</param>
        /// <returns>Dictionary mapping texture index to saved filename</returns>
        private static Dictionary<int, string> ExportTexturesToPNG(
            IList<TEX> textures,
            string outputDirectory,
            string baseFileName)
        {
            var textureFiles = new Dictionary<int, string>();

            for (int i = 0; i < textures.Count; i++)
            {
                var tex = textures[i];
                if (tex.pixelData == null || tex.pixelData.Length == 0)
                    continue;

                // Use original TEX filename if available, otherwise generate one
                string texName = !string.IsNullOrEmpty(tex.TEXfileName)
                    ? Path.GetFileNameWithoutExtension(tex.TEXfileName)
                    : $"{baseFileName}_tex{i}";

                string pngFileName = $"{texName}.png";
                string pngPath = Path.Combine(outputDirectory, pngFileName);

                var tmp = new Bitmap(tex.bitmap);
                tmp.RotateFlip(RotateFlipType.RotateNoneFlipY);
                tmp.Save(pngPath, ImageFormat.Png);
                textureFiles[i] = pngFileName;
            }

            return textureFiles;
        }

        /// <summary>
        /// Creates an Assimp Material with an optional texture reference.
        /// </summary>
        private static Material CreateMaterial(int groupIndex, string? textureFileName = null)
        {
            var material = new Material();
            material.Name = $"Material_{groupIndex}";

            // Use pure white for all material colors so vertex colors are displayed as-is
            material.ColorDiffuse = new System.Numerics.Vector4(1.0f, 1.0f, 1.0f, 1.0f);
            material.ColorAmbient = new System.Numerics.Vector4(1.0f, 1.0f, 1.0f, 1.0f);
            material.ColorEmissive = new System.Numerics.Vector4(0, 0, 0, 1);
            material.ColorSpecular = new System.Numerics.Vector4(0.1f, 0.1f, 0.1f, 1.0f);
            material.Shininess = 1.0f;
            material.Opacity = 1.0f;
            material.ShadingMode = ShadingMode.Gouraud;

            // Add texture reference if provided
            if (!string.IsNullOrEmpty(textureFileName))
            {
                var textureSlot = new TextureSlot
                {
                    FilePath = textureFileName,
                    TextureType = TextureType.Diffuse,
                    TextureIndex = 0,
                    Mapping = TextureMapping.FromUV,
                    UVIndex = 0,
                    BlendFactor = 1.0f,
                    Operation = TextureOperation.Multiply,
                    WrapModeU = TextureWrapMode.Wrap,
                    WrapModeV = TextureWrapMode.Wrap
                };

                material.AddMaterialTexture(textureSlot);
            }

            return material;
        }

        /// <summary>
        /// Converts a FieldSkeleton to an Assimp Scene for export.
        /// </summary>
        /// <param name="skeleton">The FieldSkeleton to convert</param>
        /// <param name="textureFileMap">Optional mapping of texture index to filename (for material texture references)</param>
        /// <param name="bAdjust">Whether to reverse FF7 coordinate adjustments</param>
        /// <returns>An Assimp Scene representing the skeleton hierarchy and models</returns>
        public static Scene ConvertFieldSkeletonToScene(FieldSkeleton skeleton, FieldAnimation animation, Dictionary<int, string>? textureFileMap, bool bAdjust = true)
        {
            var scene = new Scene();
            scene.RootNode = new Node(skeleton.name ?? "FieldSkeleton");

            if (skeleton.bones == null || skeleton.bones.Count == 0)
            {
                // Empty skeleton - add a default material
                scene.Materials.Add(new Material { Name = "Default" });
                return scene;
            }

            // Always create an Armature node to hold the skeleton hierarchy
            // This helps 3D applications recognize this as a rigged model
            var armatureNode = new Node("Armature");
            // When bAdjust flips X and Z, add 180° Z rotation to match coordinate system
            // Note: FF7's Y axis maps to glTF/Blender's Z axis
            if (bAdjust)
            {
                armatureNode.Transform = ToMatrix4x4(Matrix4.CreateRotationZ((float)Math.PI));
            }
            else
            {
                armatureNode.Transform = ToMatrix4x4(Matrix4.Identity);
            }
            scene.RootNode.Children.Add(armatureNode);

            // Dictionary to find nodes by joint name
            var jointNodes = new Dictionary<string, Node>();

            // Dictionary to track bone lengths by joint name (for parent length lookup)
            var boneLengths = new Dictionary<string, double>();
            foreach (var b in skeleton.bones)
            {
                boneLengths[b.joint_i] = b.len;
            }

            // Count unique root joints (joint_f values that don't appear as any bone's joint_i)
            var allJointIs = new HashSet<string>(skeleton.bones.Select(b => b.joint_i));
            var rootJointFs = skeleton.bones.Select(b => b.joint_f).Distinct()
                .Where(jf => !allJointIs.Contains(jf)).ToList();

            // Map all root joint_f names to the armature
            foreach (var jf in rootJointFs)
            {
                jointNodes[jf] = armatureNode;
            }

            // Process each bone
            for (int bi = 0; bi < skeleton.bones.Count; bi++)
            {
                var bone = skeleton.bones[bi];

                // Create node for this bone
                var boneNode = new Node(bone.joint_i);

                // In FF7's rendering, bones are positioned at their PARENT's end (after parent's bone length).
                // Root bones (parent is "root" or armature) have identity transform.
                // Non-root bones translate by their parent's bone length along -Z.
                // When bAdjust is true, mesh vertices have Z negated, so bone translation Z must also be negated.
                Matrix4 transform;
                if (rootJointFs.Contains(bone.joint_f))
                {
                    // Root bone - no parent length to apply
                    transform = Matrix4.Identity;
                }
                else if (boneLengths.TryGetValue(bone.joint_f, out double parentLen))
                {
                    // Child bone - translate by parent's bone length
                    // FF7 uses -Z, but if bAdjust flips mesh Z, we need +Z to match
                    float zDir = bAdjust ? 1 : -1;
                    transform = Matrix4.CreateTranslation(0, 0, (float)(parentLen * zDir));
                }
                else
                {
                    // Fallback - shouldn't happen but use identity if parent not found
                    transform = Matrix4.Identity;
                }
                boneNode.Transform = ToMatrix4x4(transform);

                // Find parent node
                if (jointNodes.TryGetValue(bone.joint_f, out var parentNode))
                {
                    parentNode.Children.Add(boneNode);
                }
                else
                {
                    // Parent not found, attach to armature
                    armatureNode.Children.Add(boneNode);
                }

                // Register this bone's joint_i for children to find
                jointNodes[bone.joint_i] = boneNode;

                // Process resources (models) attached to this bone
                for (int ri = 0; ri < bone.nResources; ri++)
                {
                    var resource = bone.fRSDResources[ri];
                    if (resource.Model.Polys == null || resource.Model.Header.numGroups == 0)
                        continue;

                    // Create a child node for the model with its local transform
                    var modelNode = new Node($"Model_{bi}_{ri}");

                    // Build model transform from reposition, rotation, resize
                    // FF7 rendering order: Translate -> Rotate (quaternion) -> Scale
                    // Matrix multiplication order: T * R * S
                    var modelTransform = Matrix4.Identity;
                    modelTransform *= Matrix4.CreateTranslation(
                        resource.Model.repositionX,
                        resource.Model.repositionY,
                        resource.Model.repositionZ);

                    // Apply rotation quaternion
                    var quat = new Quaternion(
                        (float)resource.Model.rotationQuaternion.X,
                        (float)resource.Model.rotationQuaternion.Y,
                        (float)resource.Model.rotationQuaternion.Z,
                        (float)resource.Model.rotationQuaternion.W);
                    modelTransform *= Matrix4.CreateFromQuaternion(quat);

                    modelTransform *= Matrix4.CreateScale(resource.Model.resizeX, resource.Model.resizeY, resource.Model.resizeZ);

                    modelNode.Transform = ToMatrix4x4(modelTransform);

                    // Convert each group to a mesh
                    for (int g = 0; g < resource.Model.Header.numGroups; g++)
                    {
                        var group = resource.Model.Groups[g];
                        if (group.HiddenQ || group.numVert == 0 || group.numPoly == 0)
                            continue;

                        // Field models use polygon colors (flat shading per face)
                        var mesh = ConvertGroupToMesh(resource.Model, g, bAdjust, usePolygonColors: true);
                        mesh.Name = $"Bone{bi}_Res{ri}_Group{g}";
                        mesh.MaterialIndex = scene.MaterialCount;

                        scene.Meshes.Add(mesh);
                        modelNode.MeshIndices.Add(scene.MeshCount - 1);

                        // Get texture filename for this group if available
                        string? textureFileName = null;
                        if (textureFileMap != null && group.texFlag == 1)
                        {
                            textureFileMap.TryGetValue(group.texID, out textureFileName);
                        }

                        // Create material with optional texture
                        var material = CreateMaterial(g, textureFileName);
                        scene.Materials.Add(material);
                    }

                    boneNode.Children.Add(modelNode);
                }
            }

            // Note: FF7 uses rigid binding - each mesh is 100% attached to its parent bone node.
            // We rely on the node hierarchy (mesh → modelNode → boneNode) rather than skinning
            // bones. Animations drive the bone nodes directly.

            // Ensure at least one material exists
            if (scene.MaterialCount == 0)
            {
                scene.Materials.Add(new Material { Name = "Default" });
            }

            // Add animation data if present
            if (animation.frames != null && animation.nFrames > 0)
            {
                AddFieldAnimationToScene(scene, skeleton, animation, bAdjust);
            }

            return scene;
        }

        /// <summary>
        /// Exports a FieldSkeleton to a file using Assimp.
        /// </summary>
        /// <param name="skeleton">The FieldSkeleton to export</param>
        /// <param name="filePath">Output file path</param>
        /// <param name="bAdjust">Whether to reverse FF7 coordinate adjustments</param>
        /// <param name="exportTextures">Whether to export textures as PNG files</param>
        /// <returns>True if export succeeded</returns>
        public static bool ExportFieldSkeleton(FieldSkeleton skeleton, FieldAnimation animation, string filePath, bool bAdjust)
        {
            string outputDir = Path.GetDirectoryName(filePath) ?? ".";
            string baseName = Path.GetFileNameWithoutExtension(filePath);

            // Export textures and build mapping if requested
            Dictionary<int, string>? textureFileMap = null;
            if (skeleton.textures_pool != null && skeleton.textures_pool.Count > 0)
            {
                textureFileMap = ExportTexturesToPNG(skeleton.textures_pool, outputDir, baseName);
            }

            var scene = ConvertFieldSkeletonToScene(skeleton, animation, textureFileMap, bAdjust);

            string formatId = GetExportFormatFromExtension(Path.GetExtension(filePath));
            if (string.IsNullOrEmpty(formatId))
            {
                throw new NotSupportedException(
                    $"Export format not supported for extension: {Path.GetExtension(filePath)}");
            }

            return _context.ExportFile(scene, filePath, formatId);
        }

        /// <summary>
        /// Converts a BattleSkeleton to an Assimp Scene for export.
        /// </summary>
        /// <param name="skeleton">The BattleSkeleton to convert</param>
        /// <param name="animation">The BattleAnimation to convert</param>
        /// <param name="weaponAnimation">The weapon animation to convert, if it exists</param>
        /// <param name="textureFileMap">Optional mapping of texture index to filename (for material texture references)</param>
        /// <param name="bAdjust">Whether to reverse FF7 coordinate adjustments</param>
        /// <returns>An Assimp Scene representing the skeleton hierarchy and models</returns>
        public static Scene ConvertBattleSkeletonToScene(BattleSkeleton skeleton,
            BattleAnimation animation, BattleAnimation? weaponAnimation,
            Dictionary<int, string>? textureFileMap, bool bAdjust)
        {
            var scene = new Scene();
            scene.RootNode = new Node(skeleton.fileName ?? "BattleSkeleton");

            // Store skeleton type metadata
            // SetFF7Metadata(scene.RootNode, "FF7_SkeletonType", (int)skeleton.skeletonType);
            // SetFF7Metadata(scene.RootNode, "FF7_IsBattleLocation", skeleton.IsBattleLocation ? 1 : 0);

            if (skeleton.bones == null || skeleton.bones.Count == 0)
            {
                scene.Materials.Add(new Material { Name = "Default" });
                return scene;
            }

            // Always create an Armature node to hold the skeleton hierarchy
            // This helps 3D applications recognize this as a rigged model
            var armatureNode = new Node("Armature");
            // When bAdjust flips X and Z, add 180° Z rotation to match coordinate system
            // Note: FF7's Y axis maps to glTF/Blender's Z axis
            if (bAdjust)
            {
                armatureNode.Transform = ToMatrix4x4(Matrix4.CreateRotationZ((float)Math.PI));
            }
            else
            {
                armatureNode.Transform = ToMatrix4x4(Matrix4.Identity);
            }
            scene.RootNode.Children.Add(armatureNode);

            // Create nodes array parallel to bones for parent lookup
            var boneNodes = new Node[skeleton.bones.Count];

            for (int bi = 0; bi < skeleton.bones.Count; bi++)
            {
                var bone = skeleton.bones[bi];

                var boneNode = new Node($"Bone_{bi}");

                // In FF7's rendering, bones are positioned at their PARENT's end (after parent's bone length).
                // Root bones (parentBone == -1) have identity transform.
                // Non-root bones translate by their parent's bone length along +Z.
                // When bAdjust is true, mesh vertices have Z negated, so bone translation Z must also be negated.
                Matrix4 transform;
                if (bone.parentBone >= 0 && bone.parentBone < skeleton.bones.Count)
                {
                    // Child bone - translate by parent's bone length
                    // FF7 uses +Z, but if bAdjust flips mesh Z, we need -Z to match
                    float parentLen = skeleton.bones[bone.parentBone].len;
                    float zDir = bAdjust ? -1 : 1;
                    transform = Matrix4.CreateTranslation(0, 0, parentLen * zDir);
                }
                else
                {
                    // Root bone - no parent length to apply
                    transform = Matrix4.Identity;
                }
                boneNode.Transform = ToMatrix4x4(transform);

                // Store FF7-specific metadata
                // SetFF7Metadata(boneNode, "FF7_BoneLength", (double)bone.len);
                // SetFF7Metadata(boneNode, "FF7_ParentIndex", bone.parentBone);
                // SetFF7Metadata(boneNode, "FF7_ResizeX", bone.resizeX);
                // SetFF7Metadata(boneNode, "FF7_ResizeY", bone.resizeY);
                // SetFF7Metadata(boneNode, "FF7_ResizeZ", bone.resizeZ);
                // SetFF7Metadata(boneNode, "FF7_HasModel", bone.hasModel);

                // Find parent
                if (bone.parentBone >= 0 && bone.parentBone < bi && boneNodes[bone.parentBone] != null)
                {
                    boneNodes[bone.parentBone].Children.Add(boneNode);
                }
                else
                {
                    // Root bone or invalid parent - attach to armature node
                    armatureNode.Children.Add(boneNode);
                }

                boneNodes[bi] = boneNode;

                // Process models attached to this bone
                if (bone.hasModel == 1 && bone.Models != null)
                {
                    for (int mi = 0; mi < bone.nModels; mi++)
                    {
                        var model = bone.Models[mi];
                        if (model.Polys == null || model.Header.numGroups == 0)
                            continue;

                        // Create a child node for the model
                        var modelNode = new Node($"Model_{bi}_{mi}");

                        // Build model transform
                        // FF7 rendering order: Translate -> Rotate (quaternion) -> Scale
                        // Matrix multiplication order: T * R * S
                        var modelTransform = Matrix4.Identity;
                        modelTransform *= Matrix4.CreateTranslation(model.repositionX, model.repositionY, model.repositionZ);
                        modelTransform *= Matrix4.CreateRotationX(MathHelper.DegreesToRadians((float)model.rotateAlpha));
                        modelTransform *= Matrix4.CreateRotationY(MathHelper.DegreesToRadians((float)model.rotateBeta));
                        modelTransform *= Matrix4.CreateRotationZ(MathHelper.DegreesToRadians((float)model.rotateGamma));
                        modelTransform *= Matrix4.CreateScale(model.resizeX, model.resizeY, model.resizeZ);

                        modelNode.Transform = ToMatrix4x4(modelTransform);

                        // Store model metadata
                        // SetFF7Metadata(modelNode, "FF7_ModelReposition",
                        //    $"{model.repositionX},{model.repositionY},{model.repositionZ}");
                        // SetFF7Metadata(modelNode, "FF7_ModelRotation",
                        //    $"{model.rotateAlpha},{model.rotateBeta},{model.rotateGamma}");
                        // SetFF7Metadata(modelNode, "FF7_ModelResize",
                        //    $"{model.resizeX},{model.resizeY},{model.resizeZ}");

                        // Convert groups to meshes
                        for (int g = 0; g < model.Header.numGroups; g++)
                        {
                            var group = model.Groups[g];
                            if (group.HiddenQ || group.numVert == 0 || group.numPoly == 0)
                                continue;

                            // Battle models use vertex colors (smooth shading interpolated across vertices)
                            var mesh = ConvertGroupToMesh(model, g, bAdjust, usePolygonColors: false);
                            mesh.Name = $"Bone{bi}_Model{mi}_Group{g}";
                            mesh.MaterialIndex = scene.MaterialCount;

                            scene.Meshes.Add(mesh);
                            modelNode.MeshIndices.Add(scene.MeshCount - 1);

                            // Get texture filename for this group if available
                            string? textureFileName = null;
                            if (textureFileMap != null && group.texFlag == 1)
                            {
                                textureFileMap.TryGetValue(group.texID, out textureFileName);
                            }

                            var material = CreateMaterial(g, textureFileName);
                            scene.Materials.Add(material);
                        }

                        boneNode.Children.Add(modelNode);
                    }
                }
            }

            // Note: FF7 uses rigid binding - each mesh is 100% attached to its parent bone node.
            // We rely on the node hierarchy (mesh → modelNode → boneNode) rather than skinning
            // bones. Animations drive the bone nodes directly.

            // Process weapon models
            // Structure: Armature -> WeaponArmature -> Weapons -> Weapon_0, Weapon_1, etc.
            // WeaponArmature is a child of Armature so it inherits the skeleton's position/rotation.
            // The weapon animation then provides the relative offset from the skeleton.
            if (skeleton.wpModels != null && skeleton.wpModels.Count > 0)
            {
                var weaponArmatureNode = new Node("WeaponArmature");
                // WeaponArmature inherits skeleton transforms, so leave static transform as identity.
                // Animation provides relative position offset.

                armatureNode.Children.Add(weaponArmatureNode);

                var weaponsNode = new Node("Weapons");
                weaponArmatureNode.Children.Add(weaponsNode);

                for (int wi = 0; wi < skeleton.wpModels.Count; wi++)
                {
                    var wpModel = skeleton.wpModels[wi];
                    if (wpModel.Polys == null || wpModel.Header.numGroups == 0)
                        continue;

                    var weaponNode = new Node($"Weapon_{wi}");

                    // Build weapon transform
                    // Apply bAdjust coordinate transformation to weapon reposition (same as skeleton)
                    float wpReposX = bAdjust ? -wpModel.repositionX : wpModel.repositionX;
                    float wpReposZ = bAdjust ? -wpModel.repositionZ : wpModel.repositionZ;
                    float wpRotAlpha = bAdjust ? -wpModel.rotateAlpha : wpModel.rotateAlpha;
                    float wpRotGamma = bAdjust ? -wpModel.rotateGamma : wpModel.rotateGamma;

                    // FF7 rendering order: Translate -> RotX -> RotY -> RotZ -> Scale (applied to vertices)
                    // Matrix multiplication order: T * Rx * Ry * Rz * S
                    var wpTransform = Matrix4.Identity;
                    wpTransform *= Matrix4.CreateTranslation(wpReposX, wpModel.repositionY, wpReposZ);
                    wpTransform *= Matrix4.CreateRotationX(MathHelper.DegreesToRadians((float)wpRotAlpha));
                    wpTransform *= Matrix4.CreateRotationY(MathHelper.DegreesToRadians((float)wpModel.rotateBeta));
                    wpTransform *= Matrix4.CreateRotationZ(MathHelper.DegreesToRadians((float)wpRotGamma));
                    wpTransform *= Matrix4.CreateScale(wpModel.resizeX, wpModel.resizeY, wpModel.resizeZ);

                    weaponNode.Transform = ToMatrix4x4(wpTransform);

                    for (int g = 0; g < wpModel.Header.numGroups; g++)
                    {
                        var group = wpModel.Groups[g];
                        if (group.HiddenQ || group.numVert == 0 || group.numPoly == 0)
                            continue;

                        // Weapon models use vertex colors (same as battle models)
                        var mesh = ConvertGroupToMesh(wpModel, g, bAdjust, usePolygonColors: false);
                        mesh.Name = $"Weapon{wi}_Group{g}";
                        mesh.MaterialIndex = scene.MaterialCount;

                        scene.Meshes.Add(mesh);
                        weaponNode.MeshIndices.Add(scene.MeshCount - 1);

                        // Get texture filename for this group if available
                        string? textureFileName = null;
                        if (textureFileMap != null && group.texFlag == 1)
                        {
                            textureFileMap.TryGetValue(group.texID, out textureFileName);
                        }

                        var material = CreateMaterial(g, textureFileName);
                        scene.Materials.Add(material);
                    }

                    weaponsNode.Children.Add(weaponNode);
                }
            }

            // Ensure at least one material exists
            if (scene.MaterialCount == 0)
            {
                scene.Materials.Add(new Material { Name = "Default" });
            }

            // Add animation data if present (includes both skeleton and weapon animation in one clip)
            if (animation.frames != null && animation.numFramesShort > 0)
            {
                AddBattleAnimationToScene(scene, skeleton, animation, weaponAnimation, bAdjust);
            }

            return scene;
        }

        /// <summary>
        /// Converts a BattleLocation to an Assimp Scene for export.
        /// Battle locations are different from battle skeletons: they have flat hierarchy
        /// where each piece is independent, with no bone-based animation.
        /// </summary>
        /// <param name="skeleton">The BattleSkeleton representing the battle location</param>
        /// <param name="textureFileMap">Optional mapping of texture index to filename</param>
        /// <param name="bAdjust">Whether to reverse FF7 coordinate adjustments</param>
        /// <returns>An Assimp Scene representing the battle location</returns>
        public static Scene ConvertBattleLocationToScene(BattleSkeleton skeleton,
            Dictionary<int, string>? textureFileMap, bool bAdjust)
        {
            var scene = new Scene();
            scene.RootNode = new Node(skeleton.fileName ?? "BattleLocation");

            if (skeleton.bones == null || skeleton.bones.Count == 0)
            {
                scene.Materials.Add(new Material { Name = "Default" });
                return scene;
            }

            // Create an Armature node to hold all pieces
            var armatureNode = new Node("Armature");
            // When bAdjust flips X and Z, add 180° Z rotation to match coordinate system
            if (bAdjust)
            {
                armatureNode.Transform = ToMatrix4x4(Matrix4.CreateRotationZ((float)Math.PI));
            }
            else
            {
                armatureNode.Transform = ToMatrix4x4(Matrix4.Identity);
            }
            scene.RootNode.Children.Add(armatureNode);

            // Battle locations have flat hierarchy - each piece is independent
            // Each "bone" is really just a piece with its own transform
            for (int bi = 0; bi < skeleton.bones.Count; bi++)
            {
                var bone = skeleton.bones[bi];

                // Create a node for this piece (direct child of Armature - no hierarchy)
                var pieceNode = new Node($"Piece_{bi}");

                // Battle location pieces have identity transform at the piece level
                // The actual positioning comes from the model's reposition/rotation properties
                pieceNode.Transform = ToMatrix4x4(Matrix4.Identity);

                armatureNode.Children.Add(pieceNode);

                // Process models attached to this piece
                if (bone.hasModel == 1 && bone.Models != null)
                {
                    for (int mi = 0; mi < bone.nModels; mi++)
                    {
                        var model = bone.Models[mi];
                        if (model.Polys == null || model.Header.numGroups == 0)
                            continue;

                        // Create a child node for the model
                        var modelNode = new Node($"Model_{bi}_{mi}");

                        // Build model transform following FF7 rendering order:
                        // 1. Bone resize (bBone.resizeX/Y/Z)
                        // 2. Model translate (repositionX/Y/Z)
                        // 3. Model rotate X->Y->Z (rotateAlpha->rotateBeta->rotateGamma)
                        // 4. Model resize (resizeX/Y/Z)
                        // Matrix multiplication order: BoneScale * Translate * RotX * RotY * RotZ * ModelScale

                        // Apply bAdjust coordinate transformation
                        float reposX = bAdjust ? -model.repositionX : model.repositionX;
                        float reposZ = bAdjust ? -model.repositionZ : model.repositionZ;
                        float rotAlpha = bAdjust ? -model.rotateAlpha : model.rotateAlpha;
                        float rotGamma = bAdjust ? -model.rotateGamma : model.rotateGamma;

                        var modelTransform = Matrix4.Identity;
                        modelTransform *= Matrix4.CreateScale(bone.resizeX, bone.resizeY, bone.resizeZ);
                        modelTransform *= Matrix4.CreateTranslation(reposX, model.repositionY, reposZ);
                        modelTransform *= Matrix4.CreateRotationX(MathHelper.DegreesToRadians((float)rotAlpha));
                        modelTransform *= Matrix4.CreateRotationY(MathHelper.DegreesToRadians((float)model.rotateBeta));
                        modelTransform *= Matrix4.CreateRotationZ(MathHelper.DegreesToRadians((float)rotGamma));
                        modelTransform *= Matrix4.CreateScale(model.resizeX, model.resizeY, model.resizeZ);

                        modelNode.Transform = ToMatrix4x4(modelTransform);

                        // Convert groups to meshes
                        for (int g = 0; g < model.Header.numGroups; g++)
                        {
                            var group = model.Groups[g];
                            if (group.HiddenQ || group.numVert == 0 || group.numPoly == 0)
                                continue;

                            // Battle locations use vertex colors
                            var mesh = ConvertGroupToMesh(model, g, bAdjust, usePolygonColors: false);
                            mesh.Name = $"Piece{bi}_Model{mi}_Group{g}";
                            mesh.MaterialIndex = scene.MaterialCount;

                            scene.Meshes.Add(mesh);
                            modelNode.MeshIndices.Add(scene.MeshCount - 1);

                            // Get texture filename for this group if available
                            string? textureFileName = null;
                            if (textureFileMap != null && group.texFlag == 1)
                            {
                                textureFileMap.TryGetValue(group.texID, out textureFileName);
                            }

                            var material = CreateMaterial(g, textureFileName);
                            scene.Materials.Add(material);
                        }

                        pieceNode.Children.Add(modelNode);
                    }
                }
            }

            // Ensure at least one material exists
            if (scene.MaterialCount == 0)
            {
                scene.Materials.Add(new Material { Name = "Default" });
            }

            // Battle locations don't have meaningful animation - skip animation export

            return scene;
        }

        /// <summary>
        /// Exports a BattleSkeleton to a file using Assimp.
        /// </summary>
        /// <param name="skeleton">The BattleSkeleton to export</param>
        /// <param name="filePath">Output file path</param>
        /// <param name="bAdjust">Whether to reverse FF7 coordinate adjustments</param>
        /// <param name="includeWeapons">Whether to include weapon models</param>
        /// <returns>True if export succeeded</returns>
        public static bool ExportBattleSkeleton(BattleSkeleton skeleton, BattleAnimation animation,
            BattleAnimation? weaponAnimation, string filePath, bool bAdjust)
        {
            string outputDir = Path.GetDirectoryName(filePath) ?? ".";
            string baseName = Path.GetFileNameWithoutExtension(filePath);

            // Export textures if they exist
            Dictionary<int, string>? textureFileMap = null;
            if (skeleton.textures != null && skeleton.textures.Count > 0)
            {
                textureFileMap = ExportTexturesToPNG(skeleton.textures, outputDir, baseName);
            }

            // Use appropriate conversion based on whether this is a battle location or battle skeleton
            Scene scene;
            if (skeleton.IsBattleLocation)
            {
                // Battle locations have flat hierarchy with no animation
                scene = ConvertBattleLocationToScene(skeleton, textureFileMap, bAdjust);
            }
            else
            {
                // Regular battle skeletons have bone hierarchy and animation
                scene = ConvertBattleSkeletonToScene(skeleton, animation, weaponAnimation,
                    textureFileMap, bAdjust);
            }

            string formatId = GetExportFormatFromExtension(Path.GetExtension(filePath));
            if (string.IsNullOrEmpty(formatId))
            {
                throw new NotSupportedException(
                    $"Export format not supported for extension: {Path.GetExtension(filePath)}");
            }

            return _context.ExportFile(scene, filePath, formatId);
        }

        /// <summary>
        /// Exports a skeleton of the specified type
        /// </summary>
        /// 
        public static void ExportSkeleton(UnifiedSkeleton skeleton, UnifiedAnimation animation,
            UnifiedAnimation? weaponAnimation, ModelType modelType, string filePath, bool bAdjust)
        {
            if (modelType == ModelType.HRCSkeleton)
            {
                var fSkeleton = skeleton.ToFieldSkeleton();
                var fAnimation = animation.ToFieldAnimation();
                ExportFieldSkeleton(fSkeleton, fAnimation, filePath, bAdjust);
            }
            else
            {
                var bSkeleton = skeleton.ToBattleSkeleton();
                var bAnimation = animation.ToBattleAnimation();
                var wpAnimation = weaponAnimation?.ToBattleAnimation();
                ExportBattleSkeleton(bSkeleton, bAnimation, wpAnimation, filePath, bAdjust);
            }
        }

        /// <summary>
        /// Adds field animation data to an Assimp scene.
        /// </summary>
        /// <param name="scene">The scene to add animation to</param>
        /// <param name="skeleton">The field skeleton for bone hierarchy info</param>
        /// <param name="animation">The field animation data</param>
        /// <param name="bAdjust">Whether coordinate adjustments are applied (affects Z direction)</param>
        /// <param name="ticksPerSecond">Animation playback rate (default 30 fps)</param>
        private static void AddFieldAnimationToScene(Scene scene, FieldSkeleton skeleton, FieldAnimation animation, bool bAdjust, double ticksPerSecond = 30.0)
        {
            if (animation.frames == null || animation.nFrames == 0)
                return;

            var anim = new Animation
            {
                Name = animation.strFieldAnimationFile ?? "FieldAnimation",
                DurationInTicks = animation.nFrames,
                TicksPerSecond = ticksPerSecond
            };

            // Z direction for bone translations: FF7 uses -Z, but bAdjust flips mesh Z so we need +Z
            float zDir = bAdjust ? 1 : -1;

            // Create channel for the Armature node (root transforms)
            var armatureChannel = new NodeAnimationChannel
            {
                NodeName = "Armature"
            };

            for (int fi = 0; fi < animation.nFrames; fi++)
            {
                var frame = animation.frames[fi];

                // Root translation key - Y is negated in FF7 rendering (see ModelDrawing.cs line 425)
                // X and Z are negated if bAdjust is true (coordinate system conversion)
                // When bAdjust, the 180° Z rotation flips Y, so don't pre-negate it
                float xDir = bAdjust ? -1 : 1;
                float yDir = bAdjust ? 1 : -1;
                armatureChannel.PositionKeys.Add(new VectorKey(fi,
                    new System.Numerics.Vector3(
                        frame.rootTranslationX * xDir,
                        frame.rootTranslationY * yDir,
                        frame.rootTranslationZ * zDir)));

                // Root rotation key
                // For 180° Y rotation (X and Z negated), negate X and Z rotations for all rotations
                float rootAlpha = frame.rootRotationAlpha;
                float rootBeta = frame.rootRotationBeta;
                float rootGamma = frame.rootRotationGamma;
                if (bAdjust)
                {
                    rootAlpha = -rootAlpha;  // Negate X rotation
                    rootGamma = -rootGamma;  // Negate Z rotation
                    // Y rotation (beta) unchanged
                }
                var rootQuat = EulerToQuaternionFF7(rootAlpha, rootBeta, rootGamma, false);

                // Apply 180° Z rotation to match coordinate system (animation overrides static transform)
                // Note: FF7's Y axis maps to glTF/Blender's Z axis
                if (bAdjust)
                {
                    var rot180Z = System.Numerics.Quaternion.CreateFromAxisAngle(
                        System.Numerics.Vector3.UnitZ, (float)Math.PI);
                    rootQuat = System.Numerics.Quaternion.Concatenate(rootQuat, rot180Z);
                }
                armatureChannel.RotationKeys.Add(new QuaternionKey(fi, rootQuat));

                // Scaling key (constant)
                armatureChannel.ScalingKeys.Add(new VectorKey(fi,
                    new System.Numerics.Vector3(1, 1, 1)));
            }

            anim.NodeAnimationChannels.Add(armatureChannel);

            // Build bone length lookup and identify root joints (same as in ConvertFieldSkeletonToScene)
            var boneLengths = new Dictionary<string, double>();
            foreach (var b in skeleton.bones)
            {
                boneLengths[b.joint_i] = b.len;
            }
            var allJointIs = new HashSet<string>(skeleton.bones.Select(b => b.joint_i));
            var rootJointFs = skeleton.bones.Select(b => b.joint_f).Distinct()
                .Where(jf => !allJointIs.Contains(jf)).ToHashSet();

            // Create channels for each bone
            for (int bi = 0; bi < skeleton.bones.Count && bi < animation.nBones; bi++)
            {
                var bone = skeleton.bones[bi];
                var channel = new NodeAnimationChannel
                {
                    NodeName = bone.joint_i
                };

                // Check if this is a root-level bone (parent is "root" or similar non-bone)
                bool isRootBone = rootJointFs.Contains(bone.joint_f);

                // Calculate bone position (must match static transform in ConvertFieldSkeletonToScene)
                // Root bones have identity, child bones translate by parent's length
                float boneZ = 0;
                if (!isRootBone && boneLengths.TryGetValue(bone.joint_f, out double parentLen))
                {
                    boneZ = (float)(parentLen * zDir);
                }

                for (int fi = 0; fi < animation.nFrames; fi++)
                {
                    var frame = animation.frames[fi];

                    // Get bone rotation for this frame
                    if (bi < frame.rotations.Count)
                    {
                        var rot = frame.rotations[bi];
                        float alpha = rot.alpha;
                        float beta = rot.beta;
                        float gamma = rot.gamma;

                        // For 180° Y rotation (X and Z negated), negate X and Z rotations for all bones
                        if (bAdjust)
                        {
                            alpha = -alpha;  // Negate X rotation
                            gamma = -gamma;  // Negate Z rotation
                            // Y rotation (beta) unchanged
                        }
                        var quat = EulerToQuaternionFF7(alpha, beta, gamma, false);
                        channel.RotationKeys.Add(new QuaternionKey(fi, quat));
                    }
                    else
                    {
                        // Default identity rotation
                        channel.RotationKeys.Add(new QuaternionKey(fi,
                            System.Numerics.Quaternion.Identity));
                    }

                    // Position key - must match static transform (parent bone length along -Z)
                    // Animation keyframes replace static transform, so we need the full position
                    channel.PositionKeys.Add(new VectorKey(fi,
                        new System.Numerics.Vector3(0, 0, boneZ)));

                    // Scale key (constant)
                    channel.ScalingKeys.Add(new VectorKey(fi,
                        new System.Numerics.Vector3(1, 1, 1)));
                }

                anim.NodeAnimationChannels.Add(channel);
            }

            scene.Animations.Add(anim);
        }

        /// <summary>
        /// Adds battle animation data to an Assimp scene.
        /// </summary>
        /// <param name="scene">The scene to add animation to</param>
        /// <param name="skeleton">The battle skeleton for bone hierarchy info</param>
        /// <param name="animation">The battle animation data</param>
        /// <param name="bAdjust">Whether coordinate adjustments are applied (affects Z direction)</param>
        /// <param name="ticksPerSecond">Animation playback rate (default 30 fps)</param>
        private static void AddBattleAnimationToScene(Scene scene, BattleSkeleton skeleton, BattleAnimation animation, BattleAnimation? weaponAnimation, bool bAdjust, double ticksPerSecond = 30.0)
        {
            if (animation.frames == null || animation.numFramesShort == 0)
                return;

            var anim = new Animation
            {
                Name = "BattleAnimation",
                DurationInTicks = animation.numFramesShort,
                TicksPerSecond = ticksPerSecond
            };

            // Z direction for bone translations: FF7 uses +Z, but bAdjust flips mesh Z so we need -Z
            float zDir = bAdjust ? -1 : 1;

            // Determine bone index offset for animation data
            // If skeleton has more than 1 bone, animation bone[0] is root rotation, bone[1+] are actual bones
            // If skeleton has only 1 bone, animation bone[0] is both root and the only bone
            int itmpbones = skeleton.nBones > 1 ? 1 : 0;

            // Create channel for the Armature node (root position from frame start + rotation from bone[0])
            var armatureChannel = new NodeAnimationChannel
            {
                NodeName = "Armature"
            };

            for (int fi = 0; fi < animation.numFramesShort && fi < animation.frames.Count; fi++)
            {
                var frame = animation.frames[fi];

                // Root translation from frame start position
                // X and Z negated if bAdjust, Y not negated when bAdjust (180° Z rotation handles it)
                float xDir = bAdjust ? -1 : 1;
                armatureChannel.PositionKeys.Add(new VectorKey(fi,
                    new System.Numerics.Vector3(
                        frame.startX * xDir,
                        frame.startY,
                        frame.startZ * zDir)));

                // Root rotation from bone[0]
                System.Numerics.Quaternion rootQuat;
                if (frame.bones != null && frame.bones.Count > 0)
                {
                    var rootBone = frame.bones[0];
                    float alpha = rootBone.alpha;
                    float beta = rootBone.beta;
                    float gamma = rootBone.gamma;

                    // For bAdjust, negate X and Z rotations
                    if (bAdjust)
                    {
                        alpha = -alpha;
                        gamma = -gamma;
                    }
                    rootQuat = EulerToQuaternionFF7(alpha, beta, gamma, false);
                }
                else
                {
                    rootQuat = System.Numerics.Quaternion.Identity;
                }

                // Apply 180° Z rotation to match coordinate system (animation overrides static transform)
                // Note: FF7's Y axis maps to glTF/Blender's Z axis
                if (bAdjust)
                {
                    var rot180Z = System.Numerics.Quaternion.CreateFromAxisAngle(
                        System.Numerics.Vector3.UnitZ, (float)Math.PI);
                    rootQuat = System.Numerics.Quaternion.Concatenate(rootQuat, rot180Z);
                }
                armatureChannel.RotationKeys.Add(new QuaternionKey(fi, rootQuat));

                // Scaling key (constant)
                armatureChannel.ScalingKeys.Add(new VectorKey(fi,
                    new System.Numerics.Vector3(1, 1, 1)));
            }

            anim.NodeAnimationChannels.Add(armatureChannel);

            // Create channels for each bone
            for (int bi = 0; bi < skeleton.bones.Count; bi++)
            {
                var bone = skeleton.bones[bi];
                var channel = new NodeAnimationChannel
                {
                    NodeName = $"Bone_{bi}"
                };

                // Calculate bone position (must match static transform in ConvertBattleSkeletonToScene)
                // Root bones have identity, child bones translate by parent's length
                float boneZ = 0;
                if (bone.parentBone >= 0 && bone.parentBone < skeleton.bones.Count)
                {
                    boneZ = skeleton.bones[bone.parentBone].len * zDir;
                }

                for (int fi = 0; fi < animation.numFramesShort && fi < animation.frames.Count; fi++)
                {
                    var frame = animation.frames[fi];

                    // Animation bone index: bone[bi + itmpbones] in animation data
                    // (matches ModelDrawing.cs DrawBattleSkeleton logic)
                    int animBoneIndex = bi + itmpbones;

                    if (frame.bones != null && animBoneIndex < frame.bones.Count)
                    {
                        var frameBone = frame.bones[animBoneIndex];
                        float alpha = frameBone.alpha;
                        float beta = frameBone.beta;
                        float gamma = frameBone.gamma;

                        // For bAdjust, negate X and Z rotations
                        if (bAdjust)
                        {
                            alpha = -alpha;
                            gamma = -gamma;
                        }
                        var quat = EulerToQuaternionFF7(alpha, beta, gamma, false);
                        channel.RotationKeys.Add(new QuaternionKey(fi, quat));
                    }
                    else
                    {
                        // Default identity rotation
                        channel.RotationKeys.Add(new QuaternionKey(fi,
                            System.Numerics.Quaternion.Identity));
                    }

                    // Position key - must match static transform (parent bone length along +Z)
                    // Animation keyframes replace static transform, so we need the full position
                    channel.PositionKeys.Add(new VectorKey(fi,
                        new System.Numerics.Vector3(0, 0, boneZ)));

                    // Scale key (constant)
                    channel.ScalingKeys.Add(new VectorKey(fi,
                        new System.Numerics.Vector3(1, 1, 1)));
                }

                anim.NodeAnimationChannels.Add(channel);
            }

            // Add weapon animation channel to the same animation (if present)
            // This ensures skeleton and weapon animate together in sync
            if (weaponAnimation != null && weaponAnimation.Value.frames != null && weaponAnimation.Value.numFramesShort > 0)
            {
                AddWeaponChannelToAnimation(anim, animation, weaponAnimation.Value, bAdjust);
            }

            scene.Animations.Add(anim);
        }

        /// <summary>
        /// Adds weapon animation channel to an existing Animation object.
        /// The weapon animation controls the WeaponArmature node's position and rotation.
        /// Weapon position is computed relative to the skeleton's root position.
        /// </summary>
        private static void AddWeaponChannelToAnimation(Animation anim, BattleAnimation skeletonAnimation, BattleAnimation weaponAnimation, bool bAdjust)
        {
            // Create channel for the WeaponArmature node
            var weaponArmatureChannel = new NodeAnimationChannel
            {
                NodeName = "WeaponArmature"
            };

            // Use the minimum frame count between skeleton and weapon animations
            int frameCount = Math.Min(
                Math.Min(weaponAnimation.numFramesShort, weaponAnimation.frames.Count),
                Math.Min(skeletonAnimation.numFramesShort, skeletonAnimation.frames.Count));

            for (int fi = 0; fi < frameCount; fi++)
            {
                var weaponFrame = weaponAnimation.frames[fi];
                var skelFrame = skeletonAnimation.frames[fi];

                // Compute weapon position RELATIVE to skeleton position
                // This keeps the weapon attached to the character as it moves
                float relX = weaponFrame.startX - skelFrame.startX;
                float relY = weaponFrame.startY - skelFrame.startY;
                float relZ = weaponFrame.startZ - skelFrame.startZ;

                // WeaponArmature is a child of Armature, which has 180° Z rotation.
                // We need to compensate by negating X and Z in the relative position.
                float xDir = bAdjust ? -1 : 1;
                float zDir = bAdjust ? -1 : 1;
                weaponArmatureChannel.PositionKeys.Add(new VectorKey(fi,
                    new System.Numerics.Vector3(
                        relX * xDir,
                        relY,
                        relZ * zDir)));

                // Weapon rotation from frame bone data
                System.Numerics.Quaternion weaponQuat;
                if (weaponFrame.bones != null && weaponFrame.bones.Count > 0)
                {
                    var weaponBone = weaponFrame.bones[0];
                    float alpha = weaponBone.alpha;
                    float beta = weaponBone.beta;
                    float gamma = weaponBone.gamma;

                    // For bAdjust, negate X and Z rotations (same as skeleton)
                    if (bAdjust)
                    {
                        alpha = -alpha;
                        gamma = -gamma;
                    }
                    weaponQuat = EulerToQuaternionFF7(alpha, beta, gamma, false);

                    // Note: WeaponArmature is a child of Armature which already has 180° Z rotation,
                    // so we don't apply it again here (the weapon inherits it from parent)
                }
                else
                {
                    weaponQuat = System.Numerics.Quaternion.Identity;
                }
                weaponArmatureChannel.RotationKeys.Add(new QuaternionKey(fi, weaponQuat));

                // Scaling key (constant)
                weaponArmatureChannel.ScalingKeys.Add(new VectorKey(fi,
                    new System.Numerics.Vector3(1, 1, 1)));
            }

            // Add to existing animation (alongside skeleton channels)
            anim.NodeAnimationChannels.Add(weaponArmatureChannel);
        }

        #endregion
    }
}
