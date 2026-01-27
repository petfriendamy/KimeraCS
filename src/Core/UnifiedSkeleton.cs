using KimeraCS.Rendering;
using OpenTK.Graphics.OpenGL.Compatibility;
using OpenTK.Mathematics;

namespace KimeraCS.Core
{
    using static FF7BattleAnimationsPack;
    using static FF7BattleSkeleton;
    using static FF7FieldAnimation;
    using static FF7FieldRSDResource;
    using static FF7FieldSkeleton;
    using static FF7PModel;
    using static FF7TEXTexture;
    using static GLRenderer;
    using static KimeraCS.Core.FF7BattleAnimation;
    using static Utils;

    /// <summary>
    /// Bone direction along Z axis - field uses negative Z, battle uses positive Z.
    /// </summary>
    public enum BoneDirection
    {
        NegativeZ,  // Field skeletons
        PositiveZ   // Battle skeletons
    }

    /// <summary>
    /// Source format of the skeleton.
    /// </summary>
    public enum SkeletonSourceType
    {
        Field,   // .HRC files
        Battle,  // ??AA files
        Magic    // .D files
    }

    /// <summary>
    /// Represents a model attached to a bone.
    /// </summary>
    public class UnifiedBoneModel
    {
        /// <summary>
        /// The PModel containing geometry.
        /// </summary>
        public PModel Model;

        /// <summary>
        /// Textures used by this model.
        /// </summary>
        public List<TEX> Textures = [];

        public int TextureCount => Textures.Count;

        /// <summary>
        /// Resource file name (RSD for field, model filename for battle).
        /// </summary>
        public string ResourceFile;

        public UnifiedBoneModel(FieldRSDResource rsd)
        {
            Model = rsd.Model;
            ResourceFile = rsd.res_file;

            // Copy textures from RSD resource
            if (rsd.textures != null)
            {
                Textures = rsd.textures;
            }
        }

        public UnifiedBoneModel(PModel model)
        {
            Model = model;
            ResourceFile = string.Empty;
        }

        public UnifiedBoneModel(UnifiedBoneModel other)
        {
            Model = other.Model;
            ResourceFile = other.ResourceFile;

            foreach (var t in other.Textures)
            {
                Textures.Add(t);
            }
        }

        /// <summary>
        /// Converts this unified bone model back to a field RSD resource.
        /// </summary>
        public FieldRSDResource ToRSDResource()
        {
            return new FieldRSDResource
            {
                ID = "@RSD940102",
                res_file = ResourceFile,
                Model = Model,
                numTextures = Textures?.Count ?? 0,
                textures = Textures ?? []
            };
        }
    }

    /// <summary>
    /// Unified bone representation that abstracts field and battle bone differences.
    /// </summary>
    public class UnifiedBone
    {
        /// <summary>
        /// Index of this bone in the skeleton's bone list.
        /// </summary>
        public int Index;

        /// <summary>
        /// Bone name. From joint_i (field) or generated "bone_N" (battle).
        /// </summary>
        public string Name;

        /// <summary>
        /// Index of parent bone, or -1 for root bones.
        /// </summary>
        public int ParentIndex = -1;

        /// <summary>
        /// Parent bone name. From joint_f (field) or looked up from parent index (battle).
        /// </summary>
        public string ParentName;

        /// <summary>
        /// Indices of child bones.
        /// </summary>
        public List<int> ChildIndices = [];

        /// <summary>
        /// Bone length along Z axis (direction determined by BoneDirection).
        /// </summary>
        public float Length;

        /// <summary>
        /// Bone scale factors.
        /// </summary>
        public Vector3 Scale = Vector3.One;

        /// <summary>
        /// Models attached to this bone.
        /// </summary>
        public List<UnifiedBoneModel> Models = [];

        /// <summary>
        /// Whether this bone has any model geometry.
        /// </summary>
        public bool HasModel => Models != null && Models.Count > 0;

        public UnifiedBone(FieldBone fb, int index, int parentIndex)
        {
            Index = index;
            ParentIndex = parentIndex;
            Name = fb.joint_i;
            ParentName = fb.joint_f;
            Length = (float)fb.len;
            Scale = new Vector3(fb.resizeX, fb.resizeY, fb.resizeZ);

            // Convert RSD resources to unified models
            if (fb.fRSDResources != null)
            {
                foreach (var rsd in fb.fRSDResources)
                {
                    Models.Add(new UnifiedBoneModel(rsd));
                }
            }
        }

        public UnifiedBone(BattleBone bb, int index)
        {
            Index = index;
            Name = $"bone_{index}";
            ParentIndex = bb.parentBone;
            ParentName = bb.parentBone >= 0 ? $"bone_{bb.parentBone}" : string.Empty;
            Length = bb.len;
            Scale = new Vector3(bb.resizeX, bb.resizeY, bb.resizeZ);

            // Convert models
            if (bb.hasModel != 0 && bb.Models != null)
            {
                foreach (var pm in bb.Models)
                {
                    Models.Add(new UnifiedBoneModel(pm));
                }
            }
        }

        public UnifiedBone(UnifiedBone other)
        {
            Index = other.Index;
            Name = other.Name;
            ParentIndex = other.ParentIndex;
            ParentName = other.ParentName;
            Length = other.Length;
            Scale = other.Scale;

            foreach (var m in other.Models)
            {
                Models.Add(new UnifiedBoneModel(m));
            }
        }

        /// <summary>
        /// Computes the bounding box for a single unified bone based on its attached models.
        /// </summary>
        public void ComputeBoundingBox(ref Vector3 p_min, ref Vector3 p_max)
        {
            if (!HasModel)
            {
                // No models attached - return zero bounding box
                p_min = Vector3.Zero;
                p_max = Vector3.Zero;
                return;
            }

            p_min.X = float.PositiveInfinity;
            p_min.Y = float.PositiveInfinity;
            p_min.Z = float.PositiveInfinity;
            p_max.X = float.NegativeInfinity;
            p_max.Y = float.NegativeInfinity;
            p_max.Z = float.NegativeInfinity;

            Vector3 p_min_part = new();
            Vector3 p_max_part = new();

            foreach (var model in Models)
            {
                if (model.Model.Verts == null || model.Model.Verts.Length == 0)
                    continue;

                ComputePModelBoundingBox(model.Model, ref p_min_part, ref p_max_part);

                if (p_max.X < p_max_part.X) p_max.X = p_max_part.X;
                if (p_max.Y < p_max_part.Y) p_max.Y = p_max_part.Y;
                if (p_max.Z < p_max_part.Z) p_max.Z = p_max_part.Z;

                if (p_min.X > p_min_part.X) p_min.X = p_min_part.X;
                if (p_min.Y > p_min_part.Y) p_min.Y = p_min_part.Y;
                if (p_min.Z > p_min_part.Z) p_min.Z = p_min_part.Z;
            }

            // If we ended up with no valid geometry, reset to zero
            if (float.IsInfinity(p_min.X))
            {
                p_min = Vector3.Zero;
                p_max = Vector3.Zero;
            }
        }

        /// <summary>
        /// Converts this unified bone back to a field bone structure.
        /// </summary>
        public FieldBone ToFieldBone()
        {
            var fieldBone = new FieldBone
            {
                joint_i = Name,
                joint_f = ParentName ?? string.Empty,
                len = Length,
                resizeX = Scale.X,
                resizeY = Scale.Y,
                resizeZ = Scale.Z,
                nResources = Models?.Count ?? 0,
                fRSDResources = []
            };

            if (Models != null)
            {
                foreach (var model in Models)
                {
                    fieldBone.fRSDResources.Add(model.ToRSDResource());
                }
            }

            return fieldBone;
        }

        /// <summary>
        /// Converts this unified bone back to a battle bone structure.
        /// </summary>
        public BattleBone ToBattleBone()
        {
            var battleBone = new BattleBone
            {
                parentBone = ParentIndex,
                len = Length,
                hasModel = (Models != null && Models.Count > 0) ? 1 : 0,
                nModels = Models?.Count ?? 0,
                resizeX = Scale.X,
                resizeY = Scale.Y,
                resizeZ = Scale.Z,
                Models = []
            };

            if (Models != null)
            {
                foreach (var model in Models)
                {
                    battleBone.Models.Add(model.Model);
                }
            }

            return battleBone;
        }

        /// <summary>
        /// Computes the diameter (diagonal of bounding box) for this bone's models.
        /// Used for field skeleton diameter calculation.
        /// </summary>
        public float ComputeDiameter()
        {
            if (!HasModel)
                return 0;

            Vector3 p_min = new();
            Vector3 p_max = new();

            p_min.X = float.PositiveInfinity;
            p_min.Y = float.PositiveInfinity;
            p_min.Z = float.PositiveInfinity;
            p_max.X = float.NegativeInfinity;
            p_max.Y = float.NegativeInfinity;
            p_max.Z = float.NegativeInfinity;

            foreach (var boneModel in Models)
            {
                var model = boneModel.Model;
                if (model.BoundingBox.max_x > p_max.X) p_max.X = model.BoundingBox.max_x;
                if (model.BoundingBox.max_y > p_max.Y) p_max.Y = model.BoundingBox.max_y;
                if (model.BoundingBox.max_z > p_max.Z) p_max.Z = model.BoundingBox.max_z;

                if (model.BoundingBox.min_x < p_min.X) p_min.X = model.BoundingBox.min_x;
                if (model.BoundingBox.min_y < p_min.Y) p_min.Y = model.BoundingBox.min_y;
                if (model.BoundingBox.min_z < p_min.Z) p_min.Z = model.BoundingBox.min_z;
            }

            if (float.IsInfinity(p_min.X))
                return 0;

            return CalculateDistance(p_max, p_min);
        }
    }

    /// <summary>
    /// Unified skeleton that works for both field and battle models.
    /// Abstracts the differences between formats while preserving the ability
    /// to save back to the original format.
    /// </summary>
    public class UnifiedSkeleton
    {
        /// <summary>
        /// Skeleton file name (without path).
        /// </summary>
        public string FileName;

        /// <summary>
        /// Skeleton name (from HRC header or generated).
        /// </summary>
        public string Name;

        /// <summary>
        /// Source format type.
        /// </summary>
        public SkeletonSourceType SourceType;

        /// <summary>
        /// Bone direction along Z axis.
        /// Field skeletons use NegativeZ, battle skeletons use PositiveZ.
        /// </summary>
        public BoneDirection BoneDirection;

        /// <summary>
        /// Whether root Y translation should be negated.
        /// True for field skeletons, false for battle.
        /// </summary>
        public bool NegateRootY;

        /// <summary>
        /// List of bones in hierarchy order.
        /// </summary>
        public List<UnifiedBone> Bones = [];

        /// <summary>
        /// Shared texture pool for all bones.
        /// </summary>
        public List<TEX> Textures = [];

        /// <summary>
        /// GPU texture IDs (populated after textures are uploaded).
        /// </summary>
        public uint[] TextureIDs = [];

        /// <summary>
        /// Weapon models (battle skeletons only).
        /// Weapons have their own animation track and are rendered separately.
        /// </summary>
        public List<PModel> Weapons = [];

        /// <summary>
        /// Quick lookup from bone name to index.
        /// </summary>
        public Dictionary<string, int> BoneNameToIndex = [];

        /// <summary>
        /// Indices of root bones (those with ParentIndex == -1).
        /// </summary>
        public List<int> RootBoneIndices = [];

        /// <summary>
        /// Number of bones.
        /// </summary>
        public int BoneCount => Bones?.Count ?? 0;

        public int WeaponCount => Weapons?.Count ?? 0;

        public int TextureCount => TextureIDs.Length;

        public bool IsBattleLocation { get; }

        #region Battle Skeleton Metadata (for round-trip conversion)
        /// <summary>
        /// Original skeleton type from battle skeleton (preserved for round-trip).
        /// </summary>
        public SkeletonType SkeletonType { get; private set; } = SkeletonType.EnemyOrSummon;

        /// <summary>
        /// Whether skeleton can have limit break (battle skeletons only).
        /// </summary>
        public bool CanHaveLimitBreak { get; private set; }

        /// <summary>
        /// Battle skeleton metadata - preserved for round-trip conversion.
        /// </summary>
        public int BattleUnk1 { get; private set; } = 1;
        public int BattleUnk2 { get; private set; } = 1;
        public int BattleUnk3 { get; private set; } = 0;
        public int BattleUnk4 { get; private set; } = 2;
        public int BattleUnk5 { get; private set; } = 0;
        public int BattleUnk6 { get; private set; } = 0;
        public int NumSkeletonAnims { get; private set; } = 0;
        public int WeaponAnimationCount { get; private set; } = 0;
        public int NumJoints { get; private set; } = 0;
        #endregion

        public UnifiedSkeleton(FieldSkeleton fs)
        {
            FileName = fs.fileName;
            Name = fs.name;
            SourceType = SkeletonSourceType.Field;
            BoneDirection = BoneDirection.NegativeZ;
            NegateRootY = true;

            // Copy texture pool
            if (fs.textures_pool != null)
            {
                Textures = fs.textures_pool;
            }

            // Build name-to-index map first (needed for parent lookup)
            var nameToIndex = new Dictionary<string, int>();
            if (fs.bones != null)
            {
                for (int i = 0; i < fs.bones.Count; i++)
                {
                    var fb = fs.bones[i];
                    if (!string.IsNullOrEmpty(fb.joint_i))
                    {
                        nameToIndex[fb.joint_i] = i;
                    }
                }

                // Convert bones
                for (int i = 0; i < fs.bones.Count; i++)
                {
                    var fb = fs.bones[i];
                    int index = -1;

                    // Resolve parent index from name
                    if (!string.IsNullOrEmpty(fb.joint_f) && nameToIndex.TryGetValue(fb.joint_f, out int parentIdx))
                    {
                        index = parentIdx;
                    }

                    Bones.Add(new UnifiedBone(fb, i, index));
                }
            }

            // Build lookup structures
            BuildLookups();
        }

        public UnifiedSkeleton(BattleSkeleton bs)
        {
            FileName = bs.fileName;
            Name = bs.fileName; // Battle skeletons don't have a separate name
            SourceType = bs.IsBattleLocation ? SkeletonSourceType.Battle :
                        (bs.fileName.EndsWith(".D", StringComparison.OrdinalIgnoreCase) ?
                            SkeletonSourceType.Magic : SkeletonSourceType.Battle);
            BoneDirection = BoneDirection.PositiveZ;
            NegateRootY = false;
            IsBattleLocation = bs.IsBattleLocation;

            // Preserve battle skeleton metadata for round-trip conversion
            SkeletonType = bs.skeletonType;
            CanHaveLimitBreak = bs.CanHaveLimitBreak;
            BattleUnk1 = bs.unk1;
            BattleUnk2 = bs.unk2;
            BattleUnk3 = bs.unk3;
            BattleUnk4 = bs.unk4;
            BattleUnk5 = bs.unk5;
            BattleUnk6 = bs.unk6;
            NumSkeletonAnims = bs.nsSkeletonAnims;
            WeaponAnimationCount = bs.nsWeaponsAnims;
            NumJoints = bs.nJoints;

            // Copy textures and texture IDs
            if (bs.textures != null)
            {
                Textures = bs.textures;
            }
            TextureIDs = bs.TexIDS;

            // Convert bones
            if (bs.bones != null)
            {
                for (int i = 0; i < bs.bones.Count; i++)
                {
                    Bones.Add(new UnifiedBone(bs.bones[i], i));
                }
            }

            // Copy weapon models
            if (bs.wpModels != null)
            {
                Weapons = bs.wpModels;
            }

            // Build lookup structures
            BuildLookups();
        }

        public UnifiedSkeleton(UnifiedSkeleton other)
        {
            FileName = other.FileName;
            Name = other.Name;
            SourceType = other.SourceType;
            BoneDirection = other.BoneDirection;
            NegateRootY = other.NegateRootY;

            foreach (var t in other.Textures)
            {
                Textures.Add(t);
            }
            TextureIDs = new uint[other.TextureIDs.Length];
            Array.Copy(other.TextureIDs, TextureIDs, TextureIDs.Length);

            foreach (var b in other.Bones)
            {
                Bones.Add(new UnifiedBone(b));
            }

            foreach (var w in other.Weapons)
            {
                Weapons.Add(w);
            }

            BuildLookups();
        }

        /// <summary>
        /// Build helper data structures (BoneNameToIndex, RootBoneIndices, ChildIndices).
        /// Call this after bones are populated.
        /// </summary>
        private void BuildLookups()
        {
            BoneNameToIndex.Clear();
            RootBoneIndices.Clear();

            // Build name-to-index lookup and clear child indices
            for (int i = 0; i < Bones.Count; i++)
            {
                var bone = Bones[i];
                bone.Index = i;
                bone.ChildIndices.Clear();

                if (!string.IsNullOrEmpty(bone.Name))
                {
                    BoneNameToIndex[bone.Name] = i;
                }

                if (bone.ParentIndex < 0)
                {
                    RootBoneIndices.Add(i);
                }
            }

            // Build child indices
            for (int i = 0; i < Bones.Count; i++)
            {
                var bone = Bones[i];
                if (bone.ParentIndex >= 0 && bone.ParentIndex < Bones.Count)
                {
                    Bones[bone.ParentIndex].ChildIndices.Add(i);
                }
            }
        }

        /// <summary>
        /// Get bone by name, or null if not found.
        /// </summary>
        public UnifiedBone? GetBoneByName(string name)
        {
            if (BoneNameToIndex.TryGetValue(name, out int index))
            {
                return Bones[index];
            }
            return null;
        }

        /// <summary>
        /// Get bone by index, or null if out of range.
        /// </summary>
        public UnifiedBone? GetBoneByIndex(int index)
        {
            if (index >= 0 && index < Bones.Count)
            {
                return Bones[index];
            }
            return null;
        }

        public string GetTextureFileName(int texIndex)
        {
            switch (SourceType)
            {
                case SkeletonSourceType.Battle:
                    return FileName.Substring(0, 2) + 'A' + Convert.ToChar('C' + texIndex);
            }
            return string.Empty;
        }

        /// <summary>
        /// Computes the world-space bounding box of the skeleton in a given animation frame.
        /// Works for both field and battle skeletons by using the unified bone structure.
        /// </summary>
        /// <param name="frame">The animation frame containing root transform and bone rotations.</param>
        /// <param name="p_min">Output: minimum corner of the bounding box.</param>
        /// <param name="p_max">Output: maximum corner of the bounding box.</param>
        public void ComputeBoundingBox(UnifiedFrame frame, ref Vector3 p_min, ref Vector3 p_max)
        {
            if (Bones == null || Bones.Count == 0 || frame == null)
            {
                p_min = Vector3.Zero;
                p_max = Vector3.Zero;
                return;
            }

            // Initialize output bounds to "empty"
            p_min.X = float.PositiveInfinity;
            p_min.Y = float.PositiveInfinity;
            p_min.Z = float.PositiveInfinity;
            p_max.X = float.NegativeInfinity;
            p_max.Y = float.NegativeInfinity;
            p_max.Z = float.NegativeInfinity;

            // Stack for parent tracking (integer-based, like battle skeleton)
            int[] parentStack = new int[Bones.Count + 1];
            Matrix4[] matrixStack = new Matrix4[Bones.Count + 2];
            int stackPtr = 0;
            int parentStackPtr = 0;

            parentStack[parentStackPtr] = -1; // Root parent is -1

            // Build initial transform from root translation and rotation
            Matrix4 currentMatrix = Matrix4.Identity;
            currentMatrix *= Matrix4.CreateTranslation(frame.RootTranslation);

            Matrix4 rotMatrix = BuildRotationMatrixWithQuaternions(
                frame.RootRotation.Alpha,
                frame.RootRotation.Beta,
                frame.RootRotation.Gamma);
            currentMatrix *= rotMatrix;

            matrixStack[stackPtr++] = currentMatrix;

            Vector3 p_min_bone = new();
            Vector3 p_max_bone = new();
            Vector3 p_min_bone_trans = new();
            Vector3 p_max_bone_trans = new();

            // Determine bone translation direction based on skeleton type
            float boneDirectionSign = (BoneDirection == BoneDirection.NegativeZ) ? -1.0f : 1.0f;

            for (int bi = 0; bi < Bones.Count; bi++)
            {
                var bone = Bones[bi];

                // Pop matrix stack until we find our parent
                while (bone.ParentIndex != parentStack[parentStackPtr] && parentStackPtr > 0)
                {
                    stackPtr--;
                    currentMatrix = matrixStack[stackPtr];
                    parentStackPtr--;
                }

                // Push current matrix for children
                matrixStack[stackPtr++] = currentMatrix;

                // Apply bone rotation
                if (bi < frame.BoneRotations.Count)
                {
                    var boneRot = frame.BoneRotations[bi];
                    rotMatrix = BuildRotationMatrixWithQuaternions(boneRot.Alpha, boneRot.Beta, boneRot.Gamma);
                    currentMatrix *= rotMatrix;
                }

                // Compute this bone's bounding box
                bone.ComputeBoundingBox(ref p_min_bone, ref p_max_bone);

                // Transform the bone's bounding box to world space
                double[] mvMatrix = Matrix4ToDoubleArray(currentMatrix);
                ComputeTransformedBoxBoundingBox(mvMatrix, ref p_min_bone, ref p_max_bone,
                                                 ref p_min_bone_trans, ref p_max_bone_trans);

                // Accumulate into overall bounds
                if (p_max.X < p_max_bone_trans.X) p_max.X = p_max_bone_trans.X;
                if (p_max.Y < p_max_bone_trans.Y) p_max.Y = p_max_bone_trans.Y;
                if (p_max.Z < p_max_bone_trans.Z) p_max.Z = p_max_bone_trans.Z;

                if (p_min.X > p_min_bone_trans.X) p_min.X = p_min_bone_trans.X;
                if (p_min.Y > p_min_bone_trans.Y) p_min.Y = p_min_bone_trans.Y;
                if (p_min.Z > p_min_bone_trans.Z) p_min.Z = p_min_bone_trans.Z;

                // Translate along bone to position for children
                currentMatrix *= Matrix4.CreateTranslation(0, 0, boneDirectionSign * bone.Length);

                // Push bone index onto parent stack
                parentStackPtr++;
                parentStack[parentStackPtr] = bi;
            }

            // If no geometry was found, return zero bounds
            if (float.IsInfinity(p_min.X))
            {
                p_min = Vector3.Zero;
                p_max = Vector3.Zero;
            }
        }

        /// <summary>
        /// Computes a diameter value for the skeleton.
        /// For field skeletons: Sum of all bone model diameters.
        /// For battle skeletons: Maximum path length through bone hierarchy.
        /// </summary>
        /// <returns>Diameter value used for camera/view calculations.</returns>
        public float ComputeDiameter()
        {
            if (Bones == null || Bones.Count == 0)
                return 0;

            // Field skeletons: sum of bone diameters (model extents)
            if (SourceType == SkeletonSourceType.Field)
            {
                float totalDiameter = 0;
                foreach (var bone in Bones)
                {
                    totalDiameter += bone.ComputeDiameter();
                }
                return totalDiameter;
            }

            // Battle location: max bone length
            if (IsBattleLocation)
            {
                float maxLength = 0;
                foreach (var bone in Bones)
                {
                    if (bone.Length > maxLength)
                        maxLength = bone.Length;
                }
                return maxLength;
            }

            // Battle skeleton with single bone and zero length: use model diameter
            if (Bones.Count == 1 && Bones[0].Length <= 0 && Bones[0].HasModel)
            {
                return Bones[0].Models[0].Model.diameter;
            }

            // Battle skeleton: find longest path through hierarchy
            float maxPath = 0;
            float currentPath = 0;
            int[] parentStack = new int[Bones.Count + 1];
            int stackPtr = 0;

            parentStack[stackPtr] = -1;

            for (int bi = 0; bi < Bones.Count; bi++)
            {
                var bone = Bones[bi];

                // Pop stack until we find parent
                while (bone.ParentIndex != parentStack[stackPtr] && stackPtr > 0)
                {
                    currentPath += Bones[parentStack[stackPtr]].Length;
                    stackPtr--;
                }

                currentPath -= bone.Length;
                if (currentPath > maxPath)
                    maxPath = currentPath;

                stackPtr++;
                parentStack[stackPtr] = bi;
            }

            return maxPath;
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
        /// Tests ray intersection with a PModel's geometry.
        /// </summary>
        private static bool RayIntersectsModel(Vector3 rayOrigin, Vector3 rayDir,
                                                PModel model, Matrix4 modelTransform,
                                                out float minDist)
        {
            minDist = float.MaxValue;
            bool hit = false;

            if (model.Polys == null) return false;

            for (int gi = 0; gi < model.Header.numGroups; gi++)
            {
                if (model.Groups[gi].HiddenQ) continue;

                int offsetVert = model.Groups[gi].offsetVert;

                for (int pi = model.Groups[gi].offsetPoly;
                     pi < model.Groups[gi].offsetPoly + model.Groups[gi].numPoly;
                     pi++)
                {
                    // Transform vertices by the model transform
                    Vector4 v0h = new Vector4(
                        model.Verts[model.Polys[pi].Verts[0] + offsetVert].X,
                        model.Verts[model.Polys[pi].Verts[0] + offsetVert].Y,
                        model.Verts[model.Polys[pi].Verts[0] + offsetVert].Z, 1.0f) * modelTransform;
                    Vector4 v1h = new Vector4(
                        model.Verts[model.Polys[pi].Verts[1] + offsetVert].X,
                        model.Verts[model.Polys[pi].Verts[1] + offsetVert].Y,
                        model.Verts[model.Polys[pi].Verts[1] + offsetVert].Z, 1.0f) * modelTransform;
                    Vector4 v2h = new Vector4(
                        model.Verts[model.Polys[pi].Verts[2] + offsetVert].X,
                        model.Verts[model.Polys[pi].Verts[2] + offsetVert].Y,
                        model.Verts[model.Polys[pi].Verts[2] + offsetVert].Z, 1.0f) * modelTransform;

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

        /// <summary>
        /// Builds a model transform matrix from a PModel's local transform properties.
        /// Uses Euler angles (ZXY order) which works for both field and battle models.
        /// </summary>
        private static Matrix4 BuildModelTransform(PModel model, Matrix4 boneTransform)
        {
            // Build model transform (pre-multiply order to match OpenGL)
            // Scale, then rotation (ZXY order), then translation, then bone transform
            return Matrix4.CreateScale(model.resizeX, model.resizeY, model.resizeZ)
                * Matrix4.CreateRotationZ(MathHelper.DegreesToRadians(model.rotateGamma))
                * Matrix4.CreateRotationX(MathHelper.DegreesToRadians(model.rotateAlpha))
                * Matrix4.CreateRotationY(MathHelper.DegreesToRadians(model.rotateBeta))
                * Matrix4.CreateTranslation(model.repositionX, model.repositionY, model.repositionZ)
                * boneTransform;
        }

        /// <summary>
        /// Tests ray intersection with all models in a unified bone.
        /// </summary>
        private static bool RayIntersectsBone(Vector3 rayOrigin, Vector3 rayDir,
                                               UnifiedBone bone, Matrix4 boneTransform,
                                               out float minDist)
        {
            minDist = float.MaxValue;
            bool hit = false;

            foreach (var boneModel in bone.Models)
            {
                var model = boneModel.Model;
                if (model.Polys == null) continue;

                Matrix4 fullTransform = BuildModelTransform(model, boneTransform);

                if (RayIntersectsModel(rayOrigin, rayDir, model, fullTransform, out float dist))
                {
                    if (dist < minDist)
                    {
                        minDist = dist;
                        hit = true;
                    }
                }
            }

            return hit;
        }

        /// <summary>
        /// Finds the bone closest to a screen point using ray casting.
        /// Works for both field and battle skeletons.
        /// </summary>
        /// <param name="frame">The current animation frame.</param>
        /// <param name="px">Screen X coordinate.</param>
        /// <param name="py">Screen Y coordinate.</param>
        /// <param name="weaponFrame">Optional weapon animation frame (for battle skeletons with weapons).</param>
        /// <param name="weaponIndex">Index of weapon to test, or -1 for no weapon testing.</param>
        /// <returns>Bone index, or BoneCount if weapon was hit, or -1 if nothing hit.</returns>
        public int GetClosestBone(UnifiedFrame frame, int px, int py,
                                   UnifiedFrame? weaponFrame = null, int weaponIndex = -1)
        {
            if (Bones == null || Bones.Count == 0 || frame == null)
                return -1;

            // Get viewport
            int[] vp = new int[4];
            GL.GetInteger(GetPName.Viewport, vp);
            int height = vp[3];

            // Get view and projection matrices from GLRenderer
            Matrix4 view = ViewMatrix;
            Matrix4 projection = ProjectionMatrix;

            // Create ray from screen coordinates
            Vector4 viewport = new Vector4(vp[0], vp[1], vp[2], vp[3]);
            float screenY = height - py;

            // Unproject to create ray (using identity for model since we transform vertices manually)
            Vector3 nearPoint = Unproject(new Vector3(px, screenY, 0.0f), Matrix4.Identity, view, projection, viewport);
            Vector3 farPoint = Unproject(new Vector3(px, screenY, 1.0f), Matrix4.Identity, view, projection, viewport);
            Vector3 rayOrigin = nearPoint;
            Vector3 rayDir = Vector3.Normalize(farPoint - nearPoint);

            // Build root transform (pre-multiply to match OpenGL)
            Matrix4 rotMatrix = BuildRotationMatrixWithQuaternions(
                frame.RootRotation.Alpha,
                frame.RootRotation.Beta,
                frame.RootRotation.Gamma);
            Matrix4 rootTransform = rotMatrix * Matrix4.CreateTranslation(frame.RootTranslation);

            // Setup matrix stack for hierarchy traversal
            int[] parentStack = new int[Bones.Count + 1];
            Matrix4[] matrixStack = new Matrix4[Bones.Count + 2];
            int stackPtr = 0;
            int parentStackPtr = 0;

            parentStack[parentStackPtr] = -1;
            Matrix4 currentMatrix = rootTransform;
            matrixStack[stackPtr++] = currentMatrix;

            int closestBone = -1;
            float closestDist = float.MaxValue;

            // Determine bone translation direction based on skeleton type
            float boneDirectionSign = (BoneDirection == BoneDirection.NegativeZ) ? -1.0f : 1.0f;

            for (int bi = 0; bi < Bones.Count; bi++)
            {
                var bone = Bones[bi];

                if (IsBattleLocation)
                {
                    // Battle location bones don't have hierarchy
                    if (RayIntersectsBone(rayOrigin, rayDir, bone, currentMatrix, out float dist))
                    {
                        if (dist < closestDist)
                        {
                            closestDist = dist;
                            closestBone = bi;
                        }
                    }
                }
                else
                {
                    // Pop matrix stack until we find our parent
                    while (bone.ParentIndex != parentStack[parentStackPtr] && parentStackPtr > 0)
                    {
                        stackPtr--;
                        currentMatrix = matrixStack[stackPtr];
                        parentStackPtr--;
                    }
                    matrixStack[stackPtr++] = currentMatrix;

                    // Apply bone rotation (pre-multiply to match OpenGL's transform order)
                    if (bi < frame.BoneRotations.Count)
                    {
                        var boneRot = frame.BoneRotations[bi];
                        rotMatrix = BuildRotationMatrixWithQuaternions(boneRot.Alpha, boneRot.Beta, boneRot.Gamma);
                        currentMatrix = rotMatrix * currentMatrix;
                    }

                    // Test ray intersection with this bone
                    if (RayIntersectsBone(rayOrigin, rayDir, bone, currentMatrix, out float dist))
                    {
                        if (dist < closestDist)
                        {
                            closestDist = dist;
                            closestBone = bi;
                        }
                    }

                    // Pre-multiply translation to match OpenGL's GL.Translate behavior
                    currentMatrix = Matrix4.CreateTranslation(0, 0, boneDirectionSign * bone.Length) * currentMatrix;

                    parentStackPtr++;
                    parentStack[parentStackPtr] = bi;
                }
            }

            // Test weapon if applicable (battle skeletons only)
            if (weaponIndex >= 0 && weaponIndex < Weapons.Count && weaponFrame != null)
            {
                var wpModel = Weapons[weaponIndex];
                if (wpModel.Polys != null)
                {
                    // Build weapon transform (pre-multiply to match OpenGL)
                    Matrix4 wpTransform = BuildRotationMatrixWithQuaternions(
                        weaponFrame.RootRotation.Alpha,
                        weaponFrame.RootRotation.Beta,
                        weaponFrame.RootRotation.Gamma)
                        * Matrix4.CreateTranslation(weaponFrame.RootTranslation);

                    // Apply model's local transforms (scale, rotation, translation in reverse order)
                    wpTransform = Matrix4.CreateScale(wpModel.resizeX, wpModel.resizeY, wpModel.resizeZ)
                        * Matrix4.CreateRotationZ(MathHelper.DegreesToRadians(wpModel.rotateGamma))
                        * Matrix4.CreateRotationX(MathHelper.DegreesToRadians(wpModel.rotateAlpha))
                        * Matrix4.CreateRotationY(MathHelper.DegreesToRadians(wpModel.rotateBeta))
                        * Matrix4.CreateTranslation(wpModel.repositionX, wpModel.repositionY, wpModel.repositionZ)
                        * wpTransform;

                    if (RayIntersectsModel(rayOrigin, rayDir, wpModel, wpTransform, out float dist))
                    {
                        if (dist < closestDist)
                        {
                            closestDist = dist;
                            closestBone = BoneCount; // Weapon is indexed after all bones
                        }
                    }
                }
            }

            return closestBone;
        }

        /// <summary>
        /// Finds the model piece (resource) closest to a screen point within a specific bone.
        /// Works for both field and battle skeletons.
        /// </summary>
        /// <param name="frame">The current animation frame.</param>
        /// <param name="boneIndex">The bone to search within.</param>
        /// <param name="px">Screen X coordinate.</param>
        /// <param name="py">Screen Y coordinate.</param>
        /// <returns>Model index within the bone, or -1 if nothing hit.</returns>
        public int GetClosestBonePiece(UnifiedFrame frame, int boneIndex, int px, int py)
        {
            if (Bones == null || boneIndex < 0 || boneIndex >= Bones.Count || frame == null)
                return -1;

            var targetBone = Bones[boneIndex];
            if (!targetBone.HasModel)
                return -1;

            // Get viewport
            int[] vp = new int[4];
            GL.GetInteger(GetPName.Viewport, vp);
            int height = vp[3];

            // Get view and projection matrices from GLRenderer
            Matrix4 view = ViewMatrix;
            Matrix4 projection = ProjectionMatrix;

            // Create ray from screen coordinates
            Vector4 viewport = new Vector4(vp[0], vp[1], vp[2], vp[3]);
            float screenY = height - py;

            // Unproject to create ray
            Vector3 nearPoint = Unproject(new Vector3(px, screenY, 0.0f), Matrix4.Identity, view, projection, viewport);
            Vector3 farPoint = Unproject(new Vector3(px, screenY, 1.0f), Matrix4.Identity, view, projection, viewport);
            Vector3 rayOrigin = nearPoint;
            Vector3 rayDir = Vector3.Normalize(farPoint - nearPoint);

            // Compute bone transform by traversing the hierarchy up to boneIndex
            Matrix4 boneTransform = ComputeBoneTransform(frame, boneIndex);

            // Test each model in the bone
            int closestModel = -1;
            float closestDist = float.MaxValue;

            for (int mi = 0; mi < targetBone.Models.Count; mi++)
            {
                var boneModel = targetBone.Models[mi];
                var model = boneModel.Model;
                if (model.Polys == null) continue;

                Matrix4 modelTransform = BuildModelTransform(model, boneTransform);

                if (RayIntersectsModel(rayOrigin, rayDir, model, modelTransform, out float dist))
                {
                    if (dist < closestDist)
                    {
                        closestDist = dist;
                        closestModel = mi;
                    }
                }
            }

            return closestModel;
        }

        /// <summary>
        /// Computes the world transform for a specific bone by traversing the hierarchy.
        /// </summary>
        /// <param name="frame">The current animation frame.</param>
        /// <param name="boneIndex">The target bone index.</param>
        /// <returns>The cumulative transform matrix for the bone.</returns>
        public Matrix4 ComputeBoneTransform(UnifiedFrame frame, int boneIndex)
        {
            if (Bones == null || boneIndex < 0 || boneIndex >= Bones.Count || frame == null)
                return Matrix4.Identity;

            // Build root transform (pre-multiply to match OpenGL)
            Matrix4 rotMatrix = BuildRotationMatrixWithQuaternions(
                frame.RootRotation.Alpha,
                frame.RootRotation.Beta,
                frame.RootRotation.Gamma);
            Matrix4 rootTransform = rotMatrix * Matrix4.CreateTranslation(frame.RootTranslation);

            // Setup matrix stack
            int[] parentStack = new int[Bones.Count + 1];
            Matrix4[] matrixStack = new Matrix4[Bones.Count + 2];
            int stackPtr = 0;
            int parentStackPtr = 0;

            parentStack[parentStackPtr] = -1;
            Matrix4 currentMatrix = rootTransform;
            matrixStack[stackPtr++] = currentMatrix;

            float boneDirectionSign = (BoneDirection == BoneDirection.NegativeZ) ? -1.0f : 1.0f;

            for (int bi = 0; bi <= boneIndex; bi++)
            {
                var bone = Bones[bi];

                // Pop matrix stack until we find our parent
                while (bone.ParentIndex != parentStack[parentStackPtr] && parentStackPtr > 0)
                {
                    stackPtr--;
                    currentMatrix = matrixStack[stackPtr];
                    parentStackPtr--;
                }
                matrixStack[stackPtr++] = currentMatrix;

                // Apply bone rotation
                if (bi < frame.BoneRotations.Count)
                {
                    var boneRot = frame.BoneRotations[bi];
                    rotMatrix = BuildRotationMatrixWithQuaternions(boneRot.Alpha, boneRot.Beta, boneRot.Gamma);
                    currentMatrix = rotMatrix * currentMatrix;
                }

                if (bi < boneIndex)
                {
                    // Pre-multiply translation
                    currentMatrix = Matrix4.CreateTranslation(0, 0, boneDirectionSign * bone.Length) * currentMatrix;
                    parentStackPtr++;
                    parentStack[parentStackPtr] = bi;
                }
            }

            return currentMatrix;
        }

        public bool AnimFileHasSameBoneCount(string other)
        {
            int boneCount;
            if (SourceType == SkeletonSourceType.Field)
            {
                boneCount = GetNumFieldBones(other);
                return (BoneCount == 1 && boneCount == 0) || BoneCount == boneCount;
            }
            else
            {
                boneCount = GetNumBattleBones(other);
                if (boneCount > 1 && SourceType == SkeletonSourceType.Magic)
                    boneCount--;
                return boneCount == BoneCount;
            }
        }

        /// <summary>
        /// Adds a PModel to a bone in the skeleton.
        /// Generates appropriate file/resource names based on the skeleton's source type.
        /// </summary>
        /// <param name="boneIndex">Index of the bone to add the model to.</param>
        /// <param name="model">The PModel to add.</param>
        public void AddBoneModel(int boneIndex, PModel model)
        {
            if (boneIndex < 0 || boneIndex >= Bones.Count)
                return;

            var bone = Bones[boneIndex];
            var boneModel = new UnifiedBoneModel(model);

            if (SourceType == SkeletonSourceType.Field)
            {
                // Field skeleton: generate RSD-style resource names
                if (bone.Models.Count > 0)
                {
                    string baseName = bone.Models[0].ResourceFile;
                    if (!string.IsNullOrEmpty(baseName) && baseName.Length >= 4)
                    {
                        baseName = baseName.Substring(0, 4).ToUpper();
                    }
                    else if (!string.IsNullOrEmpty(bone.Models[0].Model.fileName) &&
                             bone.Models[0].Model.fileName.Length >= 4)
                    {
                        baseName = bone.Models[0].Model.fileName.Substring(0, 4).ToUpper();
                    }
                    else
                    {
                        baseName = "BONE";
                    }

                    boneModel.Model.fileName = baseName + bone.Models.Count.ToString() + ".P";
                    boneModel.ResourceFile = baseName + bone.Models.Count.ToString();
                }
                else
                {
                    // First model on this bone - use model's filename if available
                    if (!string.IsNullOrEmpty(model.fileName) && model.fileName.Length >= 4)
                    {
                        boneModel.ResourceFile = model.fileName.Substring(0, 4).ToUpper();
                    }
                    else
                    {
                        boneModel.ResourceFile = "BONE";
                    }
                }
            }
            else
            {
                // Battle/Magic skeleton: generate model filenames based on type
                if (bone.Models.Count > 0)
                {
                    string baseFileName = bone.Models[0].Model.fileName;
                    if (!string.IsNullOrEmpty(baseFileName))
                    {
                        if (SourceType == SkeletonSourceType.Battle)
                        {
                            // AA skeleton style: append number directly
                            boneModel.Model.fileName = baseFileName + (bone.Models.Count - 1).ToString();
                        }
                        else
                        {
                            // Magic style: insert number before extension
                            int dotIndex = baseFileName.IndexOf('.');
                            if (dotIndex > 0)
                            {
                                boneModel.Model.fileName = baseFileName.Substring(0, dotIndex) +
                                                          ".P" + (bone.Models.Count - 1).ToString();
                            }
                            else
                            {
                                boneModel.Model.fileName = baseFileName + (bone.Models.Count - 1).ToString();
                            }
                        }
                    }
                }
            }

            bone.Models.Add(boneModel);
        }

        /// <summary>
        /// Removes a model from a bone in the skeleton.
        /// </summary>
        /// <param name="boneIndex">Index of the bone to remove the model from.</param>
        /// <param name="modelIndex">Index of the model within the bone to remove.</param>
        public void RemoveBoneModel(int boneIndex, int modelIndex)
        {
            if (boneIndex < 0 || boneIndex >= Bones.Count)
                return;

            var bone = Bones[boneIndex];
            if (modelIndex >= 0 && modelIndex < bone.Models.Count)
            {
                bone.Models.RemoveAt(modelIndex);
            }
        }

        /// <summary>
        /// Applies the current animation frame transforms to all bone models, baking the pose into vertex data.
        /// This traverses the bone hierarchy, applies rotation/scale transforms, and calls ApplyPChanges on each model.
        /// </summary>
        /// <param name="frame">The animation frame providing bone rotations.</param>
        /// <param name="weaponFrame">Optional weapon animation frame (for battle skeletons with weapons).</param>
        /// <param name="merge">If true, merges multiple models per bone into one (field skeleton behavior).</param>
        public void ApplyChanges(UnifiedFrame frame, UnifiedFrame? weaponFrame = null, bool merge = false)
        {
            if (Bones == null || Bones.Count == 0 || frame == null)
                return;

            int bi, jsp;
            int[] jointStack = new int[BoneCount + 2];
            double[] rotMat = new double[16];

            jsp = 0;
            jointStack[jsp] = -1;

            GL.MatrixMode(MatrixMode.Modelview);

            // For battle skeletons, apply root translation and rotation first
            if (SourceType != SkeletonSourceType.Field)
            {
                GL.PushMatrix();
                GL.Translated(frame.RootTranslation.X, frame.RootTranslation.Y, frame.RootTranslation.Z);

                BuildRotationMatrixWithQuaternions(frame.RootRotation.Alpha, frame.RootRotation.Beta, frame.RootRotation.Gamma, ref rotMat);
                GL.MultMatrixd(rotMat);
            }

            // Determine bone direction for translation
            float boneDirectionSign = (BoneDirection == BoneDirection.NegativeZ) ? -1.0f : 1.0f;

            for (bi = 0; bi < BoneCount; bi++)
            {
                var bone = Bones[bi];

                // Pop matrix stack until we find this bone's parent
                while (bone.ParentIndex != jointStack[jsp] && jsp > 0)
                {
                    GL.PopMatrix();
                    jsp--;
                }
                GL.PushMatrix();

                // Apply bone rotation from animation frame
                if (SourceType == SkeletonSourceType.Field)
                {
                    // Field: use Euler angles directly
                    if (bi < frame.BoneRotations.Count)
                    {
                        var rot = frame.BoneRotations[bi];
                        GL.Rotated(rot.Beta, 0, 1, 0);
                        GL.Rotated(rot.Alpha, 1, 0, 0);
                        GL.Rotated(rot.Gamma, 0, 0, 1);
                    }
                }
                else
                {
                    // Battle: use quaternion-based rotation matrix
                    if (bi < frame.BoneRotations.Count)
                    {
                        var rot = frame.BoneRotations[bi];
                        BuildRotationMatrixWithQuaternions(rot.Alpha, rot.Beta, rot.Gamma, ref rotMat);
                        GL.MultMatrixd(rotMat);
                    }
                }

                // Apply changes to each model on this bone
                if (bone.HasModel)
                {
                    ApplyBoneModelChanges(bone, merge);
                }

                // Translate along bone for children
                GL.Translated(0, 0, boneDirectionSign * bone.Length);

                jsp++;
                jointStack[jsp] = bi;
            }

            // Pop remaining matrices
            while (jsp > 0)
            {
                GL.PopMatrix();
                jsp--;
            }

            // Pop root matrix for battle skeletons
            if (SourceType != SkeletonSourceType.Field)
            {
                GL.PopMatrix();
            }

            // Handle weapons for battle skeletons
            if (SourceType != SkeletonSourceType.Field && Weapons != null && Weapons.Count > 0 && weaponFrame != null)
            {
                GL.MatrixMode(MatrixMode.Modelview);
                GL.PushMatrix();

                GL.Translated(weaponFrame.RootTranslation.X, weaponFrame.RootTranslation.Y, weaponFrame.RootTranslation.Z);
                BuildRotationMatrixWithQuaternions(weaponFrame.RootRotation.Alpha, weaponFrame.RootRotation.Beta, weaponFrame.RootRotation.Gamma, ref rotMat);
                GL.MultMatrixd(rotMat);

                for (int wi = 0; wi < Weapons.Count; wi++)
                {
                    if (Weapons[wi].Polys != null)
                    {
                        var wpModel = Weapons[wi];
                        ApplyWeaponModelChanges(ref wpModel);
                        Weapons[wi] = wpModel;
                    }
                }

                GL.PopMatrix();
            }
        }

        /// <summary>
        /// Applies transforms to all models attached to a bone, baking the current pose into vertex data.
        /// </summary>
        private void ApplyBoneModelChanges(UnifiedBone bone, bool merge)
        {
            for (int mi = 0; mi < bone.Models.Count; mi++)
            {
                var boneModel = bone.Models[mi];
                var model = boneModel.Model;

                if (model.Polys == null)
                    continue;

                // Apply vertex colors if lighting is enabled
                if (GL.IsEnabled(EnableCap.Lighting))
                {
                    ApplyCurrentVColors(ref model);
                }

                GL.MatrixMode(MatrixMode.Modelview);
                GL.PushMatrix();

                // Apply model's local transform
                SetCameraModelViewQuat(model.repositionX, model.repositionY, model.repositionZ,
                                       model.rotationQuaternion,
                                       model.resizeX, model.resizeY, model.resizeZ);

                // Apply bone scale
                GL.Scaled(bone.Scale.X, bone.Scale.Y, bone.Scale.Z);

                // Bake transforms into vertex data
                ApplyPChanges(ref model, false);

                GL.MatrixMode(MatrixMode.Modelview);
                GL.PopMatrix();

                boneModel.Model = model;
                bone.Models[mi] = boneModel;
            }

            // Merge models if requested (field skeleton behavior) or always for battle
            if (merge || SourceType != SkeletonSourceType.Field)
            {
                MergeBoneModels(bone);
            }
        }

        /// <summary>
        /// Applies transforms to a weapon model, baking the current pose into vertex data.
        /// </summary>
        private static void ApplyWeaponModelChanges(ref PModel wpModel)
        {
            if (GL.IsEnabled(EnableCap.Lighting))
            {
                ApplyCurrentVColors(ref wpModel);
            }

            GL.MatrixMode(MatrixMode.Modelview);
            GL.PushMatrix();

            SetCameraModelView(wpModel.repositionX, wpModel.repositionY, wpModel.repositionZ,
                               wpModel.rotateAlpha, wpModel.rotateBeta, wpModel.rotateGamma,
                               wpModel.resizeX, wpModel.resizeY, wpModel.resizeZ);

            GL.Scaled(wpModel.resizeX, wpModel.resizeY, wpModel.resizeZ);

            ApplyPChanges(ref wpModel, true);

            GL.MatrixMode(MatrixMode.Modelview);
            GL.PopMatrix();
        }

        /// <summary>
        /// Merges all models on a bone into a single model.
        /// </summary>
        private static void MergeBoneModels(UnifiedBone bone)
        {
            if (bone.Models.Count <= 1)
                return;

            // Merge all models into the first one
            var firstModel = bone.Models[0].Model;
            for (int mi = 1; mi < bone.Models.Count; mi++)
            {
                MergePModels(ref firstModel, bone.Models[mi].Model);
            }
            bone.Models[0].Model = firstModel;

            // Remove extra models
            while (bone.Models.Count > 1)
            {
                bone.Models.RemoveAt(bone.Models.Count - 1);
            }
        }

        /// <summary>
        /// Invalidates GPU mesh caches for all models in the skeleton, forcing recreation on next render.
        /// This is the unified equivalent of CreateDListsFromFieldSkeleton and CreateDListsFromBattleSkeleton.
        /// </summary>
        public void CreateDLists()
        {
            if (Bones == null)
                return;

            // Invalidate cache for all bone models
            foreach (var bone in Bones)
            {
                if (bone.Models == null)
                    continue;

                for (int mi = 0; mi < bone.Models.Count; mi++)
                {
                    var model = bone.Models[mi].Model;
                    CreateDListsFromPModel(ref model);
                    bone.Models[mi].Model = model;
                }
            }

            // Invalidate cache for weapon models (battle skeletons)
            if (Weapons != null)
            {
                for (int wi = 0; wi < Weapons.Count; wi++)
                {
                    var wpModel = Weapons[wi];
                    CreateDListsFromPModel(ref wpModel);
                    Weapons[wi] = wpModel;
                }
            }
        }

        private int MoveToBone(UnifiedFrame frame, int boneIndex)
        {
            int iBoneIdx, jsp;
            int[] joint_stack = new int[BoneCount * 4 + 1];
            double[] rot_mat = new double[16];

            // In unified format, BoneRotations[] uses direct indices (no offset)
            // RootRotation is separate and applied in SelectBoneAndModel
            float boneDir = BoneDirection == BoneDirection.NegativeZ ? -1f : 1f;

            jsp = 0;
            joint_stack[jsp] = -1;

            GL.MatrixMode(MatrixMode.Modelview);

            for (iBoneIdx = 0; iBoneIdx < boneIndex; iBoneIdx++)
            {
                while (Bones[iBoneIdx].ParentIndex != joint_stack[jsp] && jsp > 0)
                {
                    GL.PopMatrix();
                    jsp--;
                }
                GL.PushMatrix();

                // Apply bone rotation with bounds checking (direct index, matching DrawUnifiedSkeleton)
                if (iBoneIdx < frame.BoneRotations.Count)
                {
                    BuildRotationMatrixWithQuaternions(frame.BoneRotations[iBoneIdx].Alpha,
                                                       frame.BoneRotations[iBoneIdx].Beta,
                                                       frame.BoneRotations[iBoneIdx].Gamma,
                                                       ref rot_mat);
                    GL.MultMatrixd(rot_mat);
                }

                GL.Translated(0, 0, boneDir * Bones[iBoneIdx].Length);

                jsp++;
                joint_stack[jsp] = iBoneIdx;
            }

            while (Bones[boneIndex].ParentIndex != joint_stack[jsp] && jsp > 0)
            {
                GL.PopMatrix();
                jsp--;
            }
            GL.PushMatrix();

            // Apply final bone rotation with bounds checking
            if (boneIndex < frame.BoneRotations.Count)
            {
                BuildRotationMatrixWithQuaternions(frame.BoneRotations[boneIndex].Alpha,
                                                   frame.BoneRotations[boneIndex].Beta,
                                                   frame.BoneRotations[boneIndex].Gamma,
                                                   ref rot_mat);
                GL.MultMatrixd(rot_mat);
            }

            return jsp + 1;
        }

        public int MoveToBoneMiddle(UnifiedFrame frame, int boneIndex)
        {
            int iMoveToBattleBoneMiddleResult;

            iMoveToBattleBoneMiddleResult = MoveToBone(frame, boneIndex);
            GL.Translated(0, 0, Bones[boneIndex].Length / 2);

            return iMoveToBattleBoneMiddleResult;
        }

        public int MoveToBoneEnd(UnifiedFrame frame, int boneIndex)
        {
            int iMoveToBattleBoneEndResult;

            iMoveToBattleBoneEndResult = MoveToBone(frame, boneIndex);
            GL.Translated(0, 0, Bones[boneIndex].Length);

            return iMoveToBattleBoneEndResult;
        }

        /// <summary>
        /// Selects and highlights a bone and optionally a model piece by drawing bounding boxes.
        /// For battle skeletons with weapons, also handles weapon selection.
        /// </summary>
        /// <param name="frame">The current animation frame.</param>
        /// <param name="boneIndex">The bone index to select (-1 for none).</param>
        /// <param name="modelIndex">The model/piece index within the bone to select (-1 for none).</param>
        /// <param name="weaponFrame">Optional weapon animation frame (for battle skeletons).</param>
        /// <param name="weaponIndex">The weapon index to draw (for battle skeletons).</param>
        public void SelectBoneAndModel(UnifiedFrame frame, int boneIndex, int modelIndex,
                                       UnifiedFrame? weaponFrame = null, int weaponIndex = 0)
        {
            int i, jsp;
            double[] rot_mat = new double[16];

            GL.MatrixMode(MatrixMode.Modelview);
            GL.PushMatrix();

            // Apply root transform (matching DrawUnifiedSkeleton exactly)
            GL.Translated(frame.RootTranslation.X, frame.RootTranslation.Y, frame.RootTranslation.Z);

            BuildRotationMatrixWithQuaternions(frame.RootRotation.Alpha,
                                               frame.RootRotation.Beta,
                                               frame.RootRotation.Gamma,
                                               ref rot_mat);
            GL.MultMatrixd(rot_mat);

            if (boneIndex > -1 && boneIndex < BoneCount)
            {
                jsp = MoveToBone(frame, boneIndex);

                // Draw bone bounding box
                DrawBoneBoundingBox(Bones[boneIndex]);

                // Draw model/piece bounding box if specified
                if (modelIndex > -1 && modelIndex < Bones[boneIndex].Models.Count)
                {
                    DrawModelBoundingBox(Bones[boneIndex], modelIndex);
                }

                for (i = 0; i <= jsp; i++) GL.PopMatrix();
            }

            GL.PopMatrix();

            // Weapon bounding box is drawn outside skeleton's transform context
            // (DrawUnifiedWeapon is called independently, not relative to skeleton root)
            if (boneIndex == BoneCount)
            {
                DrawWeaponBoundingBox(weaponFrame, weaponIndex);
            }
        }

        /// <summary>
        /// Draws a bounding box around all models attached to a bone.
        /// </summary>
        private void DrawBoneBoundingBox(UnifiedBone bone)
        {
            GL.Disable(EnableCap.DepthTest);
            GL.MatrixMode(MatrixMode.Modelview);

            // Apply bone scale
            GL.Scaled(bone.Scale.X, bone.Scale.Y, bone.Scale.Z);

            if (!bone.HasModel)
            {
                // No models - draw bone line
                float boneDir = BoneDirection == BoneDirection.NegativeZ ? -1f : 1f;
                GL.Color3f(SourceType == SkeletonSourceType.Field ? 1f : 0f,
                           SourceType == SkeletonSourceType.Field ? 0f : 1f, 0);
                GL.Begin(PrimitiveType.Lines);
                GL.Vertex3f(0, 0, 0);
                GL.Vertex3f(0, 0, boneDir * bone.Length);
                GL.End();

                GL.Enable(EnableCap.DepthTest);
                return;
            }

            // Compute combined bounding box of all models
            float max_x = float.NegativeInfinity;
            float max_y = float.NegativeInfinity;
            float max_z = float.NegativeInfinity;
            float min_x = float.PositiveInfinity;
            float min_y = float.PositiveInfinity;
            float min_z = float.PositiveInfinity;

            foreach (var boneModel in bone.Models)
            {
                var model = boneModel.Model;
                if (max_x < model.BoundingBox.max_x) max_x = model.BoundingBox.max_x;
                if (max_y < model.BoundingBox.max_y) max_y = model.BoundingBox.max_y;
                if (max_z < model.BoundingBox.max_z) max_z = model.BoundingBox.max_z;

                if (min_x > model.BoundingBox.min_x) min_x = model.BoundingBox.min_x;
                if (min_y > model.BoundingBox.min_y) min_y = model.BoundingBox.min_y;
                if (min_z > model.BoundingBox.min_z) min_z = model.BoundingBox.min_z;
            }

            // Draw bone bounding box - red for field, green for battle (matches old behavior)
            if (SourceType == SkeletonSourceType.Field)
                ModelDrawing.DrawBox(max_x, max_y, max_z, min_x, min_y, min_z, 1, 0, 0);
            else
                ModelDrawing.DrawBox(max_x, max_y, max_z, min_x, min_y, min_z, 0, 1, 0);

            GL.Enable(EnableCap.DepthTest);
        }

        /// <summary>
        /// Draws a bounding box around a specific model attached to a bone.
        /// </summary>
        private void DrawModelBoundingBox(UnifiedBone bone, int modelIndex)
        {
            if (modelIndex < 0 || modelIndex >= bone.Models.Count)
                return;

            var model = bone.Models[modelIndex].Model;
            if (model.Verts == null || model.Verts.Length == 0)
                return;

            double[] rot_mat = new double[16];

            GL.Disable(EnableCap.DepthTest);
            GL.MatrixMode(MatrixMode.Modelview);

            // Apply bone scale
            GL.Scaled(bone.Scale.X, bone.Scale.Y, bone.Scale.Z);

            // Apply model translation
            GL.Translated(model.repositionX, model.repositionY, model.repositionZ);

            // Apply model rotation (quaternion for field, Euler X-Y-Z for battle)
            if (SourceType == SkeletonSourceType.Field)
            {
                BuildMatrixFromQuaternion(model.rotationQuaternion, ref rot_mat);
                GL.MultMatrixd(rot_mat);
            }
            else
            {
                // Rotation order X-Y-Z to match DrawUnifiedBoneModel
                GL.Rotated(model.rotateAlpha, 1, 0, 0);
                GL.Rotated(model.rotateBeta, 0, 1, 0);
                GL.Rotated(model.rotateGamma, 0, 0, 1);
            }

            // Apply model scale
            GL.Scaled(model.resizeX, model.resizeY, model.resizeZ);

            // Draw in green for piece/model selection (matching old behavior)
            ModelDrawing.DrawBox(model.BoundingBox.max_x, model.BoundingBox.max_y, model.BoundingBox.max_z,
                                 model.BoundingBox.min_x, model.BoundingBox.min_y, model.BoundingBox.min_z,
                                 0, 1, 0);

            GL.Enable(EnableCap.DepthTest);
        }

        /// <summary>
        /// Draws a bounding box around a weapon model (battle skeletons only).
        /// </summary>
        private void DrawWeaponBoundingBox(UnifiedFrame? weaponFrame, int weaponIndex)
        {
            if (SourceType == SkeletonSourceType.Field || Weapons == null || Weapons.Count == 0 ||
                weaponFrame == null || weaponIndex < 0 || weaponIndex >= Weapons.Count)
                return;

            var wpModel = Weapons[weaponIndex];
            if (wpModel.Verts == null || wpModel.Verts.Length == 0)
                return;

            double[] rot_mat = new double[16];

            GL.PushMatrix();

            // Apply weapon frame transform
            GL.Translated(weaponFrame.RootTranslation.X, weaponFrame.RootTranslation.Y, weaponFrame.RootTranslation.Z);
            BuildRotationMatrixWithQuaternions(weaponFrame.RootRotation.Alpha,
                                               weaponFrame.RootRotation.Beta,
                                               weaponFrame.RootRotation.Gamma, ref rot_mat);
            GL.MultMatrixd(rot_mat);

            GL.PushMatrix();

            // Apply weapon model transform (X-Y-Z rotation order to match DrawUnifiedWeapon)
            GL.Translated(wpModel.repositionX, wpModel.repositionY, wpModel.repositionZ);
            GL.Rotated(wpModel.rotateAlpha, 1, 0, 0);
            GL.Rotated(wpModel.rotateBeta, 0, 1, 0);
            GL.Rotated(wpModel.rotateGamma, 0, 0, 1);
            GL.Scaled(wpModel.resizeX, wpModel.resizeY, wpModel.resizeZ);

            // Draw weapon bounding box (yellow, matching old DrawPModelBoundingBox behavior)
            GL.Disable(EnableCap.DepthTest);
            ModelDrawing.DrawBox(wpModel.BoundingBox.max_x, wpModel.BoundingBox.max_y, wpModel.BoundingBox.max_z,
                                 wpModel.BoundingBox.min_x, wpModel.BoundingBox.min_y, wpModel.BoundingBox.min_z,
                                 1, 1, 0);
            GL.Enable(EnableCap.DepthTest);

            GL.PopMatrix();
            GL.PopMatrix();
        }

        /// <summary>
        /// Converts this unified skeleton back to a field skeleton structure.
        /// Note: This creates a new structure; textures are referenced, not deep copied.
        /// </summary>
        public FieldSkeleton ToFieldSkeleton()
        {
            var skeleton = new FieldSkeleton
            {
                fileName = FileName,
                name = Name,
                nBones = Bones?.Count ?? 0,
                bones = [],
                textures_pool = Textures ?? []
            };

            if (Bones != null)
            {
                foreach (var bone in Bones)
                {
                    skeleton.bones.Add(bone.ToFieldBone());
                }
            }

            return skeleton;
        }

        /// <summary>
        /// Converts this unified skeleton back to a battle skeleton structure.
        /// Note: This creates a new structure; textures and models are referenced, not deep copied.
        /// </summary>
        public BattleSkeleton ToBattleSkeleton()
        {
            var skeleton = new BattleSkeleton
            {
                fileName = FileName,
                skeletonType = SkeletonType,
                IsBattleLocation = IsBattleLocation,
                CanHaveLimitBreak = CanHaveLimitBreak,
                nBones = Bones?.Count ?? 0,
                nJoints = NumJoints > 0 ? NumJoints : (IsBattleLocation ? (Bones?.Count ?? 0) : 0),
                nTextures = Textures?.Count ?? 0,
                nWeapons = Weapons?.Count ?? 0,
                bones = [],
                textures = Textures ?? [],
                wpModels = Weapons ?? [],
                TexIDS = TextureIDs ?? [],

                // Preserved metadata from original (or sensible defaults)
                unk1 = BattleUnk1,
                unk2 = BattleUnk2,
                unk3 = BattleUnk3,
                unk4 = BattleUnk4,
                unk5 = BattleUnk5,
                unk6 = BattleUnk6,
                nsSkeletonAnims = NumSkeletonAnims,
                nsWeaponsAnims = WeaponAnimationCount
            };

            if (Bones != null)
            {
                foreach (var bone in Bones)
                {
                    skeleton.bones.Add(bone.ToBattleBone());
                }
            }

            return skeleton;
        }

        public void WriteSkeleton(string fileName, ModelType modelType)
        {
            if (modelType == ModelType.HRCSkeleton)
            {
                var fSkeleton = ToFieldSkeleton();
                WriteFieldSkeleton(ref fSkeleton, fileName);
            }
            else
            {
                var bSkeleton = ToBattleSkeleton();
                if (modelType == ModelType.AASkeleton)
                    WriteBattleSkeleton(ref bSkeleton, fileName);
                else 
                    WriteMagicSkeleton(ref bSkeleton, fileName);
            }
        }
    }
}
