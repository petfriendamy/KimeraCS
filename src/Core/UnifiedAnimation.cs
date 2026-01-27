using OpenTK.Mathematics;

namespace KimeraCS.Core
{
    using static FF7BattleAnimation;
    using static FF7BattleAnimationsPack;
    using static FF7FieldAnimation;
    using static Utils;

    /// <summary>
    /// Unified bone rotation for a single frame.
    /// Stores Euler angles in degrees, using FF7's Z->X->Y rotation order.
    /// </summary>
    public struct UnifiedBoneRotation
    {
        /// <summary>
        /// X-axis rotation (alpha) in degrees.
        /// </summary>
        public float Alpha;

        /// <summary>
        /// Y-axis rotation (beta) in degrees.
        /// </summary>
        public float Beta;

        /// <summary>
        /// Z-axis rotation (gamma) in degrees.
        /// </summary>
        public float Gamma;

        public UnifiedBoneRotation(float alpha, float beta, float gamma)
        {
            Alpha = alpha;
            Beta = beta;
            Gamma = gamma;
        }

        public UnifiedBoneRotation(FieldRotation rot) : this(rot.alpha, rot.beta, rot.gamma) { }

        public UnifiedBoneRotation(BattleFrameBone bone) : this(bone.alpha, bone.beta, bone.gamma) { }

        public UnifiedBoneRotation(UnifiedBoneRotation other) : this(other.Alpha, other.Beta, other.Gamma) { }

        public static UnifiedBoneRotation Zero => new UnifiedBoneRotation(0, 0, 0);

        public FieldRotation ToFieldRotation()
        {
            return new FieldRotation(Alpha, Beta, Gamma);
        }

        public BattleFrameBone ToBattleFrameBone()
        {
            return new BattleFrameBone(Alpha, Beta, Gamma);
        }
    }

    /// <summary>
    /// Unified animation frame containing root transform and per-bone rotations.
    /// </summary>
    public class UnifiedFrame
    {
        /// <summary>
        /// Root translation (world-space position of skeleton root).
        /// For field skeletons, Y is already negated during conversion.
        /// </summary>
        public Vector3 RootTranslation;

        /// <summary>
        /// Root rotation (applied before bone hierarchy).
        /// </summary>
        public UnifiedBoneRotation RootRotation;

        /// <summary>
        /// Per-bone rotations indexed by bone index.
        /// BoneRotations[i] corresponds to UnifiedSkeleton.Bones[i].
        /// </summary>
        public List<UnifiedBoneRotation> BoneRotations = [];

        public UnifiedFrame() { }

        public UnifiedFrame(FieldFrame frame)
        {
            // Convert bone rotations
            if (frame.rotations != null)
            {
                foreach (var rot in frame.rotations)
                {
                    BoneRotations.Add(new UnifiedBoneRotation(rot));
                }
            }

            RootTranslation = new Vector3(
                frame.rootTranslationX,
                -frame.rootTranslationY,
                frame.rootTranslationZ);
            RootRotation = new UnifiedBoneRotation(
                frame.rootRotationAlpha,
                frame.rootRotationBeta,
                frame.rootRotationGamma);
        }

        public UnifiedFrame(BattleFrame bf, int skeletonBoneCount)
        {
            // Battle root translation (no Y negation)
            RootTranslation = new Vector3(bf.startX, bf.startY, bf.startZ);

            // Root rotation is bones[0] for battle
            if (bf.bones != null && bf.bones.Count > 0)
            {
                RootRotation = new UnifiedBoneRotation(bf.bones[0]);

                // Bone rotations start at bones[1] for battle
                // We shift the indexing so BoneRotations[i] matches skeleton.Bones[i]
                for (int i = 0; i < skeletonBoneCount; i++)
                {
                    if (i + 1 < bf.bones.Count)
                    {
                        BoneRotations.Add(new UnifiedBoneRotation(bf.bones[i + 1]));
                    }
                    else
                    {
                        // Pad with zero rotation if animation has fewer bones
                        BoneRotations.Add(UnifiedBoneRotation.Zero);
                    }
                }
            }
        }

        public UnifiedFrame(UnifiedFrame other)
        {
            RootTranslation = other.RootTranslation;
            RootRotation = other.RootRotation;

            foreach (var rot in other.BoneRotations)
            {
                BoneRotations.Add(new UnifiedBoneRotation(rot));
            }
        }

        public FieldFrame ToFieldFrame()
        {
            var rotations = new List<FieldRotation>();
            foreach (var boneRot in BoneRotations)
            {
                rotations.Add(boneRot.ToFieldRotation());
            }

            return new FieldFrame(
                RootRotation.Alpha,
                RootRotation.Beta,
                RootRotation.Gamma,
                RootTranslation.X,
                -RootTranslation.Y,  // Negate Y back (was negated during conversion to Unified)
                RootTranslation.Z,
                rotations
            );
        }

        public BattleFrame ToBattleFrame()
        {
            var bones = new List<BattleFrameBone>
            {
                // Battle format: bones[0] = root rotation, bones[1..n] = bone rotations
                RootRotation.ToBattleFrameBone()
            };
            foreach (var boneRot in BoneRotations)
            {
                bones.Add(boneRot.ToBattleFrameBone());
            }

            return new BattleFrame
            {
                startX = (int)RootTranslation.X,
                startY = (int)RootTranslation.Y,  // No Y negation for battle animations
                startZ = (int)RootTranslation.Z,
                bones = bones
            };
        }
    }

    /// <summary>
    /// Unified animation that abstracts field and battle animation differences.
    /// </summary>
    public class UnifiedAnimation
    {
        /// <summary>
        /// Animation file name.
        /// </summary>
        public string FileName = string.Empty;

        /// <summary>
        /// Number of frames in the animation.
        /// </summary>
        public int FrameCount => Frames?.Count ?? 0;

        /// <summary>
        /// Number of bones this animation was created for.
        /// Should match the skeleton's bone count.
        /// </summary>
        public int BoneCount;

        /// <summary>
        /// Animation frames.
        /// </summary>
        public List<UnifiedFrame> Frames = [];

        public UnifiedAnimation() { }

        /// <summary>
        /// Constructor for a UnifiedAnimation, converted from a FieldAnimation.
        /// </summary>
        /// <param name="fAnimation">The FieldAnimation to convert from.</param>
        public UnifiedAnimation(FieldAnimation fAnimation)
        {
            FileName = fAnimation.strFieldAnimationFile;
            BoneCount = fAnimation.nBones;

            if (fAnimation.frames != null)
            {
                foreach (var ff in fAnimation.frames)
                {
                    Frames.Add(new UnifiedFrame(ff));
                }
            }
        }

        /// <summary>
        /// Constructor for a UnifiedAnimation, converted from a BattleAnimation.
        /// </summary>
        /// <param name="bAnimation">The BattleAnimation to convert from.</param>
        /// <param name="skeletonBoneCount">The bone count of the skeleton.</param>
        public UnifiedAnimation(BattleAnimation bAnimation, int skeletonBoneCount)
        {
            BoneCount = skeletonBoneCount;

            if (bAnimation.frames != null)
            {
                foreach (var bf in bAnimation.frames)
                {
                    Frames.Add(new UnifiedFrame(bf, skeletonBoneCount));
                }
            }
        }

        public UnifiedAnimation(UnifiedAnimation other)
        {
            FileName = other.FileName;
            BoneCount = other.BoneCount;

            foreach (var f in other.Frames)
            {
                Frames.Add(new UnifiedFrame(f));
            }
        }

        /// <summary>
        /// Get frame by index, clamped to valid range.
        /// </summary>
        public UnifiedFrame? GetFrame(int index)
        {
            if (Frames == null || Frames.Count == 0)
                return null;

            index = System.Math.Clamp(index, 0, Frames.Count - 1);
            return Frames[index];
        }

        /// <summary>
        /// Interpolates between two frames using quaternion SLERP for root rotation only.
        /// Use this for simple animations without bone hierarchy (e.g., weapon animations).
        /// </summary>
        /// <param name="frameA">First frame.</param>
        /// <param name="frameB">Second frame.</param>
        /// <param name="alpha">Interpolation factor (0 = frameA, 1 = frameB).</param>
        /// <param name="frameOut">Output frame to store the result.</param>
        public static void GetTwoFramesInterpolation(UnifiedFrame frameA, UnifiedFrame frameB,
                                                     float alpha, ref UnifiedFrame frameOut)
        {
            Quaterniond quatA, quatB, quatInterp;
            Vector3 eulerResult;
            float alphaInv = 1f - alpha;
            double[] mat = new double[16];

            // Interpolate root translation linearly
            frameOut.RootTranslation = new Vector3(
                frameA.RootTranslation.X * alphaInv + frameB.RootTranslation.X * alpha,
                frameA.RootTranslation.Y * alphaInv + frameB.RootTranslation.Y * alpha,
                frameA.RootTranslation.Z * alphaInv + frameB.RootTranslation.Z * alpha
            );

            // Interpolate root rotation using quaternion SLERP
            quatA = GetQuaternionFromEulerYXZr(frameA.RootRotation.Alpha, frameA.RootRotation.Beta, frameA.RootRotation.Gamma);
            NormalizeQuaternion(ref quatA);
            quatB = GetQuaternionFromEulerYXZr(frameB.RootRotation.Alpha, frameB.RootRotation.Beta, frameB.RootRotation.Gamma);
            NormalizeQuaternion(ref quatB);

            quatInterp = QuaternionSlerp2(ref quatA, ref quatB, alpha);
            NormalizeQuaternion(ref quatInterp);

            BuildMatrixFromQuaternion(quatInterp, ref mat);
            eulerResult = GetEulerYXZrFromMatrix(mat);

            frameOut.RootRotation = new UnifiedBoneRotation(eulerResult.Y, eulerResult.X, eulerResult.Z);

            // No bone rotations for simple/weapon animations
            frameOut.BoneRotations = [];
        }

        /// <summary>
        /// Interpolates between two frames using quaternion SLERP with proper bone hierarchy handling.
        /// This produces smoother interpolation than simple linear interpolation of Euler angles.
        /// </summary>
        /// <param name="skeleton">The skeleton, used for bone hierarchy information.</param>
        /// <param name="frameA">First frame.</param>
        /// <param name="frameB">Second frame.</param>
        /// <param name="alpha">Interpolation factor (0 = frameA, 1 = frameB).</param>
        /// <param name="frameOut">Output frame to store the result.</param>
        public static void GetTwoFramesInterpolation(UnifiedSkeleton skeleton, UnifiedFrame frameA, UnifiedFrame frameB,
                                                     float alpha, ref UnifiedFrame frameOut)
        {
            int bi, jsp;
            int numBones = skeleton.BoneCount;

            int[] parentStack = new int[numBones + 2];
            Quaterniond[] rotationsStackA = new Quaterniond[numBones + 1];
            Quaterniond[] rotationsStackB = new Quaterniond[numBones + 1];
            Quaterniond[] rotationsStackAccum = new Quaterniond[numBones + 1];

            Quaterniond quatA, quatB;
            Quaterniond quatAccumA = new();
            Quaterniond quatAccumB = new();
            Quaterniond quatAccumInverse;
            Quaterniond quatInterp;
            Quaterniond quatInterpFinal = new();

            Vector3 eulerResult;
            float alphaInv = 1f - alpha;
            double[] mat = new double[16];

            // Interpolate root rotation using quaternion SLERP
            quatA = GetQuaternionFromEulerYXZr(frameA.RootRotation.Alpha, frameA.RootRotation.Beta, frameA.RootRotation.Gamma);
            quatB = GetQuaternionFromEulerYXZr(frameB.RootRotation.Alpha, frameB.RootRotation.Beta, frameB.RootRotation.Gamma);

            quatInterp = QuaternionSlerp2(ref quatA, ref quatB, alpha);

            BuildMatrixFromQuaternion(quatInterp, ref mat);
            eulerResult = GetEulerYXZrFromMatrix(mat);

            frameOut.RootRotation = new UnifiedBoneRotation(eulerResult.Y, eulerResult.X, eulerResult.Z);

            // Interpolate root translation linearly
            frameOut.RootTranslation = new Vector3(
                frameA.RootTranslation.X * alphaInv + frameB.RootTranslation.X * alpha,
                frameA.RootTranslation.Y * alphaInv + frameB.RootTranslation.Y * alpha,
                frameA.RootTranslation.Z * alphaInv + frameB.RootTranslation.Z * alpha
            );

            // Initialize bone rotation output
            frameOut.BoneRotations = new List<UnifiedBoneRotation>();

            // Initialize stacks for hierarchical interpolation
            rotationsStackA[0] = quatA;
            rotationsStackB[0] = quatB;
            rotationsStackAccum[0] = quatInterp;

            jsp = 0;
            parentStack[jsp] = -1;

            // Interpolate each bone rotation
            for (bi = 0; bi < numBones; bi++)
            {
                var bone = skeleton.Bones[bi];

                // Pop stack until we find this bone's parent
                while (bone.ParentIndex != parentStack[jsp] && jsp > 0)
                {
                    jsp--;
                }

                // Get bone rotations from both frames
                if (bi < frameA.BoneRotations.Count && bi < frameB.BoneRotations.Count)
                {
                    var rotA = frameA.BoneRotations[bi];
                    var rotB = frameB.BoneRotations[bi];

                    // Convert to quaternions and accumulate with parent
                    quatA = GetQuaternionFromEulerYXZr(rotA.Alpha, rotA.Beta, rotA.Gamma);
                    MultiplyQuaternions(rotationsStackA[jsp], quatA, ref quatAccumA);
                    NormalizeQuaternion(ref quatAccumA);
                    rotationsStackA[jsp + 1] = quatAccumA;

                    quatB = GetQuaternionFromEulerYXZr(rotB.Alpha, rotB.Beta, rotB.Gamma);
                    MultiplyQuaternions(rotationsStackB[jsp], quatB, ref quatAccumB);
                    NormalizeQuaternion(ref quatAccumB);
                    rotationsStackB[jsp + 1] = quatAccumB;

                    // SLERP between accumulated rotations
                    quatInterp = QuaternionSlerp2(ref quatAccumA, ref quatAccumB, alpha);
                    rotationsStackAccum[jsp + 1] = quatInterp;

                    // Remove parent rotation to get local rotation
                    quatAccumInverse = GetQuaternionConjugate(ref rotationsStackAccum[jsp]);
                    MultiplyQuaternions(quatAccumInverse, quatInterp, ref quatInterpFinal);
                    NormalizeQuaternion(ref quatInterpFinal);

                    // Convert back to Euler angles
                    BuildMatrixFromQuaternion(quatInterpFinal, ref mat);
                    eulerResult = GetEulerYXZrFromMatrix(mat);

                    frameOut.BoneRotations.Add(new UnifiedBoneRotation(eulerResult.Y, eulerResult.X, eulerResult.Z));
                }
                else
                {
                    // If one frame doesn't have this bone, use zero rotation
                    frameOut.BoneRotations.Add(UnifiedBoneRotation.Zero);

                    // Still need to update stacks for hierarchy
                    rotationsStackA[jsp + 1] = rotationsStackA[jsp];
                    rotationsStackB[jsp + 1] = rotationsStackB[jsp];
                    rotationsStackAccum[jsp + 1] = rotationsStackAccum[jsp];
                }

                // Push this bone onto the stack
                jsp++;
                parentStack[jsp] = bi;
            }
        }

        /// <summary>
        /// Interpolates the animation by adding frames between existing keyframes.
        /// Uses quaternion SLERP with proper bone hierarchy handling.
        /// </summary>
        /// <param name="skeleton">The skeleton, used for bone hierarchy information.</param>
        /// <param name="numInterpolatedFrames">Number of frames to insert between each pair of keyframes.</param>
        /// <param name="isLoop">If true, also interpolates from the last frame back to the first.</param>
        public void InterpolateAnimation(UnifiedSkeleton skeleton, int numInterpolatedFrames, bool isLoop)
        {
            if (Frames == null || Frames.Count <= 1 || numInterpolatedFrames < 1)
                return;

            int nextElemDiff = numInterpolatedFrames + 1;
            int frameOffset = isLoop ? 0 : numInterpolatedFrames;

            // Calculate new frame count
            int originalFrameCount = Frames.Count;
            int newFrameCount = originalFrameCount * nextElemDiff - frameOffset;

            // Expand frames list
            while (Frames.Count < newFrameCount)
            {
                Frames.Add(new UnifiedFrame());
            }

            // Move original keyframes to their new positions (work backwards to avoid overwriting)
            for (int fi = newFrameCount - (1 + numInterpolatedFrames - frameOffset); fi >= 0; fi -= nextElemDiff)
            {
                Frames[fi] = new UnifiedFrame(Frames[fi / nextElemDiff]);
            }

            // Interpolate frames between keyframes
            for (int fi = 0; fi <= newFrameCount - (1 + nextElemDiff + numInterpolatedFrames - frameOffset); fi += nextElemDiff)
            {
                for (int ifi = 1; ifi <= numInterpolatedFrames; ifi++)
                {
                    float alpha = (float)ifi / nextElemDiff;

                    UnifiedFrame frameOut = Frames[fi + ifi];
                    GetTwoFramesInterpolation(skeleton, Frames[fi], Frames[fi + nextElemDiff], alpha, ref frameOut);
                    Frames[fi + ifi] = frameOut;
                }
            }

            // Handle looped animation - interpolate from last keyframe back to first
            if (isLoop)
            {
                int baseFinalFrame = newFrameCount - numInterpolatedFrames - 1;

                for (int ifi = 1; ifi <= numInterpolatedFrames; ifi++)
                {
                    float alpha = (float)ifi / nextElemDiff;

                    UnifiedFrame frameOut = Frames[baseFinalFrame + ifi];
                    GetTwoFramesInterpolation(skeleton, Frames[baseFinalFrame], Frames[0], alpha, ref frameOut);
                    Frames[baseFinalFrame + ifi] = frameOut;
                }
            }
        }

        /// <summary>
        /// Synchronizes this animation's frame count with another animation, then interpolates.
        /// Useful for weapon animations that need to match skeleton animation frame counts.
        /// </summary>
        /// <param name="targetFrameCount">The target frame count to match.</param>
        /// <param name="numInterpolatedFrames">Number of frames to insert between each pair of keyframes.</param>
        /// <param name="isLoop">If true, also interpolates from the last frame back to the first.</param>
        public void InterpolateAnimationToMatch(int targetFrameCount, int numInterpolatedFrames, bool isLoop)
        {
            if (Frames == null || Frames.Count == 0 || numInterpolatedFrames < 1)
                return;

            int nextElemDiff = numInterpolatedFrames + 1;
            int frameOffset = isLoop ? 0 : numInterpolatedFrames;

            // Expand frames list to target count
            while (Frames.Count < targetFrameCount)
            {
                Frames.Add(new UnifiedFrame());
            }

            // Move original keyframes to their new positions (work backwards to avoid overwriting)
            for (int fi = targetFrameCount - (1 + numInterpolatedFrames - frameOffset); fi >= 0; fi -= nextElemDiff)
            {
                Frames[fi] = new UnifiedFrame(Frames[fi / nextElemDiff]);
            }

            // Interpolate frames between keyframes
            for (int fi = 0; fi <= targetFrameCount - (1 + numInterpolatedFrames + numInterpolatedFrames - frameOffset); fi += nextElemDiff)
            {
                for (int ifi = 1; ifi <= numInterpolatedFrames; ifi++)
                {
                    float alpha = (float)ifi / nextElemDiff;

                    UnifiedFrame frameOut = Frames[fi + ifi];
                    GetTwoFramesInterpolation(Frames[fi], Frames[fi + nextElemDiff], alpha, ref frameOut);
                    Frames[fi + ifi] = frameOut;
                }

                // Handle looped animation within each segment (matching original behavior)
                if (isLoop)
                {
                    int baseFinalFrame = targetFrameCount - numInterpolatedFrames - 1;

                    for (int ifi = 1; ifi <= numInterpolatedFrames; ifi++)
                    {
                        float alpha = (float)ifi / nextElemDiff;

                        UnifiedFrame frameOut = Frames[baseFinalFrame + ifi];
                        GetTwoFramesInterpolation(Frames[baseFinalFrame], Frames[0], alpha, ref frameOut);
                        Frames[baseFinalFrame + ifi] = frameOut;
                    }
                }
            }
        }

        public FieldAnimation ToFieldAnimation()
        {
            var fAnimation = new FieldAnimation
            {
                version = 1,  // Must be 1 for FF7 to load
                nFrames = FrameCount,
                nBones = BoneCount,
                rotationOrder = [1, 0, 2],  // Y, X, Z order
                unused = 0,
                runtime_data = new int[5],
                frames = new List<FieldFrame>(),
                strFieldAnimationFile = FileName ?? "UNIFIED.A"
            };

            foreach (var uFrame in Frames)
            {
                fAnimation.frames.Add(uFrame.ToFieldFrame());
            }

            return fAnimation;
        }

        public BattleAnimation ToBattleAnimation()
        {
            var bAnimation = new BattleAnimation
            {
                nBones = BoneCount,
                numFrames = FrameCount,
                numFramesShort = (ushort)FrameCount,
                // Compression fields - will be computed when writing via WriteBattleAnimation
                blockSize = 0,
                blockSizeShort = 0,
                key = 0,
                framesRawData = Array.Empty<byte>(),
                padding4bytes = Array.Empty<byte>(),
                frames = new List<BattleFrame>()
            };

            foreach (var uFrame in Frames)
            {
                bAnimation.frames.Add(uFrame.ToBattleFrame());
            }

            return bAnimation;
        }
    }

    /// <summary>
    /// Unified animation pack for battle models that have multiple animations.
    /// Field models use a single UnifiedAnimation directly.
    /// </summary>
    public class UnifiedAnimationPack
    {
        /// <summary>
        /// Animation pack file name.
        /// </summary>
        public string FileName = string.Empty;

        /// <summary>
        /// Animation pack full file path.
        /// </summary>
        public string FilePath = string.Empty;

        /// <summary>
        /// Skeleton animations (main character/enemy animations).
        /// </summary>
        public List<UnifiedAnimation> SkeletonAnimations = [];

        /// <summary>
        /// Weapon animations (separate track for weapons).
        /// </summary>
        public List<UnifiedAnimation> WeaponAnimations = [];

        /// <summary>
        /// Whether this is a limit break animation pack.
        /// </summary>
        public bool IsLimitBreak;

        /// <summary>
        /// Number of skeleton animations.
        /// </summary>
        public int SkeletonAnimationCount => SkeletonAnimations?.Count ?? 0;

        /// <summary>
        /// Number of weapon animations.
        /// </summary>
        public int WeaponAnimationCount => WeaponAnimations?.Count ?? 0;

        public UnifiedAnimationPack() { }

        public UnifiedAnimationPack(BattleAnimationsPack bAnimationsPack, int skeletonBoneCount)
        {
            FileName = bAnimationsPack.strBattleAnimPackFileName;
            FilePath = bAnimationsPack.strAnimsPackFullFileName;
            IsLimitBreak = bAnimationsPack.IsLimit;

            // Convert skeleton animations
            if (bAnimationsPack.SkeletonAnimations != null)
            {
                foreach (var ba in bAnimationsPack.SkeletonAnimations)
                {
                    SkeletonAnimations.Add(new UnifiedAnimation(ba, skeletonBoneCount));
                }
            }

            // Convert weapon animations
            if (bAnimationsPack.WeaponAnimations != null)
            {
                foreach (var wa in bAnimationsPack.WeaponAnimations)
                {
                    // Weapon animations typically have 1 bone
                    WeaponAnimations.Add(new UnifiedAnimation(wa, 1));
                }
            }
        }

        public UnifiedAnimationPack(UnifiedAnimationPack other)
        {
            FileName = other.FileName;
            FilePath = other.FilePath;
            IsLimitBreak = other.IsLimitBreak;

            foreach (var sAnim in other.SkeletonAnimations)
            {
                SkeletonAnimations.Add(new UnifiedAnimation(sAnim));
            }

            foreach (var wAnim in other.WeaponAnimations)
            {
                WeaponAnimations.Add(new UnifiedAnimation(wAnim));
            }
        }

        /// <summary>
        /// Get skeleton animation by index, or null if out of range.
        /// </summary>
        public UnifiedAnimation? GetSkeletonAnimation(int index)
        {
            if (index >= 0 && index < SkeletonAnimationCount)
            {
                return SkeletonAnimations[index];
            }
            return null;
        }

        /// <summary>
        /// Get weapon animation by index, or null if out of range.
        /// </summary>
        public UnifiedAnimation? GetWeaponAnimation(int index)
        {
            if (index >= 0 && index < WeaponAnimationCount)
            {
                return WeaponAnimations[index];
            }
            return null;
        }

        public BattleAnimationsPack ToBattleAnimationsPack()
        {
            List<BattleAnimation> skeletonAnims = [],
                weaponAnims = [];

            foreach (var anim in SkeletonAnimations)
            {
                skeletonAnims.Add(anim.ToBattleAnimation());
            }

            foreach (var anim in WeaponAnimations)
            {
                weaponAnims.Add(anim.ToBattleAnimation());
            }

            return new BattleAnimationsPack
            {
                strBattleAnimPackFileName = FileName,
                strAnimsPackFullFileName = FilePath,
                SkeletonAnimations = skeletonAnims,
                WeaponAnimations = weaponAnims,
                IsLimit = IsLimitBreak,
                nbSkeletonAnims = skeletonAnims.Count,
                nbWeaponAnims = weaponAnims.Count,
                nAnimations = skeletonAnims.Count + weaponAnims.Count
            };
        }
    }
}
