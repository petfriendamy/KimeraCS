namespace KimeraCS.Core
{
    public enum ModelType
    {
        // This will tell to all the tool which type of skeleton/model we have loaded
        // Possible values:
        // -1: Error (file not exists, error opening file...)
        //  0: Field P Model
        //  1: Battle P Model
        //  2: Magic P Model
        //  3: Field Skeleton
        //  4: Battle Skeleton
        //  5: Magic Skeleton
        //  6: 3DS Model
        None = -1,
        PFieldModel = 0,
        PBattleModel = 1,
        PMagicModel = 2,
        HRCSkeleton = 3,
        AASkeleton = 4,
        MagicSkeleton = 5,
        ImportedModel = 6
    }
}
