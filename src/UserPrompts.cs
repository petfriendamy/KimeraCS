using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Windows.Forms;

namespace KimeraCS
{
    using static FF7BattleSkeleton;
    using static FF7FieldRSDResource;
    using static FF7FieldSkeleton;
    using static FF7PModel;
    using static FF7Skeleton;
    using static FF7TEXTexture;
    using static FileTools;

    public static class UserPrompts
    {
        // I don't expect this to ever be needed
        private const bool TEXTURE_REMOVE_CHECK = true;

        private static bool RepairPolysPrompt(ref PModel fPModel)
        {
            int polyCheck = CheckPolys(ref fPModel);
            if (polyCheck >= 0)
            {
                return (MessageBox.Show("The model: " + fPModel.fileName + " has one vertex index duplicated in " +
                    "the same triangle/poly. Do you want to fix it?\n(the triangle will be " +
                    "removed)\n\nNOTE: This answer will be used for the whole model. Selecting \"No\" will " +
                    "cancel the check of duplicated vertex coordinates.",
                    "Warning", MessageBoxButtons.YesNo) == DialogResult.Yes);
                /*
                 * "[INFO] Group:     " + iGroupIdx.ToString() + "\n" +
                 * "       Poly:      " + iPolyIdx.ToString() + "\n" +
                 * "       Vertex V0: " + iV0.ToString() + "\n" +
                 * "       Vertex V1: " + iV1.ToString() + "\n" +
                 * "       Vertex V2: " + iV2.ToString()
                 * 
                 * "[INFO] Group:     " + iGroupIdx.ToString() + "\n" +
                 * "       Poly:      " + iPolyIdx.ToString() + "\n" +
                 * "       Vertex V0: " + iV0.ToString() + "\n" +
                 * "       Vertex V1: " + iV1.ToString() + "\n" +
                 * "       Vertex V2: " + iV2.ToString() + "\n" +
                 * "       Vertex V0.x: " + Model.Verts[iV0].x.ToString() + "\n" +
                 * "       Vertex V0.y: " + Model.Verts[iV0].y.ToString() + "\n" +
                 * "       Vertex V0.z: " + Model.Verts[iV0].z.ToString() + "\n" +
                 * "       Vertex V1.x: " + Model.Verts[iV1].x.ToString() + "\n" +
                 * "       Vertex V1.y: " + Model.Verts[iV1].y.ToString() + "\n" +
                 * "       Vertex V1.z: " + Model.Verts[iV1].z.ToString() + "\n" +
                 * "       Vertex V2.x: " + Model.Verts[iV2].x.ToString() + "\n" +
                 * "       Vertex V2.y: " + Model.Verts[iV2].y.ToString() + "\n" +
                 * "       Vertex V2.z: " + Model.Verts[iV2].z.ToString() + "\n"
                 */
            }
            return false;
        }

        private static bool MissingPFilePrompt(string message)
        {
            return MessageBox.Show(message +
                                   "Do you want to continue loading the other parts of the model?",
                                   "Information", MessageBoxButtons.YesNo) == DialogResult.Yes;
        }

        private static bool RemoveTextureCoordsPrompt(int iGroupIdx)
        {
            return MessageBox.Show("Group " + iGroupIdx.ToString("00") + " seems to have Texture " +
                                   "Coordinates assigned, but it has not any Texture assigned.\nDo you " +
                                   "want -Reset to 0- the Texture Coordinates?\n\nNOTE: If you don't reset " +
                                   "them to 0, the Texture Flag will be enabled.",
                                   "Question", MessageBoxButtons.YesNo) == DialogResult.Yes;
        }

        public static int LoadGenericSkeleton(string strFileName, bool loadGeometryQ)
        {
            int result;
            try
            {
                result = LoadSkeleton(strFileName, loadGeometryQ, false, false, TEXTURE_REMOVE_CHECK);
            }
            catch (PFileNotFoundException ex) //missing P files
            {
                if (MissingPFilePrompt(ex.Message))
                    result = LoadSkeleton(strFileName, loadGeometryQ, true, false, TEXTURE_REMOVE_CHECK);
                else
                    throw new FileLoadException("File could not be loaded.", ex);
            }
            int modelType = GetSkeletonType(strFileName);
            switch (modelType)
            {
                case K_HRC_SKELETON:
                    FieldSkeletonPolyCheck(ref fSkeleton);
                    break;
                case K_AA_SKELETON:
                case K_MAGIC_SKELETON:
                    BattleSkeletonPolyCheck(ref bSkeleton);
                    break;
            }
            return result;
        }

        public static int LoadSkeletonFromDB(string strFileName, string strAnimFileName, bool loadGeometryQ)
        {
            int result;
            try
            {
                result = LoadFieldSkeletonFromDB(strFileName, strAnimFileName, loadGeometryQ, false, false, TEXTURE_REMOVE_CHECK);
            }
            catch (PFileNotFoundException ex) //missing P files
            {
                if (MissingPFilePrompt(ex.Message))
                    result = LoadSkeleton(strFileName, loadGeometryQ, true, false, TEXTURE_REMOVE_CHECK);
                else
                    throw new FileLoadException("File could not be loaded.", ex);
            }
            int modelType = GetSkeletonType(strFileName);
            switch (modelType)
            {
                case K_HRC_SKELETON:
                    FieldSkeletonPolyCheck(ref fSkeleton);
                    break;
                case K_AA_SKELETON:
                case K_MAGIC_SKELETON:
                    BattleSkeletonPolyCheck(ref bSkeleton);
                    break;
            }

            return result;
        }

        public static FieldSkeleton FieldSkeletonLoader(string strfileName, bool loadGeometryQ)
        {
            FieldSkeleton fSkeleton;
            try
            {
                fSkeleton = new FieldSkeleton(strfileName, loadGeometryQ, false, false, TEXTURE_REMOVE_CHECK);
            }
            catch (PFileNotFoundException ex) //missing P files
            {
                if (MissingPFilePrompt(ex.Message))
                    fSkeleton = new FieldSkeleton(strfileName, loadGeometryQ, true, false, TEXTURE_REMOVE_CHECK);
                else
                    throw new FileLoadException("File could not be loaded.", ex);
            }
            FieldSkeletonPolyCheck(ref fSkeleton);
            return fSkeleton;
        }

        public static FieldRSDResource FieldRSDLoader(string in_res_file, ref List<TEX> textures_pool, string strFolderName)
        {
            FieldRSDResource fRSD;
            try
            {
                fRSD = new FieldRSDResource(in_res_file, ref textures_pool, strFolderName, false, false, TEXTURE_REMOVE_CHECK);
            }
            catch (PFileNotFoundException ex) //missing P files
            {
                if (MissingPFilePrompt(ex.Message))
                    fRSD = new FieldRSDResource(in_res_file, ref textures_pool, strFolderName, true, false, TEXTURE_REMOVE_CHECK);
                else
                    throw new FileLoadException("File could not be loaded.", ex);
            }
            if (!bDontCheckRepairPolys)
            {
                if (RepairPolysPrompt(ref fRSD.Model))
                {
                    RepairPolys(ref fRSD.Model);
                }
            }
            return fRSD;
        }

        public static int FieldRSDLoader(string strRSDFolder, string strRSDName)
        {
            int result = -1;
            try
            {
                result = LoadRSDResourceModel(strRSDFolder, strRSDName, false, false, TEXTURE_REMOVE_CHECK);
            }
            catch (PFileNotFoundException ex) //missing P files
            {
                if (MissingPFilePrompt(ex.Message))
                    result = LoadRSDResourceModel(strRSDFolder, strRSDName, true, false, TEXTURE_REMOVE_CHECK);
                else
                    throw new FileLoadException("File could not be loaded.", ex);
            }
            FieldSkeletonPolyCheck(ref fSkeleton);
            return result;
        }

        private static void FieldSkeletonPolyCheck(ref FieldSkeleton fSkeleton)
        {
            if (!bDontCheckRepairPolys)
            {
                PModel? result = fSkeleton.CheckPolys();
                if (result != null)
                {
                    var p = (PModel)result;
                    if (RepairPolysPrompt(ref p))
                    {
                        fSkeleton.RepairPolys();
                    }
                }
            }
        }

        public static BattleSkeleton BattleSkeletonLoader(string strFullFileName, bool isLimitBreak)
        {
            var bSkeleton = new BattleSkeleton(strFullFileName, isLimitBreak, false, TEXTURE_REMOVE_CHECK);
            BattleSkeletonPolyCheck(ref bSkeleton);
            return bSkeleton;
        }

        public static BattleSkeleton BattleSkeletonLoader(string strFullFileName, bool isLimitBreak, bool loadGeometryQ)
        {
            var bSkeleton = new BattleSkeleton(strFullFileName, isLimitBreak, loadGeometryQ, false);
            BattleSkeletonPolyCheck(ref bSkeleton);
            return bSkeleton;
        }

        private static void BattleSkeletonPolyCheck(ref BattleSkeleton bSkeleton)
        {
            if (!bDontCheckRepairPolys)
            {
                PModel? result = bSkeleton.CheckPolys();
                if (result != null)
                {
                    var p = (PModel)result;
                    if (RepairPolysPrompt(ref p))
                    {
                        bSkeleton.RepairPolys();
                    }
                }
            }
        }

        public static void PModelLoader(ref PModel fPModel, string strPFolder, string strPFileName, bool bComputeNormals)
        {
            bool repairPolys = false;
            if (!bDontCheckRepairPolys)
            {
                repairPolys = RepairPolysPrompt(ref fPModel);
            }
            bool removeTextureCoords = TextureCoordinateCheck(ref fPModel);
            LoadPModel(ref fPModel, strPFolder, strPFileName, bComputeNormals, repairPolys, removeTextureCoords);
        }

        public static bool TextureCoordinateCheck(ref PModel pModel)
        {
            for (int i = 0; i < pModel.Groups.Length; ++i)
            {
                if (TextureCoordinateCheck(ref pModel, i))
                    return true;
            }
            return false;
        }

        public static bool TextureCoordinateCheck(ref PModel pModel, int index)
        {
            if (pModel.Groups[index].texID > 0)
                return RemoveTextureCoordsPrompt(index);
            return false;
        }
    }
}
