using System.Drawing;
using System.Drawing.Drawing2D;
using System.Windows.Forms;

namespace KimeraCS.Core
{
    public static class Utils_WinForms
    {
        public static bool FindWindowOpened(string strWindowName)
        {
            foreach (Form itmFrm in Application.OpenForms)
            {
                if (itmFrm.Name == strWindowName) return true;
            }

            return false;
        }

        //  ------------------------------------------------------------------------------------------------
        //  ======================================== BITMAP HELPERS ========================================
        //  ------------------------------------------------------------------------------------------------

        /// <summary>
        /// Fits a bitmap to a PictureBox, scaling with NearestNeighbor interpolation.
        /// </summary>
        public static Bitmap FitBitmapToPictureBox(PictureBox pbIn, int iImgWidth, int iImgHeight, Bitmap srcBitmap)
        {
            if (srcBitmap == null)
                return null;

            Bitmap tmpBMP;
            float fAspectRatio = (float)iImgWidth / (float)iImgHeight;

            // Get the size available
            float fWidth = pbIn.ClientSize.Width;
            float fHeight = pbIn.ClientSize.Height;

            // Adjust the wid/hgt ratio to match aspect_src
            if (fWidth / fHeight > fAspectRatio)
            {
                // The area is too short and wide. Make it narrower.
                fWidth = fAspectRatio * fHeight;
            }
            else
            {
                // The area is too tall and thin. Make it shorter.
                fHeight = fWidth / fAspectRatio;
            }

            // Create image at the correct size.
            if (iImgWidth < pbIn.ClientSize.Width && iImgHeight < pbIn.ClientSize.Height)
                tmpBMP = new Bitmap(pbIn.ClientSize.Width, pbIn.ClientSize.Height, System.Drawing.Imaging.PixelFormat.Format32bppArgb);
            else
                tmpBMP = new Bitmap(iImgWidth, iImgHeight, System.Drawing.Imaging.PixelFormat.Format32bppArgb);

            using (Graphics g = Graphics.FromImage(tmpBMP))
            {
                g.InterpolationMode = InterpolationMode.NearestNeighbor;
                g.PixelOffsetMode = PixelOffsetMode.Half;
                g.DrawImage(srcBitmap, 0, 0, fWidth, fHeight);
            }

            return tmpBMP;
        }
    }
}
