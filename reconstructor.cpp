#include "reconstructor.h"
unsigned char Reconstructor::array[XSIZE*YSIZE] = {0};

Reconstructor::Reconstructor()
{


}

void Reconstructor::run()
{
    // extract cols

    //output TIFF



    outputTIFF("out.tif");

}
void Reconstructor::convertTo2D()
{

}
void Reconstructor::createKernel()
{

}
void Reconstructor::convolve(double *kernel)
{

}
void Reconstructor::hist_correct()
{

}
void Reconstructor::outputTIFF( QString path)
{
    int i, j;

     for (j = 0; j < YSIZE; j++)
             for(i = 0; i < XSIZE; i++)
                     Reconstructor::array[j * XSIZE + i] = (unsigned char)(i * j);

    uint32 image_width, image_height;
    float xres, yres;
    uint16 spp, bpp, photo, res_unit;
    TIFF *out;
    uint16 page;

    out = TIFFOpen(path.toLocal8Bit(), "w");
    if (!out)
    {
        return;
         // "Can't open %s for writing\n",
    }
    image_width = XSIZE;
    image_height = YSIZE;

    spp = 1; /* Samples per pixel */
    bpp = 8; /* Bits per sample */
    photo = PHOTOMETRIC_MINISBLACK;

    for (page = 0; page < NPAGES; page++)
    {
        TIFFSetField(out, TIFFTAG_IMAGEWIDTH, image_width / spp);
        TIFFSetField(out, TIFFTAG_IMAGELENGTH, image_height);
        TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, bpp);
        TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, spp);
        TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(out, TIFFTAG_PHOTOMETRIC, photo);
        TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);
        /* It is good to set resolutions too (but it is not nesessary) */
        xres = yres = 100;
        res_unit = RESUNIT_INCH;
        TIFFSetField(out, TIFFTAG_XRESOLUTION, xres);
        TIFFSetField(out, TIFFTAG_YRESOLUTION, yres);
        TIFFSetField(out, TIFFTAG_RESOLUTIONUNIT, res_unit);

         /* We are writing single page of the multipage file */
        TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);

        /* Set the page number */
        TIFFSetField(out, TIFFTAG_PAGENUMBER, page, NPAGES);

        for (j = 0; j < image_height; j++)
            TIFFWriteScanline(out, &Reconstructor::array[j * image_width], j, 0);

        TIFFWriteDirectory(out);
    }

     TIFFClose(out);

}
