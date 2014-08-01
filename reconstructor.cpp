#include "reconstructor.h"
#include "pairfinder.h"
unsigned char Reconstructor::array[XSIZE*YSIZE] = {0};

void Reconstructor::XYZfromFilter(std::vector<PairFinder::Localization> &data)
{
    for (auto i : data)
    {
        Coordinates c;
        if ( i.xShort > min_x)
            if ( i.xShort < max_x)
                if ( i.yShort > min_y)
                    if (i.yShort < max_y)
                        if( i.zShort > min_z)
                            if ( i.zShort < max_z)
                            {
                                c.x = round(i.xShort/xy_binning);
                                c.y = round(i.yShort/xy_binning);
                                c.z = round(i.zShort/z_binning);
                            }
        xyz.push_back(c);

        if ( i.xLong > min_x)
            if ( i.xLong < max_x)
                if ( i.yLong > min_y)
                    if (i.yLong < max_y)
                        if( i.zShort > min_z)
                            if ( i.zShort < max_z)
                            {
                                c.x = round(i.xShort/xy_binning);
                                c.y = round(i.yShort/xy_binning);
                                c.z = round(i.zShort/z_binning);
                            }

        xyz.push_back(c);
    }
}

void Reconstructor::getIndexFromXYZ()
{
    // N1xN2xN3 Array
    Nl[3]={round(max_x), round(max_y), round(max_z)};
    index = 0;
    std::vector<int> image;
    for ( auto i : xyz)
    {
        for (int k = 0; k < dimensions; ++k)
        {
            for(int l = k; k < dimensions; ++l)
            {
                index*=Nl[k];
            }
            if(k==0)
                index*=i.x;
            if(k==1)
                index*=i.y;
            if(k==2)
                index*=i.z;

            image.push_back(int(index));

        }
    }
}

Reconstructor::Reconstructor()
{
    //load filter output
    //convert xyz coordinates to
    //round vectors to full nm
    // get min max nm

    //horzcat input

    // CONVERT DATA TO XYZ
    // dann kann für jeden kanal einzelm geblurrt werden


    // image pixel

    // extract each

    // blur for each channel

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
    //for ( int i = 0; i < max_z; ++i)

    uint16 spp, bpp, photo, res_unit;

    out = TIFFOpen(path.toLocal8Bit(), "w");
    if (!out)
    {
        return;
         // "Can't open %s for writing\n",
    }

    spp = 1; /* Samples per pixel */
    //3 == rgb, 4 == alpha channel
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
            TIFFWriteScanline(out, &array[j * image_width], j, 0);

        TIFFWriteDirectory(out);
    }

     TIFFClose(out);

}
