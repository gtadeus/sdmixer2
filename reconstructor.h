#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

#include <QString>
#include <tiffio.h>
#include <boost/iostreams/device/mapped_file.hpp>
#include "sdmixer.h"
#include "pairfinder.h"

class Reconstructor
{
public:
    struct Coordinates
    {
        int x;
        int y;

        int z;

        int get(int dim)
        {
            if (dim == 0)
                return x;
            if (dim == 1)
                return y;
            if (dim == 2)
                return z;
            return 0;
        }

    };

    Reconstructor(sdmixer *s);
    void run();
    void convertTo2D();
    void createKernel();
    void convolve(double *kernel);
    void hist_correct();
    void getIndexFromXYZ();
    void XYZfromFilter(
            std::vector<PairFinder::Localization>& data);

    void outputTIFF(QString path);
    void setMinMax(double min_x, double max_x,
                   double min_y, double max_y,
              double min_z, double max_z);


    uint64_t linearIndex(Coordinates c);
    void getMinMax();
    void setArray();

private:
    sdmixer *sdm;

    std::vector<Coordinates> xyz;

    //uint16 page;
    float xres = 100;
    float yres = 100;

    TIFF *out;

    uint64_t maxPixels=1;

    int xy_binning=10;
    int z_binning=10;

    int hist_correct_value;
    int hist_threshold;

    double *kernel;
    int kernel_dim;

    int dimensions=3;

    int min_x, min_y, min_z;
    int max_x, max_y, max_z;

    int image_min[3]={0};
    int image_max[3]={0};

    double min_val[3]={0};
    double max_val[3]={0};

    char * tmpfile = "tmp.file";



    //uint8* data;

};

#endif // RECONSTRUCTOR_H
