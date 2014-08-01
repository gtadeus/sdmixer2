#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

#include <QString>
#include <tiffio.h>

class Reconstructor
{
public:
    struct Coordinates
    {
        int x;
        int y;
        int z;
    };

    Reconstructor();
    void run();
    void convertTo2D();
    void createKernel();
    void convolve(double *kernel);
    void hist_correct();
    void getIndexFromXYZ();
    void XYZfromFilter(
            std::vector<PairFinder::Localization> &data);

    void outputTIFF(QString path);


private:
    std::vector<Coordinates> xyz;
    uint32 image_width, image_height;
    uint16 page;
    float xres = 100;
    float yres = 100;

    TIFF *out;



    int xy_binning;
    int z_binning;

    int hist_correct_value;
    int hist_threshold;

    double *kernel;
    int kernel_dim;

    int dimensions;

    int min_x, min_y, min_z;
    int max_x, max_y, max_z;

    char * image;

};

#endif // RECONSTRUCTOR_H
