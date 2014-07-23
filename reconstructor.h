#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

#include <gsl/gsl_matrix.h>
#include <QString>
#include <tiffio.h>

class Reconstructor
{
public:
    Reconstructor();
    void run();
    void convertTo2D();
    void createKernel();
    void convolve(double *kernel);
    void hist_correct();

    void outputTIFF(QString path);


private:
    int xy_binning;
    int z_binning;

    int hist_correct_value;
    int hist_threshold;

    double *kernel;
    int kernel_dim;
    gsl_matrix *image;

    gsl_vector *x;
    gsl_vector *y;
    gsl_vector *z;

    const static int XSIZE = 20;
    const static int YSIZE = 20;
    const static int NPAGES = 2;
    // cant be like this...
    static unsigned char array[];





};

#endif // RECONSTRUCTOR_H
