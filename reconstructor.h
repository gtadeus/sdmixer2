#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

#include <QString>
#include <tiffio.h>
#include <boost/iostreams/device/mapped_file.hpp>
//#include <math.h>
//#include <complex>
#include <fftw3.h>
#include "sdmixer.h"
//#include "pairfinder.h"
#include <QThread>
#include <QEventLoop>

class sdmixer;
class PairFinder;

class Reconstructor : public QObject
{
    Q_OBJECT
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
        void set(int dim, int val)
        {
            if (dim == 0)
                x=val;
            if (dim == 1)
                y=val;
            if (dim == 2)
                z=val;
        }

    };
    void Convolution();
    void CreateGaussianKernel();
    void setKernel();
    int pow2roundup (int x);


    struct Kernel {
        int k=1;
        float sigma_xy=2;
        float sigma_z=2;
        int size=2*k+1;
        int size_z=1;
        float* data;
        bool make3D=true;
    };

    Reconstructor(sdmixer *s,
                  std::vector<PairFinder::Localization> &PFinput_file,
                  QString filename);
    void run();
    void createKernel();

    void hist_correct();
    void getIndexFromXYZ();
    void XYZfromFilter();

    void map8bit();
    void setOutputPath();
    void outputTIFF();
    void setMinMax(PairFinder::min_max m);


    uint64_t linearIndex(Coordinates c);
    uint64_t linearIndex3DFFT( int i, int j, int k, int z_fftw, int h_fftw);
    uint64_t linearIndex2DFFT( int i, int j, int h_fftw);
    uint64_t linearIndexTIFF( int i, int j, int k, int img_sizeX, int img_sizeY);

    void getMinMax();
    void setArray();

signals:
    void finished();
    //void started();

public slots:
    void doWork();

private:

    sdmixer *sdm;
    Kernel krn;

    std::vector<Coordinates> xyz;


    float xres = 100;
    float yres = 100;

    TIFF *out;

    int hist_correct_value;
    int hist_threshold;


    int dimensions=3;

    int binning[3]={0};

    // min / max, first initialized from data
    // may then be overwritten by setMinMax with header data
    int image_min[3]={0};
    int image_max[3]={0};

    uint64_t image_size[3]={0};  // widhtx height x depth

    uint64_t maxPixels=1; // product for each i:dim   (image_size[i]+1)


    char *tiff_temp_file = "tiff_uint8.tmp";
    //char *tiff_uint8_file ="tiff_uint8.tmp";
    char *fft_src_file = "fft_src.tmp";
    char *fft_krnl_file = "fft_kernel.tmp";


    int dbl_image_min=0;
    int dbl_image_max=0;

    QString tiff_out_suffix=QString("out");
    QString tiff_out_file;
    QString input_base_name;
    QString input_file;
    QString output_dir;



};

#endif // RECONSTRUCTOR_H
