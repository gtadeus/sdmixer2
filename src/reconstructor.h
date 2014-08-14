#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

#include "sdmixer.h"

#include <tiffio.h>
#include <boost/iostreams/device/mapped_file.hpp>
#include <fftw3.h>


#include <QThread>
#include <QEventLoop>
#include <QFile>
#include <QFileInfo>
#include <QString>

#include <QtXml/QtXml>
#include <QtXml/QDomDocument>



class Reconstructor : public QObject
{
    Q_OBJECT
public:

    void Convolution2();
    void CreateGaussianKernel();
    void setKernel();

    struct Coordinates
    {
        int x;
        int y;
        int z=0;

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


    struct Kernel {
        int k=1;
        float sigma_xy=2;
        float sigma_z=2;
        int size=2*k+1;
        int size_z=1;
        float* data;
        bool make3D=true;
    };

    struct ConvPixel{
        double sum;
        int index;
    };

    Reconstructor(sdmixer *s,
                  std::vector<sdmixer::Localization> *data,
                  int current_filter);
    Reconstructor(sdmixer *s, QString xyz_file);
    void run();
    void createKernel();
    bool file_exists(char *filename);

    void hist_correct();
    void getIndexFromXYZ();
    void XYZfromFilter();

    void map8bit();
    void setOutputPath();
    void outputTIFF();
    void setMinMax(sdmixer::min_max m);


    uint64_t linearIndex(Coordinates c);
    uint64_t linearIndex3DFFT( int i, int j, int k, int z_fftw, int h_fftw);
    uint64_t linearIndex2DFFT( int i, int j, int h_fftw);
    uint64_t linearIndexTIFF( int i, int j, int k, int img_sizeX, int img_sizeY);

    void getMinMax();
    void setArray();

    void getHeader(QString header);
    void initData(int current_filter);
    void getSettingsFromGUI();
    void doWorkNow();

    void hist_eq();

signals:
    void finished();
    //void started();

public slots:
    void doWork();

private:

    sdmixer *sdm;
    Kernel krn;

    QString xyz_file_parameter;

    bool FilterInput=false;

    sdmixer::gaussian_kernel globalKernel;
    std::vector<sdmixer::gaussian_kernel> vecKernel;
    std::vector<Coordinates> xyz;
    std::vector<sdmixer::Localization> *input_data;

    sdmixer::min_max min_maxValues;

    int rawDataCols=0;

    bool runConvolution;
    bool oneConvolutionKernel;
    bool force2D;

    int NM_PER_PX;

    float xres = 100;
    float yres = 100;

    TIFF *out;

    double hist_correct_value;
    double hist_threshold;
    bool perform_hist_eq;
    bool sqrtCum;

    int xCol=0;
    int yCol=1;
    int zCol=2;

    int dimensions=0;

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
    char *convolved_image = "conv_img.tmp";


    int dbl_image_min=0;
    int dbl_image_max=0;

    QString tiff_out_suffix=QString("out");
    QString tiff_out_file;
    QString input_base_name;
    QString input_file;
    QString output_dir;

    int curr_filter=1;
    bool minMaxDefined=false;
    bool dimDefined=false;
    int max_filter=1;

    bool xyzFile = false;

    bool doWorkLater=false;



};

#endif // RECONSTRUCTOR_H
