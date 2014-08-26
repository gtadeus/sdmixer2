#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

#include "sdmixer.h"


#include <tiffio.h>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
//#include <fftw3.h>


#include <QThread>
#include <QEventLoop>
#include <QFile>
#include <QFileInfo>
#include <QString>

#include <QtXml/QtXml>
#include <QtXml/QDomDocument>

#include "nr3.h"
#include "pointbox.h"



class Reconstructor : public QObject
{
    Q_OBJECT
public:

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

    void Convolution(Kernel krn);
    void CreateGaussianKernel(Kernel &krn);
    void setKernel();

    Reconstructor(sdmixer *s,
                  std::vector<sdmixer::Localization> *data,
                  int current_filter);
    Reconstructor(sdmixer *s, QString xyz_file);
    void run();
    void createKernel();
    bool file_exists(const char *filename);

    std::vector<Kernel> all_kernels;

    void hist_correct();
    void getIndexFromXYZ();
    void XYZfromFilter();

    void setOutputPath();
    void outputTIFF();
    void setMinMax(sdmixer::min_max m);


    uint64_t linearIndex(Coordinates c);


    void getMinMax();
    void setArray();

    void getHeader(QString header);
    void initData(int current_filter);
    void getSettingsFromGUI();
    void doWorkNow();

    void hist_eq();
    void findNN(QString outputFile);

signals:
    void finished();
    //void started();

public slots:
    void doWork();

private:

    sdmixer *sdm;
    //Kernel krn;

    QString xyz_file_parameter;

    bool FilterInput=false;

    sdmixer::gaussian_kernel globalKernel;
    std::vector<sdmixer::gaussian_kernel> vecKernel;
    std::vector<Coordinates> xyz;
    vector< Point<3> > pts3D;
    vector< Point<2> > pts2D;
    //std::vector<dblCoordinates> xyz_not_rounded;
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

    sdmixer::Columns columns;
    sdmixer::input_file_t INPUT_FILE;

    int dimensions=0;

    int binning[3]={0};

    // min / max, first initialized from data
    // may then be overwritten by setMinMax with header data
    int image_min[3]={0};
    int image_max[3]={0};

    uint64_t image_size[3]={0};  // widhtx height x depth

    uint64_t maxPixels=1; // product for each i:dim   (image_size[i]+1)


    QString tiff_temp_file; // = "tiff_uint8.tmp";
    //char *tiff_uint8_file ="tiff_uint8.tmp";
    //char *fft_src_file = "fft_src.tmp";
    //char *fft_krnl_file = "fft_kernel.tmp";
    QString convolved_image; //= "conv_img.tmp";



    int dbl_image_min=0;
    int dbl_image_max=0;

    QString tiff_out_suffix=QString("out");
    QString tiff_out_file;
    QString input_base_name;
    QString input_file;
    QString output_dir;
    QString NN_output_file;

    int curr_filter=1;
    bool minMaxDefined=false;
    bool dimDefined=false;
    int max_filter=1;

    bool xyzFile = false;

    bool doWorkLater=false;

    bool ReslizeZ=false;
    int startResliceZ, endResliceZ;
    bool performNNStatistic;



};

#endif // RECONSTRUCTOR_H
