#ifndef SDMIXER_H
#define SDMIXER_H

//#define DEFAULT_SETTINGS "default_settings.txt"

#define SDMIXER_VERSION  2.01



#include <QStandardPaths>

#include <QMainWindow>
#include <QListWidget>
#include <QDebug>
#include <string>
#include <vector>
#include <QString>
#include <QTextEdit>
#include <QTextStream>
#include <ctime>
#include <sstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <fcntl.h>
/*#include "pairfinder.h"
#include "reconstructor.h"
class Reconstructor;*/

class Settings;
class PairFinder;


namespace Ui {
class sdmixer;
}

class sdmixer : public QMainWindow
{
    Q_OBJECT

public:
    struct fishing{
        bool run=false;
        double increment=0;
        int range=0;
        int subset=0;
        QString unit;
    };
    struct log {
        log(sdmixer *s) {
            this->sdm = s;
        }
        log(sdmixer *s, std::string os) {
            this->sdm = s;
            this->msg = os;
        }

        ~log() {
            /*QString qstr = sdm->timestamp();
            qstr.append(" : ");
            qstr.append(QString::fromStdString(m_SS.str()));

            sdm->writeToConsole(qstr);*/
            QString qstr = sdm->timestamp();
            qstr.append(" : ");
            qstr.append(QString::fromStdString(msg));
            sdm->writeToConsole(qstr);
            qDebug() << qstr;
        }

    public:
        // accepts just about anything
        template<class T>
        log &operator<<(const T &x) {
            m_SS << x;
            return *this;
        }

    private:
        sdmixer *sdm;
        std::ostringstream m_SS;
        std::string msg;
    };
    struct Localization {

        double xLong, yLong, zLong;
        double xShort, yShort, zShort;
        double LongIntensity, ShortIntensity;
        int frame;

        int filter = 0;

        double getShortDim(int dim)
        {
            if ( dim == 0)
                return xShort;
            if ( dim == 1)
                return yShort;
            if ( dim == 2)
                return zShort;
            return 0;
        }
        double getLongDim(int dim)
        {
            if ( dim == 0)
                return xLong;
            if ( dim == 1)
                return yLong;
            if ( dim == 2)
                return zLong;
            return 0;
        }
        void setShortDim(int dim, double val)
        {
            if ( dim == 0)
                xShort = val;
            if ( dim == 1)
                yShort = val;
            if ( dim == 2)
                zShort = val;
        }
        void setLongDim(int dim, double val)
        {
            if ( dim == 0)
                xLong = val;
            if ( dim == 1)
                yLong = val;
            if ( dim == 2)
                zLong = val;
        }
    };
    struct min_max {
        double min_x=0, max_x=0;
        double min_y=0, max_y=0;
        double min_z=0, max_z=0;

        double getMin(int dim)
        {
            if ( dim == 0)
                return min_x;
            if ( dim == 1)
                return min_y;
            if ( dim == 2)
                return min_z;
            return 0;
        }
        double getMax(int dim)
        {
            if ( dim == 0)
                return max_x;
            if ( dim == 1)
                return max_y;
            if ( dim == 2)
                return max_z;
            return 0;
        }

    };
    struct Columns{
        int xShort;
        int yShort;
        int zShort;
        int ShortAmp;
        int frame;
        int xLong;
        int yLong;
        int zLong;
        int LongAmp;
        int filter;
        int dimensions=0;

        int x, y, z;
        int Amplitude;
        int rawDataCols=0;

        int getShortCol(int dim)
        {
            if ( dim == 0)
                return xShort;
            if ( dim == 1)
                return yShort;
            if ( dim == 2)
                return zShort;
            return 0;
        }
        int getLongCol(int dim)
        {
            if ( dim == 0)
                return xLong;
            if ( dim == 1)
                return yLong;
            if ( dim == 2)
                return zLong;
            return 0;
        }
        int getXYZCol(int dim)
        {
            if ( dim == 0)
                return x;
            if ( dim == 1)
                return y;
            if ( dim == 2)
                return z;
            return 0;
        }

    };
    enum input_file_t{XYZ_FILE, PAIRS_FILE, FILTER_FILE};



    struct offset_units{
        QString xOffset;
        QString yOffset;
        QString zOffset;
        QString xEpsilon;
        QString yEpsilon;
        QString zEpsilon;
        QString getOffset(int dim){
            if(dim==0)
                return xOffset;
            if(dim==1)
                return yOffset;
            if(dim==2)
                return zOffset;
        }
        QString getEpsilon(int dim){
            if(dim==0)
                return xEpsilon;
            if(dim==1)
                return yEpsilon;
            if(dim==2)
                return zEpsilon;
        }
    };

    struct gaussian_kernel{
        double FWHM_xy=0;
        double FWHM_z=0;
        QString unitFWHM_xy="nm";
        QString unitFWHM_z="nm";
        QString filterName;
    };


    explicit sdmixer(QWidget *parent = 0);
    ~sdmixer();

    // multithread
    void nextStage();
    QString getCurrentFile(){ return InputFiles[current_file]; }

    // logging:
    QString timestamp();
    QString timestampLogFile();
    QString error_msg(QString msg);
    void writeToConsole(QString q);

    // helper functions to manipulate the GUI
    void insertItem(QString filename, QListWidget *list);
    void setStartDemixingButtonEnabled(bool val);
    void setCancelRunButtonEnabled(bool val);
    void UpdateShortPositionComboBox(QString str);
    void getColsAndRows(QString file);

    // save & load Settings
    bool getSettingsFromUI();
    void setSettingsToUI(Settings s);



    // Session & Pairfinder
    std::vector<QString> getInputFiles(){return this->InputFiles;}
    bool getRunPairfinder(){return this->runPairFinder;}
    bool getRunFilter(){return this->runFilter;}
    bool getRunReconstructor(){return this->runReconstructor;}
    bool getForce2D(){return this->force2D;}
    QString getOutputDirectory(){return this->output_directory;}

    int getPixelSize(){return this->pixelSizeNM;}
    double getOffset(int dim){
        //qDebug() << "sdmixer : offset " << dim << ": " <<offset[dim];
        return offset[dim];
    }
    double getEpsilon(int dim){
        //qDebug()<<"sdmixer: epsilon " << dim << ": " << epsilon[dim];
        return epsilon[dim];}
    fishing getFishing(){return this->fishingSettings;}
    offset_units getOffsetUnits(){return this->offsetUnits;}
    QString getCameraOrientation() { return this->CameraOrientation;}
    QString getShortChannelPosition() { return this->ShortChannelPosition;}
    bool getRunGrouping() { return runGrouping; }
    double getGroupingRadius() { return groupingRadius; }
    QString getGroupingUnits() { return groupingUnits; }

    // Filter
    std::vector<QString> getFilterFiles(){return this->FilterFiles;}
    double getMaxIntShort(){return this->maxIntensityShort;}
    double getMaxIntLong(){return this->maxIntensityLong;}
    double getPrecision(){return this->precision;}
    QString getFilterOrientation(){return this->FilterOrientation;}
    bool getPlotIntensitySpace() { return this->plotIntensitySpace;}

    // Reconstructor
    double getReconstructor_xyBinning(){return this->xyBinning;}
    double getReconstructor_zBinning(){return this->zBinning;}
    bool getRunConvolution() { return this->runConvolution; }
    bool getOneKernelForAllChannels() { return oneKernelForAllChannels;}
    gaussian_kernel getGlobalKernel() {return globalKernel; }
    std::vector<gaussian_kernel> getConvolutionKernel(){ return vec_kernel;}
    bool getNonLinearHistEq(){return this->nonLinearHistogramEqual;}
    double getHisteqCoefficient() { return histeqCoefficient;}
    double getThreshold() { return Threshold;}
    bool getSqrtCummulation() { return sqrtCummulation;}
    bool getLZWCompression() { return LZWCompression;}
    bool getResliceZ() { return ResliceZ; }
    int getStartRescliceZ(){ return startRescliceZ; }
    int getEndRescliceZ() { return endRescliceZ;}

    bool getPerformNNStatistic() { return performNNStatistic; }

    QString getTiffTempFile(){ return tiff_temp_file; }
    QString getConvImgTempFile(){ return conv_image_temp_file; }

    std::vector<Localization> * getPfOutput(){ return pf_output; }


    // save Pairfinder Data
    void pushBackLocalization(Localization l)
    {
        pf_output->push_back(l);
    }
    void clearLocalizations()
    {
        pf_output->clear();
        //std::vector<Localization>()->swap(pf_output);
    }

    // save min/max from header from pairfinder run for later
    void setPF_min_maxValues(min_max m){MinMaxValues=m;}
    min_max getPF_min_maxValues(){return MinMaxValues;}

    //void pushBackKernel(gaussian_kernel gk){ vec_kernel.push_back(gk); }

    //int getCurrentDimensions() { return current_dimensions; }

    void runStage(int stage);
    void nextStage2();

    int getCurrentDimensions(QString file);

    void writeToLogFile(QString msg, QString msg2="");

    void getHeader(QString header,
                   Columns &columns,
                   min_max &min_maxValues,
                   input_file_t &INPUT_FILE);

    void writeHeader(QTextStream &out,
                     int dimensions,
                     min_max &min_maxValues,
                     input_file_t &INPUT_FILE);


protected:
    void dragEnterEvent(QDragEnterEvent* event);
    void dropEvent(QDropEvent* event);
private slots:
    // Quit Application
    void on_actionQuit_triggered();
    // File Lists QListWidget Add Items
    void on_addFileButton_clicked();
    void on_addFilterButton_clicked();
    void on_removeFilesButton_clicked();
    void on_removeFilterButton_clicked();

    //start & cancel Demixing
    void on_startDemixing_clicked();
    void on_pushButton_CancelRun_clicked();

    //About
    void on_actionAbout_sdmixer_triggered();


    // Save & Load Settings
    void on_actionSave_Settings_as_Default_triggered();
    void on_actionLoad_Settings_triggered();
    void on_actionSave_Preferences_triggered();

    // SelectOutputDirectory
    void on_pushButton_selectOutputDirectory_clicked();


    // this signal is triggered by Pairfinder/Filter/Reconstructor
    // when respective Thread has finished
    void threadReady();


    void InputFileClicked(QListWidgetItem *item);
    void CameraOrientationChanged(QString str);
    void filterInfo(QString str);
    void loadKernel(int a);
    void SameKernelCheckBoxChanged(int state);
    void KernelChannelHighlighted(QString str);

    void HistEqCheckBoxClicked(int state);


    void on_pushButton_saveChannel_clicked();



private:

    QString DEFAULT_SETTINGS;

    std::vector<gaussian_kernel> vec_kernel;
    std::vector<Localization> *pf_output;

    int current_dimensions=0;
    static const int max_dims = 3;
    // Settings
    // Session & Pairfinder
    std::vector<QString> InputFiles;
    QString output_directory="";
    bool runPairFinder=false;
    bool runFilter=false;
    bool runReconstructor=false;
    bool force2D=false;
    int pixelSizeNM=0;

    double offset[max_dims]={0};
    double epsilon[max_dims]={0};

    fishing fishingSettings;
    QString CameraOrientation;
    QString ShortChannelPosition;
    offset_units offsetUnits;

    bool runGrouping;
    double groupingRadius;
    QString groupingUnits;

    std::vector<QString> FilterFiles;
    int maxIntensityLong=0;
    int maxIntensityShort=0;
    double precision=0;
    QString FilterOrientation;
    bool plotIntensitySpace;

    double xyBinning=0;
    double zBinning=0;
    bool runConvolution;
    bool oneKernelForAllChannels;
    gaussian_kernel globalKernel;
    bool nonLinearHistogramEqual;
    double histeqCoefficient;
    double Threshold;
    bool sqrtCummulation;

    bool LZWCompression;
    bool ResliceZ;
    int startRescliceZ;
    int endRescliceZ;

    bool performNNStatistic;


    Ui::sdmixer *ui;
    QListWidget *listWidgetInputFiles;
    QListWidget *listWidgetFilters;
    QTextEdit *console;

    min_max MinMaxValues;


    int rawDataCols, rawDataRows;
    int current_stage=0;
    int current_filter=1;
    int filter_max=0;
    std::vector<QString>::size_type current_file = 0;

    std::ofstream logFile;

    QString tiff_temp_file;
    QString conv_image_temp_file;


};


#endif // SDMIXER_H
