#ifndef SDMIXER_H
#define SDMIXER_H

#define DEFAULT_SETTINGS "default_settings.txt"

#include <QMainWindow>
#include <QListWidget>
#include <QDebug>
#include <string>
#include <vector>
#include <QString>
#include <QTextEdit>
#include <QTextStream>
#include <ctime>
#include "pairfinder.h"
#include "reconstructor.h"
class Reconstructor;

class Settings;

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

    enum camera_orientation_t{
        LEFT_RIGHT,
        TOP_BOTTOM
    };
    enum short_channel_position_t{
        TOP,
        BOTTOM,
        LEFT,
        RIGHT
    };
    enum filter_orientation_t{
        X_SHORT,
        Y_SHORT
    };

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


    explicit sdmixer(QWidget *parent = 0);
    ~sdmixer();

    // multithread
    void nextStage();

    // logging:
    QString timestamp();
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
        qDebug() << "offset " << dim << ": " <<offset[dim];
        return offset[dim];
    }
    double getEpsilon(int dim){
        qDebug()<<" epsilon " << dim << ": " << epsilon[dim];
        return epsilon[dim];}
    sdmixer::fishing getFishing(){return this->fishingSettings;}
    offset_units getOffsetUnits(){return this->offsetUnits;}
    QString getCameraOrientation() { return this->CameraOrientation;}
    QString getShortChannelPosition() { return this->ShortChannelPosition;}

    // Filter
    std::vector<QString> getFilterFiles(){return this->FilterFiles;}
    double getMaxIntShort(){return this->maxIntensityShort;}
    double getMaxIntLong(){return this->maxIntensityLong;}
    double getPrecision(){return this->precision;}
    QString getFilterOrientation(){return this->FilterOrientation;}

    // Reconstructor
    double getReconstructor_xyBinning(){return this->xyBinning;}
    double getReconstructor_zBinning(){return this->zBinning;}
    bool getRunConvolution() { return this->runConvolution; }
    bool getNonLinearHistEq(){return this->nonLinearHistogramEqual;}
    double getHisteqCoefficient() { return histeqCoefficient;}
    double getThreshold() { return Threshold;}
    bool getSqrtCummulation() { return sqrtCummulation;}
    bool getLZWCompression() { return LZWCompression;}
    bool getResliceZ() { return ResliceZ; }
    int getStartRescliceZ(){ return startRescliceZ; }
    int getEndRescliceZ() { return endRescliceZ;}


    // save Pairfinder Data
    void pushBackLocalization(PairFinder::Localization l)
    {
        pf_output.push_back(l);
    }
    // save min/max from header from pairfinder run for later
    void setPF_min_maxValues(PairFinder::min_max m){MinMaxValues=m;}
    PairFinder::min_max getPF_min_maxValues(){return MinMaxValues;}




    std::map<camera_orientation_t, QString> cam_orient;
    std::map<short_channel_position_t, QString> short_pos;
    std::map<filter_orientation_t, QString> filter_orient;




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


private:
    std::vector<PairFinder::Localization> pf_output;

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

    std::vector<QString> FilterFiles;
    int maxIntensityLong=0;
    int maxIntensityShort=0;
    double precision=0;
    QString FilterOrientation;

    double xyBinning=0;
    double zBinning=0;
    bool runConvolution;
    bool nonLinearHistogramEqual;
    double histeqCoefficient;
    double Threshold;
    bool sqrtCummulation;

    bool LZWCompression;
    bool ResliceZ;
    int startRescliceZ;
    int endRescliceZ;




    Ui::sdmixer *ui;
    QListWidget *listWidgetInputFiles;
    QListWidget *listWidgetFilters;
    QTextEdit *console;

    PairFinder::min_max MinMaxValues;


    int rawDataCols, rawDataRows;
    int current_stage=0;
    std::vector<QString>::size_type current_file = 0;




};


#endif // SDMIXER_H
