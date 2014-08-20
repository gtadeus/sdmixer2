#ifndef SETTINGS_H
#define SETTINGS_H

#include <vector>
#include <QtXml/QtXml>
#include <QtXml/QDomDocument>

#include "sdmixer.h"

class Settings
{
public:
    Settings(sdmixer *s);   // open from sdmixer
    Settings(QString file); // open file
    void writeSettingsToFile(QString fileName);
    void loadFromFile(QString file);
    void initXML();
    void retrieveElements(QDomElement root, QString tag, QString att);
    void retrieveElements2(QDomNodeList nodes);
    QDomElement appendChildNode(QDomElement parent, QString name, double val);
    QDomElement appendChildNode(QDomElement parent, QString name, QString str);
    QDomElement createField(QString name);

    //General
    void setInputFiles(std::vector<QString> str);
    void setRunPairfinder(bool run){this->runPairFinder=run;}
    void setRunFilter(bool run){this->runFilter=run;}
    void setRunReconstructor(bool run){this->runReconstructor=run;}
    void setForce2D(bool v){this->force2D=v;}
    void setPixelSizeNM(int px){this->pixelSizeNM = px;}
    //PairFinder
    void setOffset(int dim, double val){offset[dim]=val;}
    void setEpsilon(int dim, double val){epsilon[dim]= val;}
    void setFishing(sdmixer::fishing f){this->fishing=f;}
    //Filter
    void setFilterFiles(std::vector<QString> s){this->FilterFiles=s;}
    void setMaxIntLong(int val){this->maxIntensityLong=val;}
    void setMaxIntShort(int val){this->maxIntensityShort=val;}
    void setPrecision(double val){this->precision=val;}
    void setPlotIntensitySpace(bool val){this->plotIntensitySpace=val;}

    //Reconstructor
    void setXYBinning(double val){this->xyBinning=val;}
    void setZBinning(double val){this->zBinning=val;}
    void setRunConvolution(bool val){this->runConvolution=val;}



    std::vector<QString> getInputFiles() {return InputFiles;}
    int getPixelSizeNM(){return pixelSizeNM;}
    bool getRunPairFinder(){return runPairFinder;}
    bool getRunFilter(){return runFilter;}
    bool getRunReconstructor(){return runReconstructor;}
    bool getForce2D(){return force2D;}
    QString getOutputDirectory(){return output_directory;}

    double getOffset(int dim){
        if(dim < max_dims)
            return offset[dim];
        else
            return 0;
    }
    double getEpsilon(int dim){
        if(dim < max_dims)
            return epsilon[dim];
        else
            return 0;
    }
    sdmixer::offset_units getOffsetUnits(){return offsetUnits;}
    QString getCameraOrientation(){return CameraOrientation;}
    QString getShortChannelPosition(){return ShortChannelPosition;}
    sdmixer::fishing getFishing(){return fishing;}

    std::vector<QString> getFilterFiles() {return FilterFiles;}
    int getMaxIntLong(){return maxIntensityLong;}
    int getMaxIntShort(){return maxIntensityShort;}
    double getPrecision(){return precision;}
    QString getFilterOrientation(){return FilterOrientation; }
    bool getPlotIntensitySpace() {return plotIntensitySpace; }

    double getXYbinning(){return xyBinning;}
    double getZbinning(){return zBinning;}
    bool getRunConvolution() { return runConvolution;}
    bool getOneKernelForAllChannels() { return oneKernelForAllChannels; }
    sdmixer::gaussian_kernel getGlobalKernel() {return global_kernel;}
    std::vector<sdmixer::gaussian_kernel> getConvolutionKernel() { return vec_kernel;}
    bool getNonLinearHistEq() { return nonLinearHistogramEqual;}
    double getCorrectionCoefficient() {return histeqCoefficient;}
    double getThreshold() { return Threshold;}
    bool getSqrtCummulation(){ return sqrtCummulation;}
    bool getLZWCompression() { return LZWCompression;}
    bool getResliceZ(){return ResliceZ;}
    int getStartSliceZ(){return startRescliceZ;}
    int getEndSliceZ() { return endRescliceZ;}


private:
    static const int max_dims = 3;
    // Session & Pairfinder
    std::vector<QString> InputFiles;
    QString output_directory;
    bool runPairFinder=true;
    bool runFilter=true;
    bool runReconstructor=true;
    bool force2D=true;
    int pixelSizeNM=0;

    double offset[max_dims]={0};
    double epsilon[max_dims]={0};
    sdmixer::fishing fishing;
    QString CameraOrientation;
    QString ShortChannelPosition;
    sdmixer::offset_units offsetUnits;

    // FIlter
    std::vector<QString> FilterFiles;
    int maxIntensityLong=0;
    int maxIntensityShort=0;
    double precision=0.0;
    QString FilterOrientation;
    bool plotIntensitySpace;

    // Reconstructor
    double xyBinning=0;
    double zBinning=0;
    bool runConvolution;
    bool oneKernelForAllChannels;
    sdmixer::gaussian_kernel global_kernel;
    std::vector<sdmixer::gaussian_kernel> vec_kernel;
    bool nonLinearHistogramEqual;
    double histeqCoefficient;
    double Threshold;
    bool sqrtCummulation;

    bool LZWCompression;
    bool ResliceZ;
    int startRescliceZ;
    int endRescliceZ;



    QDomDocument settingsFile;

};

#endif // SETTINGS_H
