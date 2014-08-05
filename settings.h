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

    //Reconstructor
    void setXYBinning(double val){this->xyBinning=val;}
    void setZBinning(double val){this->zBinning=val;}


    int getPixelSizeNM(){return pixelSizeNM;}
    bool getRunPairFinder(){return runPairFinder;}
    bool getRunFilter(){return runFilter;}
    bool getRunReconstructor(){return runPairFinder;}
    bool getForce2D(){return force2D;}

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

    sdmixer::fishing getFishing(){return fishing;}
    int getMaxIntLong(){return maxIntensityLong;}
    int getMaxIntShort(){return maxIntensityShort;}
    double getPrecision(){return precision;}

    double getXYbinning(){return xyBinning;}
    double getZbinning(){return zBinning;}

private:

    bool runPairFinder=true;
    bool runFilter=true;
    bool runReconstructor=true;
    bool force2D=true;

    int pixelSizeNM=0;

    static const int max_dims = 3;
    double offset[max_dims]={0};
    double epsilon[max_dims]={0};

    sdmixer::fishing fishing;

    int maxIntensityLong=0;
    int maxIntensityShort=0;
    double precision=0.0;

    double xyBinning=0;
    double zBinning=0;

    std::vector<QString> InputFiles;
    std::vector<QString> FilterFiles;

    QDomDocument settingsFile;

};

#endif // SETTINGS_H
