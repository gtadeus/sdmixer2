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
    void WriteSettingsToFile(QString fileName);
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
    void setOffset(std::vector<double> val){
        for(int i=0; i < val.size(); ++i)
        {
         this->offset[i]=val[i];
        }
    }
    void setEpsilon(std::vector<double> val){
        for(int i=0; i < val.size(); ++i)
        {
         this->epsilon[i]= val[i];
        }
    }
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


private:

    bool runPairFinder;
    bool runFilter;
    bool runReconstructor;
    bool force2D;

    int pixelSizeNM;

    static const int max_dims = 3;
    double offset[max_dims]={0};
    double epsilon[max_dims]={0};

    sdmixer::fishing fishing;

    int maxIntensityLong;
    int maxIntensityShort;
    double precision;

    double xyBinning;
    double zBinning;

    std::vector<QString> InputFiles;
    std::vector<QString> FilterFiles;

    QDomDocument settingsFile;

};

#endif // SETTINGS_H
