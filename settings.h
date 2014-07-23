#ifndef SETTINGS_H
#define SETTINGS_H

#include <vector>
#include "sdmixer.h"

class Settings
{
public:
    Settings();


    void getSettingsFromUI();

    //General
    void setRunPairfinder(bool run){this->runPaiFinder=run;}
    void setRunFilter(bool run);
    void setRunReconstructor(bool run);

    void setPixelSizeNM(int px){this->pixelSizeNM = px;}

    //PairFinder
    void setOffset(std::vector<double> offset);
    void setEpsilon(std::vector<double> epsilon);

    void setPerformFishing(bool val);
    void setFishingIncrement(int val);
    void setFishingRange(int val);
    void setFishingSubset(int val);

    //Filter
    void setMaxIntLong(int val);
    void setMaxIntShort(int val);
    void setPrecision(double val);

    //Reconstructor
    void setXYBinning(double xy);
    void setZBinning(double z);

private:
    bool runPaiFinder;
    bool runFilter;
    bool runReconstructor;

    int pixelSizeNM;

    std::vector<double> offset;
    std::vector<double> epsilon;

    bool performFishing;
    int fishingIncrement;
    int fishingRange;
    int fishingSubset;

    int maxIntensityLong;
    int maxIntensityShort;
    double precision;

    double xyBinning;
    double zBinning;

};

#endif // SETTINGS_H
