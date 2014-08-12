#ifndef FILTER_H
#define FILTER_H

#include <QApplication>
#include <QImage>
#include <QColor>
#include <vector>
#include <fstream>
#include "sdmixer.h"
#include "pairfinder.h"
#include <QThread>

#include <QtXml/QtXml>
#include <QtXml/QDomDocument>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

//class sdmixer;

class Filter : public QObject
{
    Q_OBJECT
public:
    Filter(sdmixer *s, QString file);
    Filter(sdmixer *s, std::vector<sdmixer::Localization> *data);

    int getMaxIntLong() { return maxIntLong; }
    void setMaxIntLong(int val) { this->maxIntLong=val; }
    void setMaxIntShort(int val) { this->maxIntShort=val; }


    void init();
    void loadFile(QString str);

signals:
    void finished();

public slots:
    void doWork();

private:

    sdmixer *sdm;
    int maxIntLong;
    int maxIntShort;

    QString fileName, output_dir, outputFile, intensitySpaceFile;
    QString FilterSuffix = "_filter_out.txt";
    QString IntensitySpaceSuffix = "_IntensitySpace.png";
    int dimensions=0;
    int rawDataCols=0;

    double maxIntLongFromFile;
    double maxIntShortFromFile;

    std::vector<sdmixer::Localization> *input;

    //std::vector<double> pairs_out_input;

    double precision = 0.1;
    int nr_of_filters = 1;


    QString pathIntensitySpace;

    std::vector<QImage> QImgFilters;
    std::vector<QString> FilterInputFiles;
    QString output_directory;

    sdmixer::min_max min_maxValues;

};

#endif // FILTER_H
