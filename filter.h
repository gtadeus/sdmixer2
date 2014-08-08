#ifndef FILTER_H
#define FILTER_H

#include <QApplication>
#include <QImage>
#include <QColor>
#include <vector>
#include <fstream>
#include "sdmixer.h"
//#include "pairfinder.h"
#include <QThread>

//class sdmixer;

class Filter : public QObject
{
    Q_OBJECT
public:
    Filter(sdmixer *s, QString file);
    Filter(sdmixer *s, std::vector<sdmixer::Localization> *data);
    void run();
    QImage loadFilterImage(QString path);
    void initializeIntensities();
    void roundIntensityValues();
    void drawIntensitySpace();



    int getMaxIntLong() { return maxIntLong; }
    void setMaxIntLong(int val) { this->maxIntLong=val; }
    void setMaxIntShort(int val) { this->maxIntShort=val; }

signals:
    void finished();

public slots:
    void doWork();

private:

    sdmixer *sdm;
    int maxIntLong;
    int maxIntShort;

    double maxIntLongFromFile;
    double maxIntShortFromFile;

    std::vector<sdmixer::Localization> *input;

    double precision = 0.1;
    int nr_of_filters = 1;


    QString pathIntensitySpace;

    std::vector<QImage> QImgFilters;
    std::vector<QString> FilterInputFiles;
    QString output_directory;

};

#endif // FILTER_H
