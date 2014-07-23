#ifndef FILTER_H
#define FILTER_H

#include <gsl/gsl_matrix.h>
#include <QApplication>
#include <QImage>
#include <QColor>
#include <vector>
#include <fstream>

class Filter
{
public:
    Filter();
    void run();
    void loadFilterImage(QString path);
    void initializeIntensities();
    void roundIntensityValues();
    void drawIntensitySpace();

    gsl_matrix *input;
    gsl_matrix *output;

    int getMaxIntLong() { return maxIntLong; }
    void setMaxIntLong(int val) { this->maxIntLong=val; }
    void setMaxIntShort(int val) { this->maxIntShort=val; }

private:

    int maxIntLong;
    int maxIntShort;


    double precision = 0.1;
    int nr_of_filters = 1;

    QString pathIntensitySpace;

    std::vector<QImage> QImgFilters;

};

#endif // FILTER_H
