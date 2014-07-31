#ifndef PAIRFINDER_H
#define PAIRFINDER_H

#include <algorithm>
#include <vector>
#include <QString>
#include <QTextStream>
#include <fcntl.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <QtXml/QtXml>
#include <QtXml/QDomDocument>

#include <QThread>
#include <QEventLoop>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include "sdmixer.h"



class PairFinder : public QObject
{
    Q_OBJECT
public:
    struct Localization {

        double xLong, yLong, zLong;
        double xShort, yShort, zShort;
        double LongIntensity, ShortIntensity;
        int frame;

        int filter = 0;

    };

    PairFinder(sdmixer *s, QString f);
    void getHeader();
void Tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string delimiters);
    QString removeCharacters(QString input, char chars[])
    {
        std::string str = input.toStdString();
        for (unsigned int i = 0; i < strlen(chars); ++i)
        {
           str.erase (std::remove(str.begin(), str.end(), chars[i]), str.end());
        }
        return QString::fromStdString(str);
    }

    void FindPairs(int last_frame=-1);
    int loadFile();

    int getDimensions() {return dimensions;}

    int getNrOfDifferentFrames() { return NrOfDifferentFrames; }
    int getMultipleCounter() { return multiple_counter; }



    void loadInputFile();

signals:
    void finished();

public slots:
    void doWork();

private:
    sdmixer *sdm;
    QString file;
    std::vector<double> input;

    bool _working=false;
    bool _abort=false;
    QMutex mutex;

    int dimensions=0;
    double Offset[3]={0};
    double Epsilon[3]={0};


    std::vector<Localization> output_file;

    double min_x=0, max_x=0;
    double min_y=0, max_y=0;
    double min_z=0, max_z=0;

    int numpairs = 0;

    int xCol = 0;
    int yCol = 1;
    int zCol = 2;

    bool LeftRight = true;
    int ShortChannel = 1;

    std::vector<int> grouped_rows;

    int FrameColumn=4;
    int IntensityColumn=5;

    int NrOfDifferentFrames = 0;
    int multiple_counter = 0;
    bool create_output = false;
    int row_stop = 0;

    int rawDataCols = 0;
    int rawDataRows = 0;

    QString header;

};

#endif // PAIRFINDER_H
