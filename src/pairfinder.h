#ifndef PAIRFINDER_H
#define PAIRFINDER_H

#include <algorithm>
#include <vector>
#include <QString>
#include <QTextStream>


#include <QtXml/QtXml>
#include <QtXml/QDomDocument>

#include <QThread>
#include <QEventLoop>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

#include "sdmixer.h"
//class sdmixer


class PairFinder : public QObject
{
    Q_OBJECT
public:

    PairFinder(sdmixer *s, QString f);
    void getHeader(QString header_file="");
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

    void FindPairs(bool fishing, int last_frame=-1);
    int loadFile();

    int getDimensions() {return dimensions;}

    int getNrOfDifferentFrames() { return NrOfDifferentFrames; }
    int getMultipleCounter() { return multiple_counter; }
    sdmixer::min_max getMinMaxValues() { return min_maxValues; }

    void loadInputFile();
    void saveFile();
    void writeHeader(QTextStream &out, bool AppendHeader=false);

    void startGrouping();


signals:
    void finished();

public slots:
    void doWork();

private:
    struct fishing_run{
        double xOffset;
        double yOffset;
        double zOffset;
        int numpairs;
        double getOffset(int dim){
            if(dim==0)
                return xOffset;
            if(dim==1)
                return yOffset;
            if(dim==2)
                return zOffset;
            return 0;
        }
        void setOffset(int dim, double val)
        {
            if(dim==0)
                xOffset=val;
            if(dim==1)
                yOffset=val;
            if(dim==2)
                zOffset=val;
        }

    };

    //std::vector<sdmixer::Localization> output_file;
    bool canceled=false;
    // bool _working=false;
    // bool _abort=false;
     //QMutex mutex;
     //bool create_output = false;

    sdmixer::Columns columns;
    sdmixer::input_file_t INPUT_FILE;

    //internal Variables
    std::vector<fishing_run> fishing_results;

    std::vector<sdmixer::Localization> grouping_input;

    sdmixer *sdm;
    QString file, fileName, outputFile, groupoutFile;
    QString PairFinderSuffix = "_pairs_out.txt";
    QString GroupedFileSuffix = "_grouped_out.txt";
    QString output_dir;

    std::vector<double> input;
    int row_stop = 0;


    // Settings:
    double NM_PER_PX = 0;
    int dimensions=0;
    double Offset[3]={0};
    double Epsilon[3]={0}; 

    int xCol = 0;
    int yCol = 1;
    int zCol = 2;
    int FrameColumn=3;
    int IntensityColumn=4;

    sdmixer::min_max min_maxValues;
    sdmixer::fishing fishing_settings;
    sdmixer::offset_units offset_units;


    bool LeftRight = false;
    int ShortChannel = 1;


    //Output
    QString header;
    int numpairs = 0;

    std::vector<int> grouped_rows;

    int NrOfDifferentFrames = 0;
    int multiple_counter = 0;
    int rawDataCols = 0;
    int rawDataRows = 0;

    bool runGrouping;
    double groupingRadius;
    QString groupingUnits;


};

#endif // PAIRFINDER_H
