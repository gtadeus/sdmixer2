#ifndef PAIRFINDER_H
#define PAIRFINDER_H

#include <gsl/gsl_matrix.h>
#include <vector>

struct SDFile {

    double xLong, yLong, zLong;
    double xShort, yShort, zShort;
    double LongIntensity, ShortIntensity;
    int frame;
    std::vector<double> fitvalues;
    int filter;

};

class PairFinder
{
public:
    PairFinder();
    void FindPairs();


    gsl_matrix *input;
    gsl_matrix *output;

    void setDimensions(int dim) { this->dimensions = dim ; }
    void setOffset(std::vector< double > Offset)
    {
        this->Offset=Offset;
    }
    void setEpsilon(std::vector< double > Epsilon)
    {
        this->Epsilon=Epsilon;
    }
    void setFrameCol(int a) { this->FrameColumnC = a ; }

    int getNrOfDifferentFrames() { return NrOfDifferentFrames; }
    int getMultipleCounter() { return multiple_counter; }

    void setStopRow(int row_stop) { this->row_stop = row_stop ; }

private:
    int dimensions;
    std::vector< double > Offset;
    std::vector< double > Epsilon;

    int numpairs = 0;
    int FrameColumnC;
    int NrOfDifferentFrames = 0;
    int multiple_counter = 0;
    bool create_output = false;
    int row_stop = 0;

    int rawDataCols;
    int rawDataRows;

};

#endif // PAIRFINDER_H
