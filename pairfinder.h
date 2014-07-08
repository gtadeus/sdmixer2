#ifndef PAIRFINDER_H
#define PAIRFINDER_H

#include <gsl/gsl_matrix.h>

class PairFinder
{
public:
    PairFinder();
    void FindPairs();

    gsl_matrix *input;
    gsl_matrix *output;
    int dimensions;
    double Offset[];
    double Epsilon[];
    int FrameColumnC;
    int NrOfDifferentFrames;
    int multiple_counter;
    bool create_output = false;
    int row_stop = -1;
};

#endif // PAIRFINDER_H
