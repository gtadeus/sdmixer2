#include "pairfinder.h"
#include <algorithm>

PairFinder::PairFinder()
{
    //set dim
    //set Offset
    //set stoprow
    //set rows & cols
}

void PairFinder::FindPairs()
{
        int curr_row=0;

        int last_frame = -1;
        int increment = 1;

        int endrow = rawDataRows;

        std::vector< std::vector<double> > vecPairsOut;
        std::vector<int> grouped_rows;

        while (curr_row < endrow)
        {
            double EllipsoidSumR=0;
            double EllipsoidSumL=0;
            if (last_frame != gsl_matrix_get(input, curr_row, FrameColumnC))
            {
                last_frame = gsl_matrix_get(input, curr_row, FrameColumnC);
                NrOfDifferentFrames++;
            }

            if( gsl_matrix_get(input, curr_row, FrameColumnC) == gsl_matrix_get(input, curr_row + increment, FrameColumnC) )
            {
                for (int d = 0; d < dimensions; ++d)
                {
                    double tempL = ((gsl_matrix_get(input, curr_row, d) - Offset[d]) - gsl_matrix_get(input, curr_row+increment, d));
                    double tempR = ((gsl_matrix_get(input, curr_row+increment, d) - Offset[d]) - gsl_matrix_get(input, curr_row, d));
                    tempL*=tempL;
                    tempR*=tempR;
                    tempL/=(Epsilon[d]*Epsilon[d]);
                    tempR/=(Epsilon[d]*Epsilon[d]);
                    EllipsoidSumL += tempL;
                    EllipsoidSumR += tempR;
                }
                if (EllipsoidSumL <= 1)
                {
                    std::vector<double> temp;

                    for(int i = 0; i < rawDataCols; ++i)
                        temp.push_back(gsl_matrix_get(input, curr_row, i));

                    for(int i = 0; i < rawDataCols; ++i)
                        temp.push_back(gsl_matrix_get(input, curr_row+increment, i));

                    vecPairsOut.push_back(temp);

                    numpairs++;

                    grouped_rows.push_back(curr_row);
                    grouped_rows.push_back(curr_row+increment);
                }
                else if (EllipsoidSumR <= 1)
                {
                    std::vector<double> temp;

                    for(int i = 0; i < rawDataCols; ++i)
                        temp.push_back(gsl_matrix_get(input, curr_row+increment, i));

                    for(int i = 0; i < rawDataCols; ++i)
                        temp.push_back(gsl_matrix_get(input, curr_row, i));

                    vecPairsOut.push_back(temp);
                    numpairs++;

                    grouped_rows.push_back(curr_row);
                    grouped_rows.push_back(curr_row+increment);
                }

                if (increment != rawDataRows)
                    increment++;

                 //EllipsoidSumR=0;
            }
            else
            {
                curr_row++;
                increment = 1;
                //EllipsoidSumR=0;
            }
        }

    std::sort ( grouped_rows.begin(), grouped_rows.end() );
    multiple_counter = std::abs(std::distance(std::unique ( grouped_rows.begin(), grouped_rows.end()), grouped_rows.end()));

    if ( create_output )
    {
    output = gsl_matrix_alloc(numpairs, rawDataCols*2);

    for(int i = 0; i < vecPairsOut.size(); ++i)
         for(int j=0; j< vecPairsOut[i].size(); ++j)
                gsl_matrix_set(output, i, j, vecPairsOut[i][j]);

    }

}



