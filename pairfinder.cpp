#include "pairfinder.h"

PairFinder::PairFinder()
{
}
/*
int PairFinder::FindPairs(Matrix& input, Matrix& output, int dimensions,
                          double Offset[], double Epsilon[], int FrameColumnC,
                          int& NrOfDifferentFrames, int& multiple_counter, bool create_output = false, int row_stop = -1)
{
        int pair=0;
        int increment=1;
        NrOfDifferentFrames=0;
        multiple_counter=0;
        int numpairs = 0;
        double EllipsoidSum = 0;

        int last_frame = -1;

        const int raw_data_cols = input.cols();
        const int raw_data_rows = input.rows();
        int endrow = raw_data_rows;

        if(row_stop != -1)
            endrow = row_stop;


        double (*pFullData) = new double[raw_data_cols*raw_data_rows];
        double *fv_pFullData = input.fortran_vec();

        memcpy(pFullData, fv_pFullData, raw_data_cols*raw_data_rows*sizeof(double));

        std::vector< std::vector<double> > vecPairsOut;

        std::vector<int> grouped_rows;

        while (pair < endrow)
        {
            double EllipsoidSumR=0;
            double EllipsoidSumL=0;
            if (last_frame != pFullData[pair+FrameColumnC*raw_data_rows])
            {
                last_frame = pFullData[pair+FrameColumnC*raw_data_rows];
                NrOfDifferentFrames++;
            }

            if( (pFullData[pair+FrameColumnC*raw_data_rows] == pFullData[pair+increment+FrameColumnC*raw_data_rows]) )
            {
                    for (int d = 0; d < dimensions; ++d)
                    {
                        double tempL = ((pFullData[pair+d*raw_data_rows] - Offset[d]) - pFullData[pair+increment+d*raw_data_rows]);
                        double tempR = ((pFullData[pair+increment+d*raw_data_rows] - Offset[d]) - pFullData[pair+d*raw_data_rows]);

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

                    for(int i = 0; i < raw_data_cols; ++i)
                        temp.push_back(pFullData[pair+i*raw_data_rows]);

                    for(int i = 0; i < raw_data_cols; ++i)
                        temp.push_back(pFullData[pair+increment+i*raw_data_rows]);

                    vecPairsOut.push_back(temp);

                    numpairs++;

                    grouped_rows.push_back(pair);
                    grouped_rows.push_back(pair+increment);
                }
                else if (EllipsoidSumR <= 1)
                {
                    std::vector<double> temp;

                    for(int i = 0; i < raw_data_cols; ++i)
                        temp.push_back(pFullData[pair+increment+i*raw_data_rows]);

                    for(int i = 0; i < raw_data_cols; ++i)
                        temp.push_back(pFullData[pair+i*raw_data_rows]);

                    vecPairsOut.push_back(temp);
                    numpairs++;

                    grouped_rows.push_back(pair);
                    grouped_rows.push_back(pair+increment);
                }

                if (increment != raw_data_rows)
                    increment++;

                 EllipsoidSumR=0;

            }
            else
            {
                pair++;
                increment = 1;
                EllipsoidSumR=0;
            }

        }

    std::sort ( grouped_rows.begin(), grouped_rows.end() );
    multiple_counter = std::abs(std::distance(std::unique ( grouped_rows.begin(), grouped_rows.end()), grouped_rows.end()));

    if ( create_output )
    {
    output.resize(numpairs, raw_data_cols*2);

    for(int i = 0; i < vecPairsOut.size(); ++i)
         for(int j=0; j< vecPairsOut[i].size(); ++j)
               output(i,j)= vecPairsOut[i][j] ;

    }


    delete[] pFullData;

    return numpairs;
}
*/

void PairFinder::FindPairs()
{

}
