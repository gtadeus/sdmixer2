#include "pairfinder.h"


PairFinder::PairFinder()
{
    //set dim
    //set Offset
    //set stoprow
    //set rows & cols
}

void PairFinder::FindPairs()
{/*
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

    }*/

}

int PairFinder::loadFile(QString path)
{
    const char *fname = path.toStdString().c_str();
    int lines = 0;
    int cols = 0;
    double dd;

    static const int BUFFER_SIZE = 16*1024;
    int fd = open(fname, O_RDONLY);
    if(fd == -1)
        return 0;
    char buf[BUFFER_SIZE + 1];
    while(size_t bytes_read = read(fd, buf, BUFFER_SIZE))
    {
        if(bytes_read == (size_t)-1)
            return 0;//handle_error("read failed");
        if (!bytes_read)
            break;

        for(char *p = buf; (p = (char*) memchr(p, '\n', (buf + bytes_read) - p)); ++p)
           ++lines;
    }

    std::ifstream ifs(fname);
    std::string firstLine, secondLine;
    getline (ifs, firstLine);
    getline (ifs, secondLine);
    ifs.close();

    std::stringstream countColsStream(secondLine);

    while (countColsStream >> dd)
    {
        ++cols;
    }


    FILE *infile = fopen(fname, "r");
    char buffer[4096];

    double inputdata[lines*cols];

    //gsl_matrix *K = gsl_matrix_alloc(lines,cols);
   // double* a = new double[lines*cols];

    int curr_line=0;

    fseek ( infile , int(firstLine.length()), SEEK_SET );
    int counter = 0;
    while (fgets(buffer, sizeof(buffer), infile))
    {
            double d;
            std::stringstream lineStream(buffer);

            int curr_col = 0;
            while (lineStream >> d)
            {
                inputdata[counter] = d;
                //ui->textEdit->append(QString::number(d) + "  " +  QString::number(curr_line) +  "  " + QString::number(curr_col) + "\n" );
                //gsl_matrix_set(K, curr_line, curr_col, d);
                //a[curr_line*cols+curr_col]=d;
                ++counter;
                ++curr_col;
            }

            ++curr_line;

    }


    QString qs(firstLine.c_str());

    /*ui->textConsole->append(QString::number(lines));
    ui->textConsole->append(QString::number(cols));
    ui->textConsole->append(qs);
    ui->textConsole->append(QString::number(inputdata[9 * cols + 7]));*/
    //QString message;
    //QTextStream out(&message);
    //out << timestamp() << " : File " << fname << " loaded, " << lines << " lines <br>";
    //console->append(message);
    return lines;
}

