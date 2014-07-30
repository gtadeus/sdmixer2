#include "pairfinder.h"

void PairFinder::doWork() {
    // allocate resources using new here
    qDebug()<<"started";
    //loadFile(_file);
    load2();
    emit finished();
}
void PairFinder::load2(){
    std::ifstream infile;
    infile.open(_file.toStdString(), std::ifstream::in);
    std::string line;
    int counter=0;

    std::string firstLine;
    getline (infile, firstLine);
   // is.read()

/*
    while (std::getline(infile, line))
    {
        double d;
        std::istringstream lineStream(line);

        while (lineStream >> d)
        {
            input.push_back(d);
        }

        ++counter;
    }*/

    qDebug() << counter;
}

PairFinder::PairFinder(sdmixer *s, QString file)
{
    getHeader(file);
    //qDebug() << "dimensions: " << dimensions;

    if(s->getForce2D())
        dimensions=2;

    for(int i=0; i< dimensions; ++i)
    {
        Offset[i]=s->getOffset(i);
        Epsilon[i]=s->getEpsilon(i);
    }
    s->writeToConsole("loading file...");
    _file=file;
    //loadFile(file);
    //s->writeToConsole("searching Pairs");

    /*FindPairs();

    QString str; str = "pairs found ";
    str.append(QString::number(numpairs));

    s->writeToConsole(str);*/

}
void PairFinder::getHeader(QString file)
{
    std::ifstream ifs(file.toLatin1());
    std::string firstLine;
    getline (ifs, firstLine);

    ifs.close();

    header = removeCharacters(QString::fromStdString(firstLine), "#");

    QDomDocument qd;

    qd.setContent(header);
    QDomElement element = qd.documentElement();
    for(QDomNode n = element.firstChild(); !n.isNull(); n = n.nextSibling())
    {
        QDomElement e = n.toElement();
        if( e.tagName() == "field" )
        {
            if( e.attribute("identifier").contains("Position"))
            {
                ++dimensions;

                if(e.attribute("identifier").contains("1"))
                {
                    min_y = removeCharacters(e.attribute("min")," m").toDouble();
                    max_y = removeCharacters(e.attribute("max")," m").toDouble();
                }
                else if(e.attribute("identifier").contains("2"))
                {
                    min_z = removeCharacters(e.attribute("min")," m").toDouble();
                    max_z = removeCharacters(e.attribute("max")," m").toDouble();
                }
                else
                {
                    min_z = removeCharacters(e.attribute("min")," m").toDouble();
                    max_z = removeCharacters(e.attribute("max")," m").toDouble();
                }
            }
        }
    }
    /*qDebug() << min_x << "  " << max_x;
    qDebug() << min_y << "  " << max_y;
    qDebug() << min_z << "  " << max_z;*/
}

void PairFinder::FindPairs(int last_frame)
{
        int curr_row=0;
        int increment = 1;
        int endrow = 0;

        if (last_frame != -1)
            endrow = rawDataRows;
        else
            endrow = last_frame;

        while (curr_row < endrow)
        {
            double EllipsoidSumR=0;
            double EllipsoidSumL=0;
            Localization loc;

            if( input[curr_row*rawDataCols+FrameColumn] == input[(curr_row+increment)*rawDataCols+FrameColumn] )
            {
                for (int d = 0; d < dimensions; ++d)
                {
                    double tempL = ((input[curr_row*rawDataCols + d] - Offset[d]) - input[(curr_row+increment)*rawDataCols+d]);
                    double tempR = (input[(curr_row+increment)*rawDataCols+d] - Offset[d]) - input[curr_row*rawDataCols + d];
                    tempL*=tempL;
                    tempR*=tempR;
                    tempL/=(Epsilon[d]*Epsilon[d]);
                    tempR/=(Epsilon[d]*Epsilon[d]);
                    EllipsoidSumL += tempL;
                    EllipsoidSumR += tempR;
                }
                if (EllipsoidSumL <= 1 || EllipsoidSumR <= 1)
                {
                    int compare_col;
                    if(LeftRight)
                        compare_col = xCol;
                    else
                        compare_col = yCol;

                    bool factorShort;
                    if(ShortChannel == 1)
                    {
                        if(input[curr_row*rawDataCols+compare_col] < input[(curr_row+increment)*rawDataCols+compare_col])
                            factorShort = false;
                        else
                            factorShort = true;
                    }
                    else
                    {
                        if(input[curr_row*rawDataCols+compare_col] < input[(curr_row+increment)*rawDataCols+compare_col])
                            factorShort = true;
                        else
                            factorShort = false;
                    }
                    loc.xShort=input[(curr_row+factorShort*increment)*rawDataCols+xCol];
                    loc.yShort=input[(curr_row+factorShort*increment)*rawDataCols+yCol];
                    if(dimensions>2)
                        loc.zShort=input[(curr_row+factorShort*increment)*rawDataCols+zCol];
                    loc.ShortIntensity=input[(curr_row+factorShort*increment)*rawDataCols+FrameColumn];
                    loc.xLong=input[(curr_row+(!factorShort)*increment)*rawDataCols+xCol];
                    loc.yLong=input[(curr_row+(!factorShort)*increment)*rawDataCols+yCol];
                    if(dimensions>2)
                        loc.zLong=input[(curr_row+(!factorShort)*increment)*rawDataCols+zCol];
                    loc.LongIntensity=input[(curr_row+(!factorShort)*increment)*rawDataCols+FrameColumn];

                    output_file.push_back(loc);
                    grouped_rows.push_back(curr_row);
                    grouped_rows.push_back(curr_row+increment);

                }
                if (increment != rawDataRows)
                    increment++;
            }
            else
            {
                curr_row++;
                increment = 1;
            }
        }

    std::sort ( grouped_rows.begin(), grouped_rows.end() );
    multiple_counter = std::abs(std::distance(std::unique ( grouped_rows.begin(), grouped_rows.end()), grouped_rows.end()));
}

int PairFinder::loadFile(QString path)
{
    // This is the fastest way to count lines in file
    // from the UNIX tool wget
    qDebug()<<path;
    const char *fname = path.toStdString().c_str();
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
           ++rawDataRows;
    }
    // Lines counted
    std::ifstream ifs(fname);
    std::string firstLine, secondLine;
    getline (ifs, firstLine);
    getline (ifs, secondLine);

    ifs.close();
    //determine column number from second line
    std::stringstream countColsStream(secondLine);

    while (countColsStream >> dd)
    {
        ++rawDataCols;
    }
    // open file again...

    FILE *infile = fopen(fname, "r");
    char buffer[4096];
    input.reserve(rawDataRows*rawDataCols);

    //skip first Line
    fseek ( infile , int(firstLine.length()), SEEK_SET );

    while (fgets(buffer, sizeof(buffer), infile))
    {
            double d;
            std::stringstream lineStream(buffer);

            while (lineStream >> d)
                input.push_back(d);
    }
    qDebug() << input[rawDataRows*rawDataCols];
    return 1;

}

