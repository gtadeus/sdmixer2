#include "pairfinder.h"
#include "reconstructor.h"

void PairFinder::doWork() {
    // allocate resources using new here
    qDebug()<<"started file loading in new thread";
    loadInputFile();
    qDebug()<<"searching for pairs...";

    FindPairs();

    Reconstructor r(sdm, output_file);
    qDebug()<<"set min max";
    r.setMinMax(min_x, max_x, min_y, max_y, min_z, max_z);

    //r.XYZfromFilter();
    //std::vector<PairFinder::Localization>().swap(output_file);
    //r.getMinMax();
    r.setArray();
    r.setKernel();

    //r.Convolution();
    r.map8bit();
    r.outputTIFF();

    sdm->setStartDemixingButtonEnabled(true);
    emit finished();
}
void PairFinder::Tokenize(const std::string& str,
                      std::vector<std::string>& tokens,
                      const std::string delimiters)
{
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

void PairFinder::loadInputFile()
{
    sdm->writeToConsole("loading file...");
    QFile f(file);
    f.open(QIODevice::ReadOnly| QIODevice::Text);

    QTextStream in(&f);
    QString line = in.readLine();

    while (!in.atEnd())
    {
        line = in.readLine();

        std::vector<std::string> v;
        std::string str = line.toStdString();
        boost::split(v, str, boost::is_any_of("\t "));

        for (auto i : v)
        {
            input.push_back(strtod(i.c_str(), NULL));
        }
    }

    qDebug() << "total: " << input.size() ;
    rawDataRows=input.size()/rawDataCols;
    qDebug() <<  "rows : " << rawDataRows ;
    qDebug() <<  input[0] <<"  " << input[input.size()-1] ;

}

PairFinder::PairFinder(sdmixer *s, QString f)
{
    sdm = s;
    file = f;
    getHeader();
    qDebug() << "dimensions: " << dimensions;
    qDebug() << "columns: " << rawDataCols;

    if(s->getForce2D())
        dimensions=2;

    NM_PER_PX = sdm->getPixelSize();

    for(int i=0; i< dimensions; ++i)
    {
        Offset[i]=sdm->getOffset(i)*NM_PER_PX;
        Epsilon[i]=sdm->getEpsilon(i)*NM_PER_PX;

    }

    //loadFile(file);
    //s->writeToConsole("searching Pairs");

    /*FindPairs();

    QString str; str = "pairs found ";
    str.append(QString::number(numpairs));

    s->writeToConsole(str);*/
}
void PairFinder::getHeader()
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
            ++rawDataCols;
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
                    min_x = removeCharacters(e.attribute("min")," m").toDouble();
                    max_x = removeCharacters(e.attribute("max")," m").toDouble();
                }
            }
        }
    }
    qDebug() << "max Values from config";
    qDebug() << min_x << "  " << max_x;
    qDebug() << min_y << "  " << max_y;
    qDebug() << min_z << "  " << max_z;
}

void PairFinder::FindPairs(int last_frame)
{

        int curr_row=0;
        int increment = 1;
        int endrow = 0;

        if (last_frame == -1)
            endrow = rawDataRows;
        else
            endrow = last_frame;

        while (curr_row < endrow)
        {


            Localization loc;

            //qDebug() << input[curr_row*rawDataCols+FrameColumn];
            //qDebug() << input[(curr_row+increment)*rawDataCols+FrameColumn];

            if( input[curr_row*rawDataCols+FrameColumn] == input[(curr_row+increment)*rawDataCols+FrameColumn] )
            {
                double EllipsoidSumR=0;
                double EllipsoidSumL=0;
                //qDebug() << curr_row;
                //qDebug() << curr_row+increment;

                for (int d = 0; d < dimensions; ++d)
                {
                    double tempL = ((input[curr_row*rawDataCols + d] - Offset[d]) - input[(curr_row+increment)*rawDataCols+d]);
                    double tempR = (input[(curr_row+increment)*rawDataCols+d] - Offset[d]) - input[curr_row*rawDataCols + d];
                    //qDebug() <<"tempL: "<< tempL;
                    //qDebug() << "tempR: " << tempR;
                    tempL*=tempL;
                    tempR*=tempR;
                    tempL/=(Epsilon[d]*Epsilon[d]);
                    tempR/=(Epsilon[d]*Epsilon[d]);
                    EllipsoidSumL += tempL;
                    EllipsoidSumR += tempR;

                }
                if (EllipsoidSumL <= 1 || EllipsoidSumR <= 1)
                {
                    ++numpairs;
                    //qDebug() << EllipsoidSumL << " @ " << curr_row << " & " << (curr_row+increment);
                    //qDebug() << EllipsoidSumR << " @ " << curr_row << " & " << (curr_row+increment);
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
                if (increment <= rawDataRows)
                    increment++;
                else
                {
                    curr_row++;
                    increment = 1;
                }
            }
            else
            {
                curr_row++;
                increment = 1;
            }
        }

    std::sort ( grouped_rows.begin(), grouped_rows.end() );
    multiple_counter = std::abs(std::distance(std::unique ( grouped_rows.begin(), grouped_rows.end()), grouped_rows.end()));

    qDebug() << "found " << numpairs << " pairs." ;
}


