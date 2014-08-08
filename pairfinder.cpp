#include "pairfinder.h"
//#include "reconstructor.h"
#include "sdmixer.h"

void PairFinder::doWork() {

    qDebug()<<"started file loading in new thread";

    loadInputFile();
    if(canceled)
    {
        sdm->setStartDemixingButtonEnabled(true);
        emit finished();
        return;
    }

    if(fishing_settings.run)
    {
        qDebug() << "hier";
        /*if(fishing_settings.range > min_maxValues.max_x)
        {
            fishing_settings.range = min_maxValues.max_x/10;
            if(fishing_settings.range > min_maxValues.max_y)
                fishing_settings.range = min_maxValues.max_y/10;
        }*/
            double minOffsetX = Offset[0]-(fishing_settings.range/2);
            double minOffsetY = Offset[1]-(fishing_settings.range/2);
            double maxOffsetX = Offset[0]+(fishing_settings.range/2);
            double maxOffsetY = Offset[1]+(fishing_settings.range/2);
            int counterX = 0;
            int counterY = 0;
            while(Offset[0]<maxOffsetX)
            {
                Offset[0] = minOffsetX + fishing_settings.increment*counterX;
                while(Offset[1]<maxOffsetY)
                {

                    Offset[1] = minOffsetY + fishing_settings.increment*counterY;
                    //qDebug() << "current Offset : " << Offset[0] << " " << Offset[1];
                    FindPairs(true, fishing_settings.subset);
                    counterY++;
                }
                Offset[1]=minOffsetY;
                counterY=0;
                counterX++;
            }


        int index = 0;
        fishing_run max_pairs;
        std::vector<fishing_run>::iterator it;
        for( it = fishing_results.begin(); it != fishing_results.end(); ++it )
        {
            fishing_run current_result = *it;
            if (index == 0)
            {
                max_pairs = current_result;
            }
            if(max_pairs.numpairs < current_result.numpairs)
                max_pairs = current_result;

            ++index;
        }

        for (int i = 0; i < dimensions; ++i)
        {
            //new Offset:
            Offset[i] = max_pairs.getOffset(i);
            qDebug() << "new Offset : " << Offset[i];
        }
    }

    FindPairs(false);

    qDebug() << "found " << numpairs << " pairs";

    sdm->setPF_min_maxValues(min_maxValues);

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
    sdmixer::log(sdm, "loading file...");
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

    //qDebug() << "total: " << input.size() ;
    rawDataRows=input.size()/rawDataCols;
    //qDebug() <<  "rows : " << rawDataRows ;
    //qDebug() <<  "firstelement: " << input[0] <<"  last element: " << input[input.size()-1] ;

}

PairFinder::PairFinder(sdmixer *s, QString f)
{
    sdm = s;
    file = f;
    getHeader();
    //qDebug() << "dimensions: " << dimensions;
    //qDebug() << "columns: " << rawDataCols;

    if(s->getForce2D())
        dimensions=2;

    NM_PER_PX = sdm->getPixelSize();

    for(int i=0; i< dimensions; ++i)
    {
        offset_units = sdm->getOffsetUnits();

        if(!offset_units.getOffset(i).compare("px"))
            Offset[i]=sdm->getOffset(i)*NM_PER_PX;
        else
            Offset[i]=sdm->getOffset(i);

        if(!offset_units.getEpsilon(i).compare("px"))
            Epsilon[i]=sdm->getEpsilon(i)*NM_PER_PX;
        else
            Epsilon[i]=sdm->getEpsilon(i);

        //qDebug() << "offset : " << Offset[i];
        //qDebug() << "epsilon : " << Epsilon[i];

    }
    fishing_settings = sdm->getFishing();

    if( !fishing_settings.unit.compare("px") )
    {
        fishing_settings.increment*=NM_PER_PX;
        fishing_settings.range*=NM_PER_PX;

    }

    if(!sdm->getCameraOrientation().compare("Top-Bottom"))
        LeftRight = false;
    else
        LeftRight = true;

    if(!sdm->getShortChannelPosition().compare("Top") ||
            !sdm->getShortChannelPosition().compare("Left"))
        ShortChannel=1;
    else
        ShortChannel=2;

    QFile qf(file);
    QFileInfo fi(qf);

    if(!sdm->getOutputDirectory().isEmpty())
        output_dir = sdm->getOutputDirectory();
    else
        output_dir=fi.path();

    fileName = fi.baseName();
    outputFile = output_dir;
    outputFile.append("/");
    outputFile.append(fileName);
    outputFile.append(PairFinderSuffix);
    //qDebug() << output_dir;

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
                    min_maxValues.min_y = removeCharacters(e.attribute("min")," m").toDouble();
                    min_maxValues.max_y = removeCharacters(e.attribute("max")," m").toDouble();
                }
                else if(e.attribute("identifier").contains("2"))
                {
                    min_maxValues.min_z = removeCharacters(e.attribute("min")," m").toDouble();
                    min_maxValues.max_z = removeCharacters(e.attribute("max")," m").toDouble();
                }
                else
                {
                    min_maxValues.min_x = removeCharacters(e.attribute("min")," m").toDouble();
                    min_maxValues.max_x = removeCharacters(e.attribute("max")," m").toDouble();
                }
            }
        }
    }
    /*qDebug() << "max Values from config";
    qDebug() << min_maxValues.min_x << "  " << min_maxValues.max_x;
    qDebug() << min_maxValues.min_y << "  " << min_maxValues.max_y;
    qDebug() << min_maxValues.min_z << "  " << min_maxValues.max_z;*/
}

void PairFinder::FindPairs(bool fishing, int last_frame)
{
    //qDebug()<<"searching for pairs...";
    numpairs=0;
    int curr_row=0;
    int increment = 1;
    int endrow = 0;

    if (last_frame == -1 || last_frame >= rawDataRows)
        endrow = rawDataRows;
    else
        endrow = last_frame;

    while (curr_row < endrow)
    {
        sdmixer::Localization loc;

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
                loc.ShortIntensity=input[(curr_row+factorShort*increment)*rawDataCols+IntensityColumn];
                loc.xLong=input[(curr_row+(!factorShort)*increment)*rawDataCols+xCol];
                loc.yLong=input[(curr_row+(!factorShort)*increment)*rawDataCols+yCol];
                if(dimensions>2)
                    loc.zLong=input[(curr_row+(!factorShort)*increment)*rawDataCols+zCol];
                loc.LongIntensity=input[(curr_row+(!factorShort)*increment)*rawDataCols+IntensityColumn];

                if(!fishing)
                    sdm->pushBackLocalization(loc);
                else
                {
                    fishing_run current_run;
                    for(int ii = 0; ii < dimensions; ++ii)
                        current_run.setOffset(ii, Offset[ii]);


                    current_run.numpairs = numpairs;
                    fishing_results.push_back(current_run);
                }
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
    if(numpairs!=0)
    {
        std::sort ( grouped_rows.begin(), grouped_rows.end() );
        multiple_counter = std::abs(std::distance(std::unique ( grouped_rows.begin(), grouped_rows.end()), grouped_rows.end()));
    }

    //std::ostringstream os; os << "found " << numpairs << " pairs." ;
    //sdmixer::log(sdm, os.str());
    //qDebug() << "found " << numpairs << " pairs." ;
}


