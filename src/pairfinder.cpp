#include "pairfinder.h"
//#include "reconstructor.h"
#include "sdmixer.h"

void PairFinder::startGrouping()
{
    qDebug() << "starting grouping";
    qDebug() << "radius: " << groupingRadius;

    sdm->clearLocalizations();
    std::vector< std::vector<int> > local_nodes;
    std::vector<sdmixer::Localization> centered;
    //std::vector< std::vector<double> > centers;
    int j = 0;
    std::vector<int> active_old;
    std::vector<int> active_new;
    bool found = false;
    std::vector<sdmixer::Localization>::iterator it;

    int i=0;
    for ( it = grouping_input.begin(); it < grouping_input.end(); ++it)
    {
        sdmixer::Localization loc = *it;
        std::vector<sdmixer::Localization>::iterator prev_it;
        prev_it = it;
        --prev_it;
        sdmixer::Localization prev_loc = *prev_it;
        if ( it == grouping_input.begin() || loc.frame != prev_loc.frame )
        {
            active_old = active_new;
            active_new.clear();
        }
        found = false;
        for_each(active_old.begin(), active_old.end(), [&](int k)
        {
            if(found == false)
            {
                double nodes_weightShort;
                double nodes_weightLong;
                double *dxShort = new double[dimensions];
                double *dxLong = new double[dimensions];
                //double dx[dimensions];
                double EllipsoidSumShort=0;
                double EllipsoidSumLong=0;

                for (int d = 0; d < dimensions; ++d)
                {
                    //dx[d] = centers[k][d]-pData[i+d*rows];
                    dxShort[d] = centered[k].getShortDim(d)-loc.getShortDim(d);
                    dxLong[d] = centered[k].getLongDim(d)-loc.getLongDim(d);
                }

                for (int d = 0; d < dimensions; ++d)
                {
                    EllipsoidSumShort += (dxShort[d]/groupingRadius)*(dxShort[d]/groupingRadius);
                    EllipsoidSumLong += (dxLong[d]/groupingRadius)*(dxLong[d]/groupingRadius);
                }
                //qDebug() << "dxShort[0]: " << dxShort[0];
                //qDebug() << "EllipsoidSumLong: " << EllipsoidSumLong;
                if (EllipsoidSumShort < 1 && grouping_input[local_nodes[k][local_nodes[k].size()-1]].frame
                        - loc.frame<= 1)
                {
                    //qDebug() << "EllipsoidSumShort: " << EllipsoidSumShort;
                    //qDebug() << "EllipsoidSumLong: " << EllipsoidSumLong;
                    double *tempShort = new double[dimensions];
                    double *tempLong = new double[dimensions];
                    nodes_weightShort = 0;
                    nodes_weightLong = 0;

                    for_each(local_nodes[k].begin(), local_nodes[k].end(), [&](int nodeIdx)
                    {
                        nodes_weightShort += grouping_input[nodeIdx].ShortIntensity;
                        nodes_weightLong += grouping_input[nodeIdx].LongIntensity;
                    });
                    for (int d = 0; d < dimensions; ++d)
                    {
                        tempShort[d]=nodes_weightShort*centered[k].getShortDim(d);
                        tempShort[d]+= (loc.ShortIntensity * loc.getShortDim(d));
                        tempShort[d]/=(nodes_weightShort + loc.ShortIntensity);
                        centered[k].setShortDim(d, tempShort[d]);

                        tempLong[d]=nodes_weightLong*centered[k].getLongDim(d);
                        tempLong[d]+= (loc.LongIntensity * loc.getLongDim(d));
                        tempLong[d]/=(nodes_weightShort + loc.LongIntensity);

                        //centered[k].setLongDim(d, tempLong[d]);
                    }
                    local_nodes[k].push_back(i);
                    active_new.push_back(k);
                    found = true;
                }
            }
        });
        if(found == false)
        {
            std::vector<int> tmp;
            tmp.push_back(i);
            local_nodes.push_back(tmp);
            active_new.push_back(j);

            centered.push_back(loc);

            j++;
        }
        ++i;
    }
    for ( std::vector<  std::vector<int> >::size_type u = 0; u < local_nodes.size(); u++)
    {
        double Navg=0;
        double Nsum_second = 0;
        int tmp_groupSize=0;
        for ( std::vector<int>::size_type v = 0; v < local_nodes[u].size(); v++)
        {
            Navg+=grouping_input[local_nodes[u][v]].ShortIntensity;
                    //PairsOut(local_nodes[u][v], yColIntC);
            Nsum_second += grouping_input[local_nodes[u][v]].LongIntensity;
                    //PairsOut(local_nodes[u][v], xColIntC);
            tmp_groupSize+=1;
        }
        centered[u].ShortIntensity = Navg;
        centered[u].LongIntensity = Nsum_second;
        sdm->pushBackLocalization(centered[u]);
        /*vecIntLeft.push_back(Navg);
        vecIntRight.push_back(Nsum_second);
        groupSize.push_back(tmp_groupSize);
        if(tmp_groupSize>1)
            realGroupSize.push_back(tmp_groupSize);*/


    }
    QFile group_out(groupoutFile);
    QTextStream out(&group_out);
    group_out.open(QIODevice::WriteOnly | QIODevice::Text);
    writeHeader(out);

    //std::vector<sdmixer::Localization>::iterator it;
    for ( it = centered.begin(); it < centered.end(); ++it)
    {
        sdmixer::Localization loc = *it;
        for (int i = 0; i < dimensions; ++i)
            out << loc.getShortDim(i) << " ";
        out << loc.ShortIntensity << " ";
        out << loc.frame << " ";
        for (int i = 0; i < dimensions; ++i)
            out << loc.getLongDim(i) << " ";
        out << loc.LongIntensity << "\n";

    }
    group_out.close();

    qDebug() << "grouped " << grouping_input.size()-centered.size() << " coordinates";
    qDebug() << "now there are " << centered.size() << " pairs";

    std::stringstream ss;
    ss << "grouped into " << centered.size() << " pairs";
    sdm->writeToLogFile(QString::fromStdString(ss.str()));

    qDebug() << "finished grouping";



}

void PairFinder::doWork() {

    qDebug()<<"started file loading in new thread";

    loadInputFile();
    std::stringstream ss1, ss2;
    ss1 << "Offset (nm): ";
    ss2 << "Offset (px): ";

    for (int i = 0; i < dimensions; ++i)
    {
        ss1 << Offset[i] << " ";
        ss2 << Offset[i]/NM_PER_PX << " ";
    }
    sdm->writeToLogFile(QString::fromStdString(ss1.str()));
    sdm->writeToLogFile(QString::fromStdString(ss2.str()));

    if(canceled)
    {
        sdm->setStartDemixingButtonEnabled(true);
        emit finished();
        return;
    }

    if(fishing_settings.run)
    {
        sdm->writeToLogFile("searching optimal offset");
        qDebug() << "starting fishing";
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

        std::stringstream ss1, ss2;
        ss1 << "new Offset (nm): ";
        ss2 << "new Offset (px): ";

        for (int i = 0; i < dimensions; ++i)
        {
            //new Offset:
            Offset[i] = max_pairs.getOffset(i);    
            qDebug() << "new Offset : " << Offset[i];
            ss1 << Offset[i] << " ";
            ss2 << Offset[i]/NM_PER_PX << " ";
        }
        sdm->writeToLogFile(QString::fromStdString(ss1.str()));
        sdm->writeToLogFile(QString::fromStdString(ss2.str()));
    }

    FindPairs(false);

    qDebug() << "found " << numpairs << " pairs";
    if ( numpairs == 0)
    {
        //QMessageBox msgBox;
       // msgBox.setText("");

       // msgBox.critical(0,"SDmixer Error","No pairs found! Please check input file and offset settings.");
        /*int ret = msgBox.exec();
        if(ret == QMessageBox::Ok)*/
        emit error("No pairs found! Please check input file and offset settings.");
        return;
    }


    std::ostringstream os; os << "found " << numpairs << " pairs." ;
    sdm->writeToLogFile(QString::fromStdString(os.str()));

    if(runGrouping)
    {
        qDebug() << "starting grouping";
        sdm->writeToLogFile("starting Grouping");
        startGrouping();
    }

    sdm->setPF_min_maxValues(min_maxValues);

    //saveFile();

    emit finished();

}
void PairFinder::writeHeader(QTextStream &out, bool AppendFilter)
{
    QDomDocument pf_header;

    QDomElement paired_localizations = pf_header.createElement("PairedLocalizations");
    paired_localizations.setAttribute("name", "sdmixer");
    paired_localizations.setAttribute("version", QString::number(SDMIXER_VERSION) );
    pf_header.appendChild(paired_localizations);


    QString current_dim;

    for(int i = 0; i < dimensions; ++i)
    {
        QDomElement shortLoc = pf_header.createElement("field");
        current_dim = QString::number(i);
        QString shortPosition = "Short Position-";
        shortPosition.append(current_dim);
        shortLoc.setAttribute("identifier", shortPosition);
        QString min = QString::number(min_maxValues.getMin(i));
        min.append("m");
        QString max = QString::number(min_maxValues.getMax(i));
        max.append("m");
        shortLoc.setAttribute("min", min);
        shortLoc.setAttribute("max", max);
        paired_localizations.appendChild(shortLoc);
    }

    QDomElement ShortIntensity = pf_header.createElement("field");
    ShortIntensity.setAttribute("identifier", "Short Amplitude");
    paired_localizations.appendChild(ShortIntensity);

    QDomElement frame = pf_header.createElement("field");
    frame.setAttribute("identifier", "ImageNumber");
    paired_localizations.appendChild(frame);

    for(int i = 0; i < dimensions; ++i)
    {
        QDomElement longLoc = pf_header.createElement("field");
        current_dim = QString::number(i);
        QString longPosition = "Long Position-";
        longPosition.append(current_dim);
        longLoc.setAttribute("identifier", longPosition);
        QString min = QString::number(min_maxValues.getMin(i));
        min.append("m");
        QString max = QString::number(min_maxValues.getMax(i));
        max.append("m");
        longLoc.setAttribute("min", min);
        longLoc.setAttribute("max", max);
        paired_localizations.appendChild(longLoc);
    }

    QDomElement LongIntensity = pf_header.createElement("field");
    LongIntensity.setAttribute("identifier", "Long Amplitude");
    paired_localizations.appendChild(LongIntensity);

    if(AppendFilter)
    {
        QDomElement Filter = pf_header.createElement("field");
        Filter.setAttribute("identifier", "Filter");
        paired_localizations.appendChild(Filter);
    }


    QString header = pf_header.toString();
    header.replace("\n", "");
    header.insert(0, "#");

    out << header << "\n";
}

void PairFinder::saveFile()
{
    QFile pairs_out(outputFile);
    pairs_out.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&pairs_out);
    sdmixer::Localization loc;

    for( int i = 0; i < numpairs; ++i)
    {
        for(int j = 0; j < dimensions; ++j)
        {
            out << loc.getShortDim(i) << " ";
        }
        out << loc.ShortIntensity << " ";
        out << loc.frame << " ";
        for(int j = 0; j < dimensions; ++j)
        {
            out << loc.getLongDim(i) << " ";
        }
        out << loc.LongIntensity << "\n";
    }

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

    int row = 0;

    while (!in.atEnd())
    {
        line = in.readLine();

        std::vector<std::string> v;
        std::string str = line.toStdString();
        str.erase(str.begin(),
                      std::find_if(str.begin(), str.end(),
                                   std::bind1st(std::not_equal_to<char>(), ' ')));
        boost::split(v, str, boost::is_any_of("\t "));

        std::vector<double> lineVec;

        for (auto i : v)
        {
            lineVec.push_back(strtod(i.c_str(), NULL));
        }

        input.push_back(lineVec);
        row++;
    }
    std::sort(input.begin(), input.end(),
              [&](const std::vector<double>& a, const std::vector<double>& b) {
      return a[FrameColumn] < b[FrameColumn];
    });

    qDebug() << "total: " << input.size() ;
    rawDataRows=input.size();
    //rawDataRows=row;
    std::stringstream ss;
    ss << "found " << rawDataRows << " localizations.";
    sdm->writeToLogFile(QString::fromStdString(ss.str()));

    //rawDataRows = input.size();
    qDebug() <<  "rows : " << rawDataRows ;
    //qDebug() <<  "firstelement: " << input[0] <<"  last element: " << input[input.size()-1] ;

}

PairFinder::PairFinder(sdmixer *s, QString f)
{
    sdm = s;
    file = f;

    QFile qfile(f);
    qfile.open(QIODevice::ReadOnly| QIODevice::Text);

    QTextStream in(&qfile);
    header = in.readLine();
    qfile.close();

    sdm->getHeader(header, columns, min_maxValues, INPUT_FILE);
    dimensions = columns.dimensions;
    rawDataCols = columns.rawDataCols;

    if(sdm->getForce2D())
        dimensions=2;

    xCol = columns.x;
    yCol = columns.y;
    if(dimensions>2)
        zCol = columns.z;
    FrameColumn = columns.frame;
    IntensityColumn = columns.Amplitude;

    //getHeader();
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

    fileName = fi.completeBaseName();
    outputFile = output_dir;
    outputFile.append("/");
    outputFile.append(fileName);
    outputFile.append(PairFinderSuffix);
    qDebug() << output_dir;

    runGrouping = sdm->getRunGrouping();
    groupingRadius = sdm->getGroupingRadius();
    groupingUnits = sdm->getGroupingUnits();
    if( !groupingUnits.compare("px") )
    {
        groupingRadius*=NM_PER_PX;
    }
    groupoutFile = output_dir;
    groupoutFile.append("/");
    groupoutFile.append(fileName);
    groupoutFile.append(GroupedFileSuffix);
    qDebug() << "group_out file: " << groupoutFile;


    UnpairedOut = sdm->getUnpairedOut();
    unpairedFile = output_dir;
    unpairedFile.append("/");
    unpairedFile.append(fileName);
    unpairedFile.append(UnpairedFileSuffix);
}

void PairFinder::getHeader(QString header_file)
{
    if(header_file.isEmpty())
        header_file = file;
    std::ifstream ifs(header_file.toLatin1());
    std::string firstLine;
    getline (ifs, firstLine);

    ifs.close();

    header = removeCharacters(QString::fromStdString(firstLine), "#");

    QDomDocument qd;

    qd.setContent(header);

    QDomElement element = qd.documentElement();
    dimensions=0;

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
    if(dimensions > 3)
        dimensions/=2;
    qDebug() << "PairFinder::getHeader() dimensions: " << dimensions;
    qDebug() << "PairFinder::getHeader() rawDataCols: " << rawDataCols;
    qDebug() << "PairFinder::getHeader() max Values from header";
    qDebug() << min_maxValues.min_x << "  " << min_maxValues.max_x;
    qDebug() << min_maxValues.min_y << "  " << min_maxValues.max_y;
    qDebug() << min_maxValues.min_z << "  " << min_maxValues.max_z;
}

void PairFinder::FindPairs(bool fishing, int last_frame)
{

    std::vector<int> grouped_rows;

    QFile pairs_out(outputFile);
    QTextStream out(&pairs_out);

    if(!fishing)
    {
        qDebug()<<"searching for pairs...";
        pairs_out.open(QIODevice::WriteOnly | QIODevice::Text);
        writeHeader(out);
    }


    int number_unpaired = 0;
    std::vector<UnpairedLocalization> unpaired_vec;

    numpairs=0;
    int curr_row=0;
    int increment = 1;
    int endrow = 0;

    if (last_frame == -1 || last_frame >= rawDataRows)
        endrow = rawDataRows;
    else
        endrow = last_frame;

    std::vector<int> pair_row_id;

    while (curr_row < endrow-1)
    {
        /*qDebug() << "curr_row: " << curr_row;
        qDebug() << "increment: " << increment;*/
        sdmixer::Localization loc;

            //qDebug() << input[curr_row*rawDataCols+FrameColumn];
            //qDebug() << input[(curr_row+increment)*rawDataCols+FrameColumn];
        //if( input[curr_row*rawDataCols+FrameColumn] == input[(curr_row+increment)*rawDataCols+FrameColumn] )
        //if(curr_row+increment < rawDataRows)
        if( input[curr_row][FrameColumn] == input[curr_row+increment][FrameColumn] )
        {
            double EllipsoidSumR=0;
            double EllipsoidSumL=0;
                //qDebug() << curr_row;
                //qDebug() << curr_row+increment;
            for (int d = 0; d < dimensions; ++d)
            {
                //double tempL = ((input[curr_row*rawDataCols + d] - Offset[d]) - input[(curr_row+increment)*rawDataCols+d]);
                //double tempR = (input[(curr_row+increment)*rawDataCols+d] - Offset[d]) - input[curr_row*rawDataCols + d];
                double tempL = ((input[curr_row][d] - Offset[d]) - input[curr_row+increment][d]);
                double tempR = (input[curr_row+increment][d] - Offset[d]) - input[curr_row][d];

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
                pair_row_id.push_back(curr_row);
                pair_row_id.push_back(curr_row+increment);
                //qDebug() << "saved " << curr_row << " & " << curr_row+increment;

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
                        //if(input[curr_row*rawDataCols+compare_col] < input[(curr_row+increment)*rawDataCols+compare_col])
                    if(input[curr_row][compare_col] < input[curr_row+increment][compare_col])
                            factorShort = false;
                        else
                            factorShort = true;
                }
                else
                {
                    //if(input[curr_row*rawDataCols+compare_col] < input[(curr_row+increment)*rawDataCols+compare_col])
                    if(input[curr_row][compare_col] < input[curr_row+increment][compare_col])
                        factorShort = true;
                    else
                        factorShort = false;
                }
                /*loc.xShort=input[(curr_row+factorShort*increment)*rawDataCols+xCol];
                loc.yShort=input[(curr_row+factorShort*increment)*rawDataCols+yCol];
                if(dimensions>2)
                    loc.zShort=input[(curr_row+factorShort*increment)*rawDataCols+zCol];
                loc.ShortIntensity=input[(curr_row+factorShort*increment)*rawDataCols+IntensityColumn];
                loc.frame = input[curr_row*rawDataCols+FrameColumn];
                loc.xLong=input[(curr_row+(!factorShort)*increment)*rawDataCols+xCol];
                loc.yLong=input[(curr_row+(!factorShort)*increment)*rawDataCols+yCol];
                if(dimensions>2)
                    loc.zLong=input[(curr_row+(!factorShort)*increment)*rawDataCols+zCol];
                loc.LongIntensity=input[(curr_row+(!factorShort)*increment)*rawDataCols+IntensityColumn];
*/
                loc.xShort=input[curr_row+factorShort*increment][xCol];
                loc.yShort=input[curr_row+factorShort*increment][yCol];
                if(dimensions>2)
                    loc.zShort=input[curr_row+factorShort*increment][zCol];
                loc.ShortIntensity=input[curr_row+factorShort*increment][IntensityColumn];
                loc.frame = input[curr_row][FrameColumn];
                loc.xLong=input[curr_row+(!factorShort)*increment][xCol];
                loc.yLong=input[curr_row+(!factorShort)*increment][yCol];
                if(dimensions>2)
                    loc.zLong=input[curr_row+(!factorShort)*increment][zCol];
                loc.LongIntensity=input[curr_row+(!factorShort)*increment][IntensityColumn];

                if(!fishing)
                {
                    sdm->pushBackLocalization(loc);
                    if(runGrouping)
                        grouping_input.push_back(loc);

                    for(int ii = 0; ii < dimensions; ++ii)
                    {
                        out << loc.getShortDim(ii) << " ";
                    }
                    out << loc.ShortIntensity << " ";
                    out << loc.frame << " ";
                    for(int jj = 0; jj < dimensions; ++jj)
                    {
                        out << loc.getLongDim(jj) << " ";
                    }
                    out << loc.LongIntensity << "\n";
                }
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
            if ((curr_row+increment) < rawDataRows - 1)
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
    if(numpairs!=0 && !fishing)
    {
        std::sort ( grouped_rows.begin(), grouped_rows.end() );
        multiple_counter = std::abs(std::distance(std::unique ( grouped_rows.begin(), grouped_rows.end()), grouped_rows.end()));
        std::stringstream ss;
        ss << "found " << multiple_counter << " multiple pairs.";

        sdm->writeToLogFile(QString::fromStdString(ss.str()));
        qDebug() << "found " << multiple_counter << " multiple pairs." ;
    }
    if(!fishing)
        pairs_out.close();

    if(UnpairedOut && !fishing)
    {
        for (int nrow = 0; nrow < rawDataRows; ++nrow)
        {
            if(std::find(pair_row_id.begin(), pair_row_id.end(), nrow)==pair_row_id.end())
            {
                //qDebug() << nrow;

                UnpairedLocalization uloc;
                for(int ii = 0; ii < dimensions; ++ii)
                    uloc.set(ii, input[nrow][ii]);
                uloc.frame=input[nrow][FrameColumn];
                uloc.intensity=input[nrow][IntensityColumn];

                unpaired_vec.push_back(uloc);
                number_unpaired++;
            }
        }
        std::stringstream ss;
        ss << "found " << number_unpaired << " unpaired localizations.";

        sdm->writeToLogFile(QString::fromStdString(ss.str()));
        qDebug() << "found " << number_unpaired << " unpaired localizations.";

        QFile unpaired_out(unpairedFile);
        unpaired_out.open(QIODevice::WriteOnly | QIODevice::Text);
        QTextStream out(&unpaired_out);
        out << header << "\n";

        UnpairedLocalization loc;

        for( int i = 0; i < number_unpaired; ++i)
        {
            loc = unpaired_vec[i];
            for(int j = 0; j < dimensions; ++j)
            {
                out << loc.get(j) << " ";
            }
            out << loc.frame << " ";
            out << loc.intensity << "\n";
        }
        unpaired_out.close();
    }


    //sdmixer::log(sdm, os.str());

}


