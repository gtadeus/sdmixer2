#include "filter.h"
//#include "pairfinder.h"
#include <QDebug>
#include <QPainter>
#include <QBrush>

Filter::Filter(sdmixer *s, std::vector<sdmixer::Localization> *data)
{
    qDebug() << "inititalizing Filter";
    this->sdm = s;
    /*for(auto i : data)
    {
        input.push_back(i);
    }*/

    input=data;
    init();


}
void Filter::init()
{
    PairFinder p(sdm, sdm->getCurrentFile());
    dimensions = p.getDimensions();
    qDebug() << "Filter dimensions: " << dimensions;
    if (sdm->getForce2D())
        dimensions = 2;

    precision =sdm->getPrecision();
    maxIntLong = sdm->getMaxIntLong();
    maxIntShort = sdm->getMaxIntShort();

    QString file = sdm->getCurrentFile();
    QFile qf(file);
    QFileInfo fi(qf);

    if(!sdm->getOutputDirectory().isEmpty())
        output_dir = sdm->getOutputDirectory();
    else
        output_dir=fi.path();

    fileName = fi.completeBaseName();
    if(fileName.contains("_pairs_out"))
        fileName.replace("_pairs_out", "");
    if(fileName.contains("_filter_out"))
        fileName.replace("_filter_out", "");
    if(fileName.contains("_grouped_out"))
        fileName.replace("_grouped_out", "");

    outputFile = output_dir;
    outputFile.append("/");
    outputFile.append(fileName);
    outputFile.append(FilterSuffix);
    qDebug() << output_dir;

    intensitySpaceFile = output_dir;
    intensitySpaceFile.append("/");
    intensitySpaceFile.append(fileName);
    intensitySpaceFile.append(IntensitySpaceSuffix);
    qDebug() << "IntSpaceFile : " << intensitySpaceFile;




}
void Filter::loadFile(QString str)
{
    qDebug() << "loading pairs_out " << str;


    QFile f(str);
    f.open(QIODevice::ReadOnly| QIODevice::Text);

    QTextStream in(&f);
    QString line = in.readLine();
    sdm->getHeader(line, columns, min_maxValues, INPUT_FILE);
    dimensions = columns.dimensions;
    //if(INPUT_FILE != sdmixer::FILTER_FILE)
        //qDebug() << "ERROR, invalid Filter File"
    int index_g=0;
    int v_size=0;

    while (!in.atEnd())
    {
        line = in.readLine();

        std::vector<std::string> v;
        std::vector<double> dbl_vec;
        std::string strLine = line.toStdString();
        strLine.erase(strLine.begin(),
                      std::find_if(strLine.begin(), strLine.end(),
                                   std::bind1st(std::not_equal_to<char>(), ' ')));

        boost::split(v, strLine, boost::is_any_of("\t "));

        int index = 0;
        for (auto i : v)
        {
            index++;
            dbl_vec.push_back(strtod(i.c_str(), NULL));
            //qDebug() << "index-1: " << index-1 <<  "dbl_vec: " << dbl_vec[index-1];
        }
        index_g = index;
        v_size = v.size();

        sdmixer::Localization loc;
        for(int i = 0; i < dimensions; ++i)
        {
            int indL = columns.getLongCol(i);
            int indS = columns.getShortCol(i);
            loc.setLongDim(i, dbl_vec[indL]);
            loc.setShortDim(i, dbl_vec[indS]);
        }
        if(INPUT_FILE == sdmixer::FILTER_FILE)
            loc.filter = dbl_vec[columns.filter];

        loc.ShortIntensity = dbl_vec[columns.ShortAmp];
        loc.LongIntensity = dbl_vec[columns.LongAmp];
        loc.frame = dbl_vec[columns.frame];

        sdm->pushBackLocalization(loc);
    }
    //qDebug() << "index_g " << index_g;
    //qDebug() << "v.size() " << v_size;
    input = sdm->getPfOutput();

    qDebug() << "loaded ";


}

Filter::Filter(sdmixer *s, QString str)
{
    qDebug() << "inititalizing Filter: " << str;
    this->sdm = s;
    force2D = sdm->getForce2D();
    sdm->clearLocalizations();
    init();
    doWorkLaterParameter = str;
    doWorkLater=true;


}

void Filter::doWork()
{
    qDebug() << "starting Filter in new Thread";
    qDebug()<< "filter dimensions: " <<dimensions;

    if(doWorkLater)
        loadFile(doWorkLaterParameter);

    FilterInputFiles = sdm->getFilterFiles();

    // TODO: write header
    QFile filter_out(outputFile);
    QTextStream out(&filter_out);
    filter_out.open(QIODevice::WriteOnly | QIODevice::Text);
    PairFinder p(sdm, sdm->getCurrentFile());
    p.writeHeader(out, true);

    int IntSpace_width = 0;
    int IntSpace_height = 0;
    int roundedY = 0;

    bool x_Short;
    if(sdm->getFilterOrientation().contains("x=Short"))
    {
        IntSpace_width = round(maxIntShort*precision);
        IntSpace_height = round(maxIntLong*precision);
        roundedY = round(maxIntLong*precision);
        x_Short=true;

    }
    else
    {
        IntSpace_width = round(maxIntLong*precision);
        IntSpace_height = round(maxIntShort*precision);
        roundedY = round(maxIntShort*precision);
        x_Short=false;
    }

    int ImageMaxX;
    int ImageMaxY;


    qDebug() << "roundedY: " << roundedY;

    uint32_t IntensitySpaceNPixel = IntSpace_width*IntSpace_height;
    float *IntensitySpaceArray = new float[IntensitySpaceNPixel];

    // initialize IntensitySpace
    for ( int i = 0; i < IntensitySpaceNPixel; ++i)
        IntensitySpaceArray[i]=0;

    int current_filter=1;

    float max_value_intSpace = 0;
    int max_value_index=0;
    float min_value_intSpace = 99999;
    int min_value_index=0;
    QImage IntSpaceOut(IntSpace_width, IntSpace_height, QImage::Format_RGB32);
    IntSpaceOut.fill(Qt::white);

    std::vector<QString>::iterator ff;

    for ( ff = FilterInputFiles.begin(); ff < FilterInputFiles.end(); ++ff)
    {
        //sumChannels[current_filter-1]=0;

        QString i = *ff;
        QImage img(i);

        qDebug() << "filter file " << i;
        qDebug() << "current_filter: " << current_filter;

        //int sum = 0;

        std::vector<sdmixer::Localization>::iterator it;
        qDebug() << "maxIntLong " << maxIntLong;
        qDebug() << "maxIntShort " << maxIntShort;

        for( it = input->begin(); it != input->end(); ++it )
        {
            sdmixer::Localization loc = *it;

            if(it == input->begin())
            {
                maxIntLongFromFile = loc.LongIntensity;
                maxIntShortFromFile = loc.ShortIntensity;
            }
            if (maxIntLongFromFile < loc.LongIntensity)
                maxIntLongFromFile = loc.LongIntensity;
            if(maxIntShortFromFile < loc.ShortIntensity)
                maxIntShortFromFile = loc.ShortIntensity;

            int longVal = loc.LongIntensity;
            int shortVal = loc.ShortIntensity;

            if( longVal >= maxIntLong )
                longVal = maxIntLong;
            if( shortVal >= maxIntShort )
                shortVal = maxIntShort;

            longVal = round(longVal*precision);
            shortVal = round(shortVal*precision);
            int x_val, y_val;

            if (x_Short)
            {
                x_val = shortVal;
                y_val = longVal;
            }
            else
            {
                x_val = longVal;
                y_val = shortVal;
            }
            if(x_val >= IntSpace_width)
                x_val = IntSpace_width-1;
            //if(y_val > IntSpace_height)
              //  y_val = IntSpace_height;

            y_val = (roundedY-y_val);

            //qDebug() << "long: " << loc.LongIntensity << " rounded: " << longVal;
            //qDebug() << "short: " << loc.ShortIntensity << " rounded: " << shortVal;

            QColor clrPixel(img.pixel(x_val, y_val));
            //int ind = shortVal*bin_intSpace + IntSpace_width*(roundedY-longVal)*bin_intSpace;
            int ind = x_val + IntSpace_width*y_val;

            IntensitySpaceArray[ind]+=1;
            //qDebug() << "shortval: "<< shortVal << " longval: "<< longVal<< " ind: " << ind << "IntensitySpaceArray[ind]: " <<IntensitySpaceArray[ind];
            if(ind != IntSpace_width-1)
            {
                if ( IntensitySpaceArray[ind] > max_value_intSpace )//&&
                     //shortVal < (IntSpace_width-2) &&
                     //longVal < ((roundedY-longVal)-2))
                {
                    max_value_intSpace = IntensitySpaceArray[ind];
                    max_value_index=ind;
                }
                if ( IntensitySpaceArray[ind] < min_value_intSpace )
                {
                    min_value_intSpace = IntensitySpaceArray[ind];
                    min_value_index=ind;
                }
            }
            if (clrPixel.black() == 0)
                it->filter = current_filter;
            /*if(it->filter == 0)
                sum_channel_0++;
            if(it->filter == current_filter && it->filter != 0)
            {
                sumChannels[current_filter-1]++;
                sum++;
            }

            ++index;*/
        }
       // pairs_in_channel.push_back(sum);
        ++current_filter;
    }
    qDebug() << "maxIntLongFromFile " << maxIntLongFromFile;
    qDebug() << "maxIntShortFromFile " << maxIntShortFromFile;

    int total_localizations=0;
    int sum_channel_0=0;
    int *sumChannels = new int[FilterInputFiles.size()];
    for (int i = 0; i < FilterInputFiles.size(); ++i)
        sumChannels[i] = 0;

    std::vector<sdmixer::Localization>::iterator it;
    for( it = input->begin(); it != input->end(); ++it )
    {
        sdmixer::Localization loc = *it;
        if(loc.filter!=0)
        {
            sumChannels[loc.filter-1]++;
        }
        else
            sum_channel_0++;
        total_localizations++;
    }
    qDebug() << "total_localizations: " << total_localizations;

    //std::vector<sdmixer::Localization>::iterator it;

    QString channel_results = "";

    int ii=0;
    for ( auto i : FilterInputFiles)
    {
        channel_results="found ";
        //channel_results.append(QString::number(pairs_in_channel[ii]));
        channel_results.append(QString::number(sumChannels[ii]));
        channel_results.append(" (");
        //channel_results.append(QString::number(100*pairs_in_channel[ii]/index));
        channel_results.append(QString::number(100*sumChannels[ii]/total_localizations));
        channel_results.append(" %) pairs in channel ");
        channel_results.append(QString::number(ii+1));
        sdm->writeToLogFile(channel_results);
        qDebug() << channel_results;
        ii++;
    }

    channel_results = "";
    channel_results.append(QString::number(sum_channel_0));
    channel_results.append(" (");
    channel_results.append(QString::number(100*sum_channel_0/total_localizations));
    channel_results.append(" %) pairs crosstalk (channel=0)");
    sdm->writeToLogFile(channel_results);
    qDebug() << channel_results;


    for( it = input->begin(); it != input->end(); ++it )
    {
        sdmixer::Localization loc = *it;
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
        out << loc.LongIntensity << " ";
        out << loc.filter << "\n";
    }
    filter_out.close();

    qDebug() << "min_val intspace" << min_value_intSpace << " @ind: " <<min_value_index;
    qDebug() << "max_val intspace" << max_value_intSpace << " @ind: " <<max_value_index;
    if(sdm->getPlotIntensitySpace())
    {
    int range_x = round(0.001*IntSpace_width+2.2);
    if (range_x < 2)
        range_x = 2;
    if(range_x % 2 != 0)
        range_x++;
    int range_y = range_x;
    int half_range_x = range_x/2;
    int half_range_y = half_range_x;

    qDebug() << "bins, range_x = " << range_x << " half_range: " << half_range_x;

    int minValAfterBin=0;
    int maxValAfterBin=0;


    for ( int i = 0; i < IntSpace_width; i+=range_x)
        for(int j = 0; j < IntSpace_height; j+=range_y)
        {

            int ind = i+j*IntSpace_width;
            int SUM = 0;

            if ( i < half_range_x || i >= IntSpace_width-half_range_x ||
                 j < half_range_y || j >= IntSpace_height-half_range_y )
            {
                SUM=IntensitySpaceArray[ind];
            }
            else
            {
                for (int I = -half_range_x; I < half_range_x; ++I)
                {
                    for(int J = -half_range_y; J < half_range_y; ++J)
                    {
                        ind = (i+I)+(j+J)*IntSpace_width;
                        SUM += IntensitySpaceArray[ind];
                    }
                }

                for (int I = -half_range_x; I < half_range_x; ++I)
                {
                    for(int J = -half_range_y; J < half_range_y; ++J)
                    {
                        ind = (i+I)+(j+J)*IntSpace_width;
                        IntensitySpaceArray[ind] = SUM;
                    }
                }
            }
            if (IntensitySpaceArray[ind] < minValAfterBin)
                minValAfterBin = IntensitySpaceArray[ind];
            if (IntensitySpaceArray[ind] > maxValAfterBin)
                maxValAfterBin = IntensitySpaceArray[ind];
        }
    qDebug() << "minValAfterBin :" << minValAfterBin << " maxValAfterBin: " <<  maxValAfterBin;

    for ( int i = 0; i < IntSpace_width; ++i)
        for(int j = 0; j < IntSpace_height; ++j)
        {
            int ind = i+j*IntSpace_width;
            double d = (double)IntensitySpaceArray[ind];
            double norm = (double)(maxValAfterBin-minValAfterBin);
            double val = d*(1.0/norm);
            if(val-0.0 > 1e-8)
            {
                QColor qc;

                if(val > 1)
                    val = 1;
                if(val < 0)
                    val = 0;

                //val-=1;+
                //qc.setHsvF(0.5, 0.5, val);
                qc.setHslF(val, 1.0, 0.5);

                /*QPainter painter(&IntSpaceOut);
                painter.setPen(qc);
                QPen pen = painter.pen();
                pen.setWidth(range_x);
                painter.setPen(pen);
                painter.drawEllipse(i, j, range_x, range_y);
                painter.end();*/

                IntSpaceOut.setPixel(i, j, qc.rgb());
            }
        }

    qDebug() << "starting edge detection" ;

    for ( ff = FilterInputFiles.begin(); ff < FilterInputFiles.end(); ++ff)
    {
        QString ii = *ff;
        qDebug() << "filter file: " << ii;
        QImage img(ii);
        QImage out(img.width(), img.height(), QImage::Format_RGB32);
        qDebug() << "out.width: "<< out.width() << " height: " <<out.height();

        EdgeDetectionSobel(img, out);

        for(int i = 0; i < IntSpace_width; ++i)
            for(int j = 0; j < IntSpace_height; ++j)
            {
                QColor clrPixel(out.pixel(i, j));
                if(clrPixel.black() != 0)
                {
                    if ( i < half_range_x || i >= IntSpace_width-half_range_x ||
                         j < half_range_y || j >= IntSpace_height-half_range_y )
                    {
                        IntSpaceOut.setPixel(i, j, qRgb(0,0,0));
                    }
                    else
                    {
                        for (int I = -half_range_x; I < half_range_x; ++I)
                        {
                            for(int J = -half_range_y; J < half_range_y; ++J)
                            {
                                IntSpaceOut.setPixel((i+I), (j+J), qRgb(0,0,0));
                            }
                        }

                    }
                }
            }
    }
    qDebug() << "finished edge detection" ;
    IntSpaceOut.save(intensitySpaceFile);
}

    qDebug() <<" finished filter";
    emit finished();
}


// Qt version:
/* Given image in source place Sobel edges in dest.
Grayscale sort of, with (255,255,255) as brightest edge.
sobelDestination should be same size and depth as source.
*/
void Filter::EdgeDetectionSobel(QImage &source,
                        QImage &sobelDestination)
    {
    int GX[3][3];
    int GY[3][3];
    /* 3x3 GX Sobel mask.  Ref: www.cee.hw.ac.uk/hipr/html/sobel.html */
    GX[0][0] = -1; GX[0][1] = 0; GX[0][2] = 1;
    GX[1][0] = -2; GX[1][1] = 0; GX[1][2] = 2;
    GX[2][0] = -1; GX[2][1] = 0; GX[2][2] = 1;

    /* 3x3 GY Sobel mask.  Ref: www.cee.hw.ac.uk/hipr/html/sobel.html */
    GY[0][0] =  1; GY[0][1] =  2; GY[0][2] =  1;
    GY[1][0] =  0; GY[1][1] =  0; GY[1][2] =  0;
    GY[2][0] = -1; GY[2][1] = -2; GY[2][2] = -1;

    int width = source.width();
    int height = source.height();

    qDebug() << "EdgeDetection: " << width << " x " << height;

    int I, J;
    long sumX, sumY;
    int SUM;
    QRgb rawColour;

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            if( y == 0 || y >= height-1 || x == 0 || x >= width-1 )
            {
                SUM = 0;
            }
            else
            {
                sumX = 0;
                sumY = 0;
                /*-------X and Y GRADIENT APPROXIMATION------*/
                for(I=-1; I<=1; I++)
                {
                    for(J=-1; J<=1; J++)
                    {
                        rawColour = source.pixel(x+I, y+J);
                        sumX = sumX + ((qRed(rawColour) + qGreen(rawColour) + qBlue(rawColour))/3) * GX[I+1][J+1];
                        sumY = sumY + ((qRed(rawColour) + qGreen(rawColour) + qBlue(rawColour))/3) * GY[I+1][J+1];
                    }
                }
                SUM = abs(sumX) + abs(sumY); /*---GRADIENT MAGNITUDE APPROXIMATION (Myler p.218)----*/
                //qDebug() << "rawColour: " << rawColour << " qRed " << qRed(rawColour) << " qGreen " << qGreen(rawColour) << " qBlue " << qBlue(rawColour) << " SUM " << SUM;
                if (SUM > 255)
                    SUM = 255;

                }
            //qDebug() << "x: " << x << " y: " << y;
            //sobelDestination.setPixel(x, y, qRgb(0, 0, 0));
            sobelDestination.setPixel(x, y, qRgb(SUM, SUM, SUM));
            /**dp = qRgb(SUM, SUM, SUM);
            ++p;
            ++dp;*/
            }
    }
    sobelDestination.invertPixels();
}

