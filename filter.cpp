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

    fileName = fi.baseName();
    if(fileName.contains("_pairs_out"))
        fileName.replace("_pairs_out", "");
    if(fileName.contains("_filter_out"))
        fileName.replace("_filter_out", "");

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

    while (!in.atEnd())
    {
        line = in.readLine();

        std::vector<std::string> v;
        std::vector<double> dbl_vec;
        std::string strLine = line.toStdString();

        boost::split(v, strLine, boost::is_any_of("\t "));

        int index = 0;
        for (auto i : v)
        {
            index++;
            dbl_vec.push_back(strtod(i.c_str(), NULL));
        }

        sdmixer::Localization loc;
        loc.xShort = dbl_vec[0];
        loc.yShort = dbl_vec[1];
        loc.zShort = dbl_vec[2];
        loc.ShortIntensity = dbl_vec[3];
        loc.frame = dbl_vec[4];
        loc.xLong = dbl_vec[5];
        loc.yLong =dbl_vec[6];
        loc.zLong = dbl_vec[7];
        loc.LongIntensity = dbl_vec[8];
        //qDebug() << loc.xShort;
        //qDebug() << loc.LongIntensity;


        sdm->pushBackLocalization(loc);
    }
    input = sdm->getPfOutput();
    qDebug() << "loaded ";


}

Filter::Filter(sdmixer *s, QString str)
{
    qDebug() << "inititalizing Filter: " << str;
    this->sdm = s;
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
    int roundedY = round(maxIntLong*precision);
    qDebug() << roundedY;

    QFile filter_out(outputFile);
    QTextStream out(&filter_out);
    filter_out.open(QIODevice::WriteOnly | QIODevice::Text);
    PairFinder p(sdm, sdm->getCurrentFile());
    p.writeHeader(out, true);

    int IntSpace_width = round(maxIntShort*precision);
    int IntSpace_height = round(maxIntLong*precision);

    uint32_t IntensitySpaceNPixel = IntSpace_width*IntSpace_height;
    double bin_intSpace = (maxIntShort*precision)/IntSpace_width;
    qDebug() << "bin int space :" << bin_intSpace;

    float *IntensitySpaceArray = new float[IntensitySpaceNPixel];

    for ( int i = 0; i < IntensitySpaceNPixel; ++i)
        IntensitySpaceArray[i]=0;

    int current_filter=1;
    std::vector<QString>::iterator ff;
    float max_value_intSpace = 0;
    int max_value_index=0;
    QImage IntSpaceOut(IntSpace_width, IntSpace_height, QImage::Format_RGB32);
    IntSpaceOut.fill(Qt::white);
    qDebug() << "average";
    for ( ff = FilterInputFiles.begin(); ff < FilterInputFiles.end(); ++ff)
    {
        QString ii = *ff;
        QImage img(ii);
        for(int i = 0; i < IntSpace_width; ++i)
            for(int j = 0; j < IntSpace_height; ++j)
            {
                QColor clrPixel(img.pixel(i, j));
                QColor clrPainter(IntSpaceOut.pixel(i, j));

                int average = ((clrPixel.red()+clrPainter.red())/2+
                              (clrPixel.green()+clrPainter.green())/2+
                              (clrPixel.blue()+clrPainter.blue())/2);

                IntSpaceOut.setPixel(i, j, qRgb(average, average, average));
            }
    }
     qDebug() << "average ready";

    int index = 0;

    std::vector<int> pairs_in_channel;


    int sum_channel_0=0;


    for ( ff = FilterInputFiles.begin(); ff < FilterInputFiles.end(); ++ff)
    {

        QString i = *ff;

        QImage img(i);
        qDebug() << i;

        int sum = 0;



        std::vector<sdmixer::Localization>::iterator it;

        for( it = input->begin(); it != input->end(); ++it )
        {
            sdmixer::Localization loc = *it;

            if(index == 0)
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

            /*if(longVal >= maxIntLong*precision)
                longVal--;*/
            if(shortVal >= maxIntShort*precision)
                shortVal--;


            QColor clrPixel(img.pixel(shortVal, (roundedY-longVal)));
            //int ind = shortVal*bin_intSpace + IntSpace_width*(roundedY-longVal)*bin_intSpace;
            int ind = shortVal + IntSpace_width*((roundedY-longVal));

            IntensitySpaceArray[ind]+=1;
            //qDebug() << "shortval: "<< shortVal << " longval: "<< longVal<< " ind: " << ind << "IntensitySpaceArray[ind]: " <<IntensitySpaceArray[ind];

            if ( IntensitySpaceArray[ind] > max_value_intSpace &&
                 shortVal < (IntSpace_width-2) &&
                 longVal < ((roundedY-longVal)-2))
            {
                max_value_intSpace = IntensitySpaceArray[ind];
                max_value_index=ind;
            }

            if (clrPixel.black() == 0)
                it->filter= current_filter;
            if(it->filter == 0)
                sum_channel_0++;
            else
                sum++;


            //qDebug() << shortVal << "  " << longVal << "  " << clrPixel.black() << "  " << it->filter;

            ++index;
        }
        pairs_in_channel.push_back(sum);
        ++current_filter;
    }
    std::vector<sdmixer::Localization>::iterator it;

    QString channel_results = "";

    int ii=0;
    for ( auto i : FilterInputFiles)
    {
        channel_results="found ";
        channel_results.append(QString::number(pairs_in_channel[ii]));
        channel_results.append(" (");
        channel_results.append(QString::number(100*pairs_in_channel[ii]/index));
        channel_results.append(" %) pairs in channel ");
        channel_results.append(QString::number(ii+1));
        sdm->writeToLogFile(channel_results);
        ii++;
    }

    channel_results = "";
    channel_results.append(QString::number(sum_channel_0));
    channel_results.append(" (");
    channel_results.append(QString::number(100*sum_channel_0/index));
    channel_results.append(" %) pairs crosstalk (channel=0)");
    sdm->writeToLogFile(channel_results);


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
    qDebug() << " max_val intspace" << max_value_intSpace << " @ind: " <<max_value_index;


    for ( int i = 0; i < IntSpace_width; ++i)
        for(int j = 0; j < IntSpace_height; ++j)
        {
            QColor qc;
            int ind = i+j*IntSpace_width;

            double val = (( IntensitySpaceArray[ind])*(1.0/(max_value_intSpace)));
            if(val-0.0 > 1e-5)
            {
                if(val > 1)
                    val= 1;

                qc.setHsvF(0.5, 0.5, val);

            IntSpaceOut.setPixel(i, j, qc.rgb());
            }
        }


    IntSpaceOut.save(intensitySpaceFile);

    qDebug() <<" finished filter";
    emit finished();
}


