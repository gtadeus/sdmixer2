#include "filter.h"
//#include "pairfinder.h"
#include <QDebug>

Filter::Filter(sdmixer *s, std::vector<sdmixer::Localization> *data)
{
    qDebug() << "inititalizing Filter";
    this->sdm = s;
    input=data;
    precision =sdm->getPrecision();
    maxIntLong = sdm->getMaxIntLong();
    maxIntShort = sdm->getMaxIntShort();




}
Filter::Filter(sdmixer *s, QString str)
{
    qDebug() << "inititalizing Filter: " << str;
    this->sdm = s;

}

void Filter::doWork()
{
    qDebug() << "starting Filter in new Thread";
    FilterInputFiles = sdm->getFilterFiles();


    for ( auto i : FilterInputFiles)
    {
        //QImgFilters.push_back(q);
        int index = 0;
        QImage img(i);
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

        int longVal = round(loc.LongIntensity*precision);
        int shortVal = round(loc.ShortIntensity*precision);
        if( longVal >= maxIntLong )
            longVal = maxIntLong;
        if( shortVal >= maxIntShort )
            shortVal = maxIntShort;

        int current_filter=0;

            ++current_filter;

            QColor clrPixel(img.pixel(shortVal, longVal));

            if (clrPixel.black() != 0)
                loc.filter=current_filter;

            //qDebug() << shortVal << "  " << longVal << "  " << clrPixel.black() << "  " << loc.filter;


        ++index;
    }
}
    qDebug() <<" finished filter";
    emit finished();
}

QImage Filter::loadFilterImage(QString path)
{

    QImage img(path);
    /*std::vector< std::vector < int > > filterImg;

    if ( !img.isNull())
    {
        for(int row = 0; row < img.height(); ++row)
        {
            std::vector<int> row_vec;
            for(int col = 0; col < img.width(); ++col)
            {
                QColor clrPixel(img.pixel(row, col));
                //arr[row][col] = clrPixel.black();
                row_vec.push_back(clrPixel.black());
            }
            filterImg.push_back(row_vec);
        }
    }*/
    return img;
}

void Filter::initializeIntensities()
{
    std::vector<sdmixer::Localization> PFOutput;

    for (auto &i : PFOutput )
    {
        //round Intesity values to Integers
        i.LongIntensity = round(i.LongIntensity);
        i.ShortIntensity = round(i.ShortIntensity);


        // cut off if Intensity values are > maxInt from config
       if( i.LongIntensity > maxIntLong)
           i.LongIntensity = maxIntLong;

       if( i.ShortIntensity > maxIntShort)
           i.ShortIntensity = maxIntShort;

       // divide Intensity values by 10 for example.
       // so we get a 8000x3000 px Image instead of 80000x30000
       i.LongIntensity *= precision;
       i.ShortIntensity *= precision;



    }
}


void Filter::run()
{

    loadFilterImage("test.png");
}


