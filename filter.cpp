#include "filter.h"
#include "pairfinder.h"
#include <QDebug>

Filter::Filter()
{
    // get input file
    // load filters
    // scale down intensities
    // foreach filter check if localization is inside filter
    // --> write filter number

}

void Filter::loadFilterImage(QString path)
{
    //QImgFilters[0].load(path);

    QImage img(path);
    int arr[img.height()][img.width()];
    std::vector< std::vector < int > > filterImg;

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
    }

    std::ofstream out("io3.txt");
    /*for(int i = 0; i < img.height(); ++i)
    {
        for(int j = 0; j < img.width(); ++j)
        {
            out << arr[i][j];
            out << "  ";
        }
    out << std::endl;
    }*/
    for(auto i = filterImg.begin(); i != filterImg.end(); ++i)
    {
        for(auto j = i->begin(); j != i->end(); ++j)
        {
            out << *j;
            out << "  ";
        }
        out << std::endl;
    }



}

void Filter::initializeIntensities()
{
    std::vector<PairFinder::Localization> PFOutput;

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


