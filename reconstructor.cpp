#include "reconstructor.h"
#include "sdmixer.h"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

bool Reconstructor::file_exists (char *filename)
{
    if (FILE *file = fopen(filename, "r"))
    {
        fclose(file);
        return true;
    }
    else
        return false;

}

void Reconstructor::Convolution2()
{
    qDebug() << " started convolution 2";
    uint64_t img_sizeX = image_size[0];
    uint64_t img_sizeY = image_size[1];
    uint64_t img_sizeZ;
    if(dimensions==3 && !force2D)
        img_sizeZ = image_size[2];
    else
        img_sizeZ = 1;

    uint64_t M=img_sizeX;
    uint64_t N=img_sizeY;
    uint64_t Z=img_sizeZ;

    uint64_t NPixel=img_sizeX*img_sizeY*img_sizeZ;

    boost::iostreams::mapped_file mmap_image_in;
    boost::iostreams::mapped_file_params params_in;
    params_in.mode = (std::ios::out | std::ios::in);
    params_in.path = tiff_temp_file;
    mmap_image_in.open(params_in);
    uint16 *image_in =static_cast<uint16*>((void*)mmap_image_in.data());

    boost::iostreams::mapped_file conv_out_file;
    boost::iostreams::mapped_file_params params_out;
    if(file_exists(convolved_image))
        remove(convolved_image);
    params_out.path = convolved_image;
    params_out.mode = (std::ios::out | std::ios::in);
    params_out.new_file_size =  NPixel*sizeof(uint16);
    conv_out_file.open(params_out);
    uint16 *image_out = static_cast<uint16*>((void*)conv_out_file.data());


    qDebug() << "kernel_:size: " << krn.size;
    qDebug() << "kernel_3D: " << krn.make3D;

    std::vector<double> sum_vec;
    std::vector<ConvPixel> conv_res_vec;


    if(dimensions == 3 && !force2D)
    {
        for(uint64_t i = 0; i < M; ++i)
        {
            for ( uint64_t j = 0 ; j < N ; ++j)
            {
                for ( uint64_t k = 0; k < Z; ++k)
                {
                    float sum = 0;

                    for( uint64_t m = 0; m < krn.size; ++m)
                    {
                        uint64_t mm = krn.size - 1 - m;
                        for( uint64_t n = 0; n < krn.size; ++n)
                        {
                            uint64_t nn = krn.size - 1 - n;
                            for ( uint64_t p = 0; p < krn.size; ++p)
                            {
                                uint64_t pp = krn.size - 1 - p;

                            uint64_t rowIndex = j + m - krn.k;
                            uint64_t colIndex = i + n - krn.k;
                            uint64_t zIndex = k + p - krn.k;


                            if(rowIndex >=0 && rowIndex < N &&
                                    colIndex >= 0 && colIndex < M &&
                                    zIndex >= 0 && zIndex < Z)
                            {

                                sum+=(double)image_in[M*N*zIndex+M*rowIndex+colIndex]*
                                        krn.data[krn.size*krn.size*pp+krn.size*mm + nn];

                            }
                            }
                        }
                    }
                    ConvPixel cp;

                    //int val = round(sum*(65536/255));
                    double val = sum;
                    if(val>65535)
                        val=65535;
                    //image_out[M*N*k+M*j+i] = val;
                    cp.index = M*N*k+M*j+i;
                    cp.sum = val;
                    conv_res_vec.push_back(cp);
                }
            }
        }
    }
    else
    {
        for(uint64_t i = 0; i < M; ++i)
        {
            for ( uint64_t j = 0 ; j < N ; ++j)
            {
                float sum = 0;
                for( uint64_t m = 0; m < krn.size; ++m)
                {
                    uint64_t mm = krn.size - 1 - m;
                    for( uint64_t n = 0; n < krn.size; ++n)
                    {
                        uint64_t nn = krn.size - 1 - n;

                        uint64_t rowIndex = j + m - krn.k;
                        uint64_t colIndex = i + n - krn.k;

                        if(rowIndex >=0 && rowIndex < N &&
                                colIndex >= 0 && colIndex < M)
                        {
                            sum+=(double)image_in[M*rowIndex+colIndex]*
                                    krn.data[krn.size*mm + nn];
                        }
                    }
                }
                ConvPixel cp;
                //int val = round(sum*(65536/255));
                double val = sum;
                //sum_vec.push_back(sum);
                if(val>65535)
                    val=65535;
                //image_out[M*j+i] = val;
                cp.index = M*j+i;
                cp.sum = val;
                conv_res_vec.push_back(cp);
            }
        }
    }

    std::vector<ConvPixel>::iterator it;
    double max_sum = 0;
    for(it = conv_res_vec.begin(); it<conv_res_vec.end(); ++it)
    {
        ConvPixel c = *it;
        if(c.sum>max_sum)
            max_sum=c.sum;
    }
    qDebug() << "max sum: " << max_sum;

    for(it = conv_res_vec.begin(); it<conv_res_vec.end(); ++it)
    {
        ConvPixel c = *it;
        image_out[c.index] = round(c.sum*(65535.0/max_sum));
        //qDebug() << c.index << "  " << c.sum << "  " << image_out[c.index];
    }

    mmap_image_in.close();
    conv_out_file.close();
    qDebug() << "convolution ready";

}
void Reconstructor::hist_eq()
{
    qDebug() << "start hist_eq";
    boost::iostreams::mapped_file_sink file;
    boost::iostreams::mapped_file_params params;
    if(runConvolution)
        params.path = convolved_image;
    else
        params.path = tiff_temp_file;

    file.open(params);

    if(!file.is_open())
    {
        qDebug() << "no temp file to create OutputTIFF";
        return;
    }
    uint16 *image_out =static_cast<uint16*>((void*)file.data());

    double histo_correct = hist_correct_value;
    double histo_treshold = hist_threshold;
    bool square_eq = sqrtCum;

    double NPixel=maxPixels;
    int bit = 16;
    const int mp=pow(2, bit)-1;


    double *hist = new double[mp+1];
    double *PDFwt = new double[mp+1];
    double *cdf = new double[mp+1];

    for ( int i = 0; i< mp+1; ++i)
    {
        hist[i] = 0;
        PDFwt[i] = 0;
        cdf[i] = 0;
    }
    for (int i = 0; i < NPixel; ++i)
    {
        int index = image_out[i];
        if(index < mp)
            hist[index]+=1;
    }
    PDFwt[0]=(double)hist[0];
    PDFwt[mp]=(double)hist[mp];
    double Pmax = hist[0]/NPixel;
    double Pl = hist[0]/NPixel;

    for (int g=0; g < mp; ++g)              //calculate normal PDF [0,1]
    {
        PDFwt[g]=hist[g]/NPixel;
        if(PDFwt[g]<Pl)
            Pl = hist[g]/NPixel;
        if(PDFwt[g]>Pmax)
            Pmax = hist[g]/NPixel;
    }
    double Pu = histo_treshold * Pmax;                 //thresholding the PDF
    for (int g=0; g < mp; ++g)              //calculate weightet PDF
    {
        PDFwt[g]=pow( (PDFwt[g]-Pl)/(Pu-Pl), histo_correct) * Pu;
        PDFwt[g]*=NPixel;                   //transform to histogram range [0,NPixel]
    }

    if (square_eq)
    {
        cdf[0]=sqrt(PDFwt[0]);
        for (long b=1; b<(mp+1); ++b)   //from PDF to CDF with square root of hist values
        {
            cdf[b]=cdf[b-1]+sqrt(PDFwt[b]);
        }
    }
    else
    {
        cdf[0]=PDFwt[0];
        for (long b=1; b<(mp+1); ++b)   //from PDF to CDF, classic linear approach
        {
            cdf[b]=cdf[b-1]+PDFwt[b];
        }
    }

    for (int i = 0; i < NPixel; ++i)         //for each pixel remap
    {
        int index = image_out[i];
        int val;
        if(index < mp)
            val = round(  mp * (cdf[index]-cdf[0]) / (cdf[mp]-cdf[0]));
        else
            val=mp;

        image_out[i] = val;
    }

    file.close();
    qDebug() << "end hist_eq";
}

void Reconstructor::doWork()
{
    qDebug() <<" started Reconstructor worker";

    if(doWorkLater)
        doWorkNow();

    bool repeat =true;

    if(tiff_out_file.isEmpty())
    {
        input_file = sdm->getCurrentFile();
        QFile qf(input_file);
        QFileInfo fi(qf);

        if(!sdm->getOutputDirectory().isEmpty())
            output_dir = sdm->getOutputDirectory();
        else
            output_dir=fi.path();

        input_base_name = fi.baseName();
        if(input_base_name.contains("_pairs_out"))
            input_base_name.replace("_pairs_out", "");
        if(input_base_name.contains("_filter_out"))
            input_base_name.replace("_filter_out", "");

        tiff_out_file = output_dir.append("/");
        tiff_out_file = tiff_out_file.append(input_base_name);
        if(FilterInput)
        {
            tiff_out_file = tiff_out_file.append("_ch");
            tiff_out_file = tiff_out_file.append(QString::number(curr_filter));
        }

        tiff_out_file.append(".tif");
    }

    while(repeat && curr_filter<=max_filter)
    {
        qDebug() << "curr_filter: " << curr_filter;
        repeat=false;
        setArray();
        sdm->writeToLogFile("initialized image data");
        if(runConvolution)
        {
            sdm->writeToLogFile("starting convolution");

            Convolution2();
            sdm->writeToLogFile("finished convolution");
        }
        if(perform_hist_eq)
        {
            hist_eq();
        }
            //Convolution();
        //map8bit();
        sdm->writeToLogFile("starting TIFF Output");
        outputTIFF();
        sdm->writeToLogFile("finished TIFF Output");
        if(FilterInput)
        {
            repeat=true;
            curr_filter++;
            initData(curr_filter);
        }
    }

    if(file_exists(tiff_temp_file))
        remove(tiff_temp_file);
    if(file_exists(fft_src_file))
        remove(fft_src_file);
    if(file_exists(fft_krnl_file))
        remove(fft_krnl_file);
    if(file_exists(convolved_image))
        remove(convolved_image);

    emit finished();
}


void Reconstructor::setArray()
{
    qDebug() << "start set array";
    std::stringstream ss;
    ss << "image_size (px) : ";
    for(int i = 0; i < dimensions; ++i)
    {
        maxPixels*=image_size[i];

        ss << image_size[i] << "  ";
        qDebug()<<"set Array image_size[i]: " << image_size[i];
        //maxPixels*=(image_max[i]+1);
    }
    sdm->writeToLogFile(QString::fromStdString(ss.str()));

    std::vector<Coordinates>::iterator it;

    int counter=0;
    for( it = xyz.begin(); it != xyz.end(); ++it )
    {
        Coordinates c = *it;
        for(int i = 0; i < dimensions; ++i)
        {
            //qDebug() << c.get(i);
            it->set(i, c.get(i)-image_min[i]);
            if(it->get(i)>=image_size[i])
            {
                counter++;
                xyz.erase(it);
                --it;
            }
            //qDebug() << c.get(i);
        }

    }
    qDebug() << "deleted " << counter << " coordinates";
    qDebug() << "maxPixels: " << maxPixels;
    boost::iostreams::mapped_file file;
    boost::iostreams::mapped_file_params params;

    if(file_exists(tiff_temp_file))
        remove(tiff_temp_file);
    params.path = tiff_temp_file;

    params.mode = (std::ios::out | std::ios::in);
    params.new_file_size =  maxPixels*sizeof(uint16);

    file.open(params);

    uint16 *array =static_cast<uint16*>((void*)file.data());

    if (file.is_open())
    {
        qDebug() << "initializing... max_pixel = " << maxPixels;
       for ( int i = 0; i < maxPixels; ++i)
       {
           array[i]=0;
       }
       qDebug() << "populating array...";
       uint64_t temp_max = 0;
       std::vector<Coordinates>::iterator it;
       std::stringstream ssb;
       ssb << "found " << xyz.size() << " coordinates";
       sdm->writeToLogFile(QString::fromStdString(ssb.str()));

       qDebug()<< "found " << xyz.size() << " coordinates";

       for( it = xyz.begin(); it < xyz.end(); ++it)
       {
           Coordinates c = *it;
           uint64_t lin_index = linearIndex(c);
           //qDebug() << lin_index << "  " << array[lin_index];
           if(lin_index < maxPixels)
           {
               if(array[lin_index]<65535)
                   array[lin_index]+=1;
               else
                   array[lin_index] = 65535;
           }

           if (it == xyz.begin())
           {
               temp_max = lin_index;
               dbl_image_max=array[lin_index];
               dbl_image_min=array[lin_index];
           }
           else{
               if(lin_index > temp_max)
                   temp_max = lin_index;
               if(dbl_image_max<array[lin_index])
                   dbl_image_max=array[lin_index];
               if(dbl_image_min>array[lin_index])
                   dbl_image_min=array[lin_index];
           }
       }
       qDebug() << "max array index from data: " << temp_max;
       qDebug() << "array init ready!";

    }
    file.close();

    qDebug() << "finish array";
}


void Reconstructor::map8bit()
{
    qDebug() << "map8 bit";

    boost::iostreams::mapped_file_sink image;
    boost::iostreams::mapped_file_params params;
    params.path = tiff_temp_file;
    image.open(params);
    uint8 *img_in =static_cast<uint8*>((void*)image.data());

/*
    uint64_t NPixel = 1;
    for (int i = 0; i < dimensions; ++i)
        NPixel*=image_size[i];*/

    //params_file_out.new_file_size =  maxPixels*sizeof(uint8);


    for(uint64_t i = 0; i < maxPixels; ++i)
    {
        img_in[i] =
                (uint8_t) ( (double) img_in[i] - dbl_image_min)*
                (255.0/(dbl_image_max-dbl_image_min));
    }
    image.close();
    qDebug() << "map8 bit finished";

}

void Reconstructor::CreateGaussianKernel()
{
    float f;
    int index = 0;
    float sum = 0;
    float e = 2.71828;

    // kernel[2k+1][2k+1][2k+1]

    if (krn.make3D)
    {
        krn.size_z=krn.size;
        for(int i = -krn.k; i <= krn.k; ++i)
        {
            for(int j = -krn.k; j <= krn.k ; ++j)
            {
                for(int k = -krn.k; k <= krn.k; ++k)
                {
                    f = (float)( i*i/(2.0f * krn.sigma_xy * krn.sigma_xy) +
                                 j*j/(2.0f * krn.sigma_xy * krn.sigma_xy) +
                                 k*k/(2.0f * krn.sigma_z * krn.sigma_z));
                    krn.data[index] = pow(e,-f);
                    sum += krn.data[index];
                    ++index;
                }
            }
        }
    }
    else
    {
        for(int i = -krn.k; i <= krn.k; ++i)
        {
            for(int j = -krn.k; j <= krn.k ; ++j)
            {
                f = (float)( i*i/(2.0f * krn.sigma_xy * krn.sigma_xy) +
                             j*j/(2.0f * krn.sigma_xy * krn.sigma_xy));
                krn.data[index] = pow(e,-f);
                sum += krn.data[index];
                ++index;

            }
        }
    }

    for(int i = 0; i < index; ++i)
        krn.data[i] = krn.data[i]/sum;

}
void Reconstructor::setKernel()
{
    krn.sigma_xy=2;
    krn.sigma_z=2;
    krn.make3D=true ;
    krn.data = new float[krn.size*krn.size*krn.size];
    CreateGaussianKernel();
}

uint64_t Reconstructor::linearIndex(Coordinates c)
{
    //qDebug() << "linIndex of: " << c.get(0) << "  " << c.get(1) << " " << c.get(2);
    //qDebug() << (uint64_t) c.get(0) + (uint64_t) (image_size[0])*( c.get(1) + (uint64_t) (image_size[1])*c.get(2));
    return (uint64_t) c.get(0) + (uint64_t) (image_size[0])*( c.get(1) + (uint64_t) (image_size[1])*c.get(2));
    //return c.get(2) + (uint64_t) (image_max[2]+1)*( c.get(1) + (uint64_t) (image_max[1]+1)*c.get(0));
    //return c.get(0) + (uint64_t) (image_max[0]+1)*c.get(1) + (uint64_t) (image_max[1]+1)*(image_max[0]+1)*c.get(2);
    //return image_max[0]*image_max[1]*c.get(0)+image_max[1]*c.get(1)+c.get(2);

}

uint64_t Reconstructor::linearIndex3DFFT( int i, int j, int k, int z_fftw, int h_fftw)
{
    return (uint64_t) k + ((uint64_t)(z_fftw+2))*((uint64_t)(j+h_fftw*i));
}
uint64_t Reconstructor::linearIndex2DFFT( int i, int j, int h_fftw)
{
    return (uint64_t)j+((uint64_t)(h_fftw+2))*((uint64_t)i);
}
uint64_t Reconstructor::linearIndexTIFF( int i, int j, int k, int img_sizeX, int img_sizeY)
{
    return (uint64_t)i+((uint64_t)img_sizeX)*((uint64_t) j+ ((uint64_t)img_sizeY)*((uint64_t)k) );
}


void Reconstructor::setMinMax(sdmixer::min_max m)
{
    // min max, convert from m to nm
    int min_val[3]={0};
    int max_val[3]={0};

    qDebug() << m.getMax(0);

    min_val[0]=round(m.min_x * 1e9);
    max_val[0]=round(m.max_x * 1e9);

    min_val[1]=round(m.min_y * 1e9);
    max_val[1]=round(m.max_y * 1e9);

    min_val[2]=round(m.min_z * 1e9);
    max_val[2]=round(m.max_z * 1e9);
    qDebug() << min_val[0] << "  " << max_val[0];
    qDebug() << min_val[1] << "  " << max_val[1];
    qDebug() << min_val[2] << "  " << max_val[2];

    qDebug() << " image max/min: ";
    qDebug() << " dimensions " << dimensions;

    for (int i = 0; i < dimensions; ++i)
    {
        qDebug() << i;
      if(max_val[i] != 0)
      {
          //qDebug() << "binning: " << binning[i];
          image_max[i] = round(max_val[i]/binning[i]);
          image_min[i] = round(min_val[i]/binning[i]);
          image_size[i] = (image_max[i]-image_min[i]) + 1;
          qDebug() << "Reconstructor setMinMax min: " << image_min[i] << "  " << image_max[i] << " size: " <<image_size[i];
      }



    }
    minMaxDefined=true;
    qDebug() << "reconstructor min max ready";
    //image_max[i] = round(max_val[i])/xy_binning;
    //image_max[2] = round(max_val[2])/z_binning;
}


void Reconstructor::getHeader(QString header)
{
    header = header.replace("#", "");

    QDomDocument qd;
    qd.setContent(header);

    QDomElement element = qd.documentElement();

    if(element.tagName().contains("localizations"))
        xyzFile = true;
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
                    min_maxValues.min_y = e.attribute("min").replace("m", "").toDouble();
                    min_maxValues.max_y = e.attribute("max").replace("m", "").toDouble();
                }
                else if(e.attribute("identifier").contains("2"))
                {
                    min_maxValues.min_z = e.attribute("min").replace("m", "").toDouble();
                    min_maxValues.max_z = e.attribute("max").replace("m", "").toDouble();
                }
                else
                {
                    min_maxValues.min_x = e.attribute("min").replace("m", "").toDouble();
                    min_maxValues.max_x = e.attribute("max").replace("m", "").toDouble();
                }
            }
            if( e.attribute("identifier").contains("Filter"))
            {
                FilterInput = true;
            }
        }

    }
    if(dimensions > 3)
        dimensions/=2;

    dimDefined=true;
    if(force2D)
        dimensions = 2;
    qDebug() << "FilterInput: " << FilterInput;
    qDebug() << "xyzFile: " << xyzFile;
    qDebug() << "dimensions: " << dimensions;
    qDebug() << "max Values from config";
    qDebug() << min_maxValues.min_x << "  " << min_maxValues.max_x;
    qDebug() << min_maxValues.min_y << "  " << min_maxValues.max_y;
    qDebug() << min_maxValues.min_z << "  " << min_maxValues.max_z;

    //if(FilterInput)
    {
        //dimensions--;
        setMinMax(min_maxValues);
    }

    minMaxDefined=true;

}
void Reconstructor::getSettingsFromGUI()
{
    maxPixels=1;
    binning[0] = sdm->getReconstructor_xyBinning();
    binning[1] = sdm->getReconstructor_xyBinning();
    binning[2] = sdm->getReconstructor_zBinning();
    if(!dimDefined)
    {
        dimensions = sdm->getCurrentDimensions();
        qDebug()<<"dimensions_ " << dimensions;
    }
    qDebug() << "minMaxDefined: " << minMaxDefined;
    if(!minMaxDefined)
    {
        //qDebug() << (sdm->getPF_min_maxValues()).;
        setMinMax(sdm->getPF_min_maxValues());

    }
    else
        setMinMax(min_maxValues);

    runConvolution=sdm->getRunConvolution();

    qDebug() << "get ConvKernel";
    oneConvolutionKernel = sdm->getOneKernelForAllChannels();
    qDebug() <<oneConvolutionKernel;
    NM_PER_PX = sdm->getPixelSize();
    force2D = sdm->getForce2D();

    perform_hist_eq = sdm->getNonLinearHistEq();
    sqrtCum = sdm->getSqrtCummulation();
    hist_correct_value = sdm->getHisteqCoefficient();
    hist_threshold = sdm->getThreshold();


    if(oneConvolutionKernel)
    {
        globalKernel = sdm->getGlobalKernel();
        if(!globalKernel.unitFWHM_xy.compare("px"))
        {
            krn.sigma_xy=globalKernel.FWHM_xy*NM_PER_PX / 2.355;
            krn.sigma_z=globalKernel.FWHM_z*NM_PER_PX / 2.355;
        }
        else
        {
            krn.sigma_xy=globalKernel.FWHM_xy / 2.355;
            krn.sigma_z=globalKernel.FWHM_xy / 2.355;
        }
        if(dimensions==3 && !force2D)
        {
            krn.make3D=true;
            krn.data = new float[krn.size*krn.size*krn.size];
        }
        else
        {
            krn.make3D=false;
            krn.data = new float[krn.size*krn.size];
        }
        CreateGaussianKernel();
        //qDebug() << krn.sigma_xy;
    }
    else
    {
        vecKernel = sdm->getConvolutionKernel();
        //if(getKernelVector)
        /*{
            vecKernel = sdm->getConvolutionKernel();
            if(!vecKernel.empty())
            {
                sdmixer::gaussian_kernel gk = vecKernel[current_filter];
                if(!gk.unitFWHM_xy.compare("px"))
                {
                    krn.sigma_xy=gk.FWHM_xy*NM_PER_PX / 2.355;
                    krn.sigma_z=gk.FWHM_z*NM_PER_PX / 2.355;
                }
                else
                {
                    krn.sigma_xy=gk.FWHM_xy / 2.355;
                    krn.sigma_z=gk.FWHM_z / 2.355;

                }
                if(dimensions==3 && !force2D)
                {
                    krn.make3D=true;
                    krn.data = new float[krn.size*krn.size*krn.size];
                }
                else
                {
                    krn.make3D=false;
                    krn.data = new float[krn.size*krn.size];
                }
                CreateGaussianKernel();
            }
        }*/
    }

    qDebug()<< "reconstructor init ready";



}
void Reconstructor::doWorkNow()
{
    QString xyz_file=xyz_file_parameter;

    getSettingsFromGUI();
    QFile f(xyz_file);
    f.open(QIODevice::ReadOnly| QIODevice::Text);

    QTextStream in(&f);
    QString line = in.readLine();

    getHeader(line);

    if( ! xyzFile )
    {
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
            //qDebug() << "v.size() " << v.size();


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
            if(v.size()>9)
            {
                loc.filter = dbl_vec[9];
                if(loc.filter>max_filter)
                    max_filter=loc.filter;
            }

            //qDebug() << loc.xShort;
            //qDebug() << loc.LongIntensity;


            sdm->pushBackLocalization(loc);
        }
        input_data = sdm->getPfOutput();
        qDebug() << "loaded ";

        qDebug() << " is this a filter_file?: " << FilterInput;
        if(FilterInput)
        {
            initData(1);

        }
        else
            initData(0);
    }
    else // this is xyz File
    {
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
            //qDebug() << "v.size() " << v.size();


            Coordinates c;
            for(int i = 0; i < dimensions; ++i)
            {
                c.set(i, round(dbl_vec[i]/binning[i]));
            }
            xyz.push_back(c);
        }
        // get File Name
        input_file = sdm->getCurrentFile();

        QFile qf(input_file);
        QFileInfo fi(qf);

        if(!sdm->getOutputDirectory().isEmpty())
            output_dir = sdm->getOutputDirectory();
        else
            output_dir=fi.path();

        input_base_name = fi.baseName();

        tiff_out_file = output_dir.append("/");
        tiff_out_file = tiff_out_file.append(input_base_name);

        tiff_out_file.append(".tif");

        getSettingsFromGUI();
    }

}

Reconstructor::Reconstructor(sdmixer *s, QString xyz_file)
{
    qDebug() << "Reconstructor loading file " << xyz_file;
    this->sdm=s;
    xyz_file_parameter=xyz_file;
    doWorkLater=true;

}
void Reconstructor::initData(int current_filter)
{
    if(force2D)
    {
        dimensions = 2;
        image_size[2]=0;
    }

    int limage_min[3]={0};
    int limage_max[3]={0};
    int limage_size[3]={0};

    int index = 0;
    std::vector<sdmixer::Localization>::iterator it;

    for( it = input_data->begin(); it != input_data->end(); ++it )
    {
        sdmixer::Localization loc = *it;
        Coordinates c;
        if(loc.filter == current_filter)
        {
            for(int i = 0; i < dimensions; ++i)
            {
                int short_val = round(loc.getShortDim(i)/binning[i]);
                if (index == 0)
                {
                    limage_min[i]=short_val;
                    limage_max[i]=short_val;
                }
                c.set(i, short_val);
                if(limage_min[i]>short_val)
                        limage_min[i]=short_val;
                if(limage_max[i]<short_val)
                    limage_max[i]=short_val;

                ++index;
            }
            xyz.push_back(c);
            for(int i = 0; i < dimensions; ++i)
            {
                int long_val = round(loc.getLongDim(i)/binning[i]);
                if(limage_min[i]>long_val)
                        limage_min[i]=long_val;
                if(limage_max[i]<long_val)
                    limage_max[i]=long_val;

                c.set(i, long_val);
            }
            xyz.push_back(c);
        }
    }

    for ( int i = 0; i < dimensions; ++i)
    {
        limage_size[i] = (limage_max[i] - limage_min[i]) + 1;
qDebug() << "Reconstructor (image size, from input data): " << limage_size[i];
qDebug() << "Reconstructor (image size, from minmax): " << image_size[i];
    }
    if(!minMaxDefined)
    {
        for ( int i = 0; i < dimensions; ++i)
        {
            image_size[i] = limage_size[i];
            image_max[i] = limage_max[i];
            image_min[i] = limage_min[i];

        }

    }
}

Reconstructor::Reconstructor(sdmixer *s,
                             std::vector<sdmixer::Localization> *data,
                             int current_filter)
{
    //if(current_filter == 0)

    this->sdm = s;

    input_data = data;
    min_maxValues = sdm->getPF_min_maxValues();
    minMaxDefined=true;

    input_file = sdm->getCurrentFile();
    QFile qf(input_file);
    QFileInfo fi(qf);

    if(!sdm->getOutputDirectory().isEmpty())
        output_dir = sdm->getOutputDirectory();
    else
        output_dir=fi.path();

    input_base_name = fi.baseName();
    if(input_base_name.contains("_pairs_out"))
        input_base_name.replace("_pairs_out", "");
    if(input_base_name.contains("_filter_out"))
        input_base_name.replace("_filter_out", "");

    tiff_out_file = output_dir.append("/");
    tiff_out_file = tiff_out_file.append(input_base_name);
    if(current_filter != 0)
    {
        tiff_out_file = tiff_out_file.append("_ch");
        tiff_out_file = tiff_out_file.append(QString::number(current_filter));
    }
    tiff_out_file.append(".tif");


    getSettingsFromGUI();

    initData(current_filter);

}


void Reconstructor::outputTIFF()
{
    uint16 spp, bpp, photo, res_unit;

    spp = 1; /* Samples per pixel */
    //3 == rgb, 4 == alpha channel
    bpp = 16; /* Bits per sample */
    photo = PHOTOMETRIC_MINISBLACK;

    qDebug() << "started OutputTIFF";



    boost::iostreams::mapped_file_sink file;
    boost::iostreams::mapped_file_params params;
    if(runConvolution)
        params.path = convolved_image;
    else
        params.path = tiff_temp_file;

    file.open(params);

    if(!file.is_open())
    {
        qDebug() << "no temp file to create OutputTIFF";
        return;
    }
    qDebug() << "array";
    uint16 *array =static_cast<uint16*>((void*)file.data());
    sdm->writeToLogFile("writing TIFF to ", tiff_out_file);
    qDebug() << "tiff out file" << tiff_out_file;

    int sum_pages=0;
    int number_image=0;
    QString current_image;
    int current_page=0;

    int img_sizeZ;
    if(dimensions==3 && !force2D)
        img_sizeZ = image_size[2];
    else
        img_sizeZ = 1;

    out = TIFFOpen(tiff_out_file.toLocal8Bit(), "w");
    if (!out){return;
         qDebug() << "Can't open %s for writing\n";
    }


    for (int page = 0; page < img_sizeZ; page++)
    {
            /*qDebug() << "page: " << page;
            qDebug() << "sum_pages: " << sum_pages;*/

            uint64_t left =  sum_pages * image_size[0] * image_size[1];
            uint64_t right =  std::numeric_limits<uint32_t>::max();
            /*qDebug() << left;
            qDebug() << right;*/
        if( left > 0.8*right) // 80% of the theoritally allowed max file_size
        {
            qDebug() << "current_page = 0 ";
            current_page=0;
            sum_pages=0;
            number_image++;

            TIFFClose(out);

            current_image = tiff_out_file ;
            current_image.replace(".tif", "");
            current_image += (QString::number(number_image) + ".tif");

            //qDebug()<<"TIFF " << number_image << " ready!";

            out = TIFFOpen(current_image.toLocal8Bit(), "w");
            if (!out){return;
                 // "Can't open %s for writing\n"
            }

        }
        if(img_sizeZ!=1)
        {
            if(page%(img_sizeZ/4) == 0)
                qDebug() << "writing page (total): " << page
                         << "  current: " << current_page;
        }

        TIFFSetField(out, TIFFTAG_IMAGEWIDTH, image_size[0] / spp);
        TIFFSetField(out, TIFFTAG_IMAGELENGTH, image_size[1]);
        TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, bpp);
        TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, spp);
        TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(out, TIFFTAG_PHOTOMETRIC, photo);
        TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);
        /* It is good to set resolutions too (but it is not nesessary) */
        xres = yres = 1;//sdm->getPixelSize()*1e-7;
        res_unit = RESUNIT_CENTIMETER;
        TIFFSetField(out, TIFFTAG_XRESOLUTION, xres);
        TIFFSetField(out, TIFFTAG_YRESOLUTION, yres);
        TIFFSetField(out, TIFFTAG_RESOLUTIONUNIT, res_unit);
        TIFFSetField(out, TIFFTAG_COMPRESSION, COMPRESSION_LZW);

         /* We are writing single page of the multipage file */
        TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);

        /* Set the page number */
        TIFFSetField(out, TIFFTAG_PAGENUMBER, current_page, img_sizeZ);

        for (int j = 0; j < image_size[1]; j++)
            TIFFWriteScanline(out, &array[image_size[0]*(j+image_size[1]*page)], j, 0);

                //TIFFWriteScanline(out, &array[j * image_width + image_width*image_height*page], j, 0);
            //


        TIFFWriteDirectory(out);
        sum_pages++;
        current_page++;
    }
    TIFFClose(out);
    file.close();
    qDebug()<<"TIFF out ready!!";

}
