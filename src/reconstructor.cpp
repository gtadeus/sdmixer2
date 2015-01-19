//#include "sdmixer.h"
#include "reconstructor.h"
#include "kdtree.h"




bool Reconstructor::file_exists (const char *filename)
{
    if (FILE *file = fopen(filename, "r"))
    {
        fclose(file);
        return true;
    }
    else
        return false;

}

void Reconstructor::Convolution(Kernel krn)
{
    qDebug() << "started convolution";
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
    params_in.path = tiff_temp_file.toLocal8Bit().data();
    mmap_image_in.open(params_in);
    uint16 *image_in =static_cast<uint16*>((void*)mmap_image_in.data());

    boost::iostreams::mapped_file conv_out_file;
    boost::iostreams::mapped_file_params params_out;
    if(file_exists(convolved_image.toLocal8Bit().data()))
        remove(convolved_image.toLocal8Bit().data());
    params_out.path = convolved_image.toLocal8Bit().data();
    params_out.mode = (std::ios::out | std::ios::in);
    params_out.new_file_size =  NPixel*sizeof(uint16);
    conv_out_file.open(params_out);
    uint16 *image_out = static_cast<uint16*>((void*)conv_out_file.data());


    qDebug() << "kernel.size " << krn.size;
    qDebug() << "kernel.make3D " << krn.make3D;
    qDebug() << "kernel.data[0] " <<krn.data[0];

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
    qDebug() << "convolution main loop ready";
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
        params.path = convolved_image.toLocal8Bit().data();
    else
        params.path = tiff_temp_file.toLocal8Bit().data();

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
void Reconstructor::findNN(QString outputFile)
{
    qDebug() << "starting Nearest-Neighbor Statistics";
    sdm->writeToLogFile("calculating Nearest-Neighbor Statistics");
    qDebug() << "dimensions:" << dimensions;
    QFile NNout(outputFile);
    NNout.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&NNout);
    //KDTree *kd = new KDTree(xyz_not_rounded, dimensions);


    /*std::vector<dblCoordinates>::iterator it_NN;
    for( it_NN = xyz_not_rounded.begin(); it_NN != xyz_not_rounded.end(); ++it_NN )
    {
        dblCoordinates coord = *it_NN;
        Point<3> p;
        for(int i = 0; i < dimensions; ++i)
            p.x[i] = coord.get(i);
        pts.push_back(p);

    }*/
    //NNout.close();
    std::vector<double> results;

    if(dimensions == 2 || force2D)
    {
        KDtree<2> kd(pts2D);
        for (int i = 0; i < pts2D.size(); ++i)
        {
            Int nn[1];
            Doub dn[1];
            kd.nnearest(i, nn, dn, 1);
            //qDebug() << "nn: " << nn[0] ;
            //qDebug() << "dn: " << dn[0] ;
            results.push_back((double)dn[0]);

            out << (double)dn[0] << "\n";
        }
    }
    else
    {
        KDtree<3> kd(pts3D);
        for (int i = 0; i < pts3D.size(); ++i)
        {
            Int nn[1];
            Doub dn[1];
            kd.nnearest(i, nn, dn, 1);
            //qDebug() << "nn: " << nn[0] ;
            //qDebug() << "dn: " << dn[0] ;
            results.push_back((double)dn[0]);

            out << (double)dn[0] << "\n";
        }
    }

    std::sort(results.begin(), results.end());
    qDebug() << "min: " << results[0];
    int size = results.size();
    double median;
    if (size  % 2 == 0)
      {
          median = (results[size / 2 - 1] + results[size / 2]) / 2;
      }
      else
      {
          median = results[size / 2];
      }
    qDebug() << "median: " << median;
    qDebug() << "max: " << results[size-1];

    std::stringstream ss;
    ss << "min = " << results[0] << " (nm), median = " << median << " (nm), max = " << results[size-1] << " (nm)";

    sdm->writeToLogFile(QString::fromStdString(ss.str()));
    //sdm->writeToLogFile("finished Nearest-Neighbor Statistics");
    qDebug() << "finished Nearest-Neighbor Statistics";
}

void Reconstructor::doWork()
{
    qDebug() <<" started Reconstructor worker";

    if(doWorkLater)
        doWorkNow();

    bool repeat =true;

    qDebug() << "doWork max_filter: " << max_filter;

    while(repeat)
    {
        input_file = sdm->getCurrentFile();
        QFile qf(input_file);
        QFileInfo fi(qf);

        if(!sdm->getOutputDirectory().isEmpty())
            output_dir = sdm->getOutputDirectory();
        else
            output_dir=fi.path();

        input_base_name = fi.completeBaseName();
        if(input_base_name.contains("_pairs_out"))
            input_base_name.replace("_pairs_out", "");
        if(input_base_name.contains("_filter_out"))
            input_base_name.replace("_filter_out", "");
        if(input_base_name.contains("grouped_out"))
            input_base_name.replace("grouped_out", "");

        tiff_out_file = output_dir.append("/");
        tiff_out_file = tiff_out_file.append(input_base_name);

        NN_output_file = tiff_out_file;
        NN_output_file.append("_NN_Statistics");

        if(INPUT_FILE == sdmixer::XYZ_FILE && !sdm->getRunFilter())
        {

        }
        else
        {
            NN_output_file = NN_output_file.append("_ch");
            NN_output_file = NN_output_file.append(QString::number(curr_filter));

            tiff_out_file = tiff_out_file.append("_ch");
            tiff_out_file = tiff_out_file.append(QString::number(curr_filter));
        }

        tiff_out_file.append(".tif");
        NN_output_file.append(".txt");

        qDebug() << "NN_output_file: " << NN_output_file;
        qDebug() << "tiff_out_file: " << tiff_out_file;
        qDebug() << "curr_filter: " << curr_filter;

        if (performNNStatistic)
            findNN(NN_output_file);

        repeat=false;
        setArray();

        sdm->writeToLogFile("initialized image data");
        if(runConvolution)
        {
            sdm->writeToLogFile("starting convolution");
            qDebug() << "using Kernel " << curr_filter-1;
            Kernel krn;
            if(!oneConvolutionKernel)
                krn = all_kernels[curr_filter-1];
            else
                krn = all_kernels[0];

            std::stringstream ss;
            ss << "Kernel: sigma_xy = " << krn.sigma_xy << " nm, sigma_z = " << krn.sigma_z << " nm";
            sdm->writeToLogFile(QString::fromStdString(ss.str()));
            Convolution(krn);
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
        if(INPUT_FILE == sdmixer::FILTER_FILE)
        {
            if(curr_filter < max_filter)
            {
                repeat=true;
                curr_filter++;
                initData(curr_filter);
            }
        }
    }

    if(file_exists(tiff_temp_file.toLocal8Bit().data()))
        remove(tiff_temp_file.toLocal8Bit().data());
    /*if(file_exists(fft_src_file))
        remove(fft_src_file);
    if(file_exists(fft_krnl_file))
        remove(fft_krnl_file);*/
    if(file_exists(convolved_image.toLocal8Bit().data()))
        remove(convolved_image.toLocal8Bit().data());

    emit finished();
}


void Reconstructor::setArray()
{
    qDebug() << "start set array";
    maxPixels=1;

    ReslizeZ=sdm->getResliceZ();
    if(ReslizeZ)
    {
        startResliceZ = sdm->getStartRescliceZ();
        endResliceZ = sdm->getEndRescliceZ();
        std::stringstream ss;
        ss << "Reslice Z: start = " << startResliceZ << " , end = " << endResliceZ;
        sdm->writeToLogFile(QString::fromStdString(ss.str()));
        //image_min[2] += startResliceZ;
        //image_max[2] = image_max[2] - (image_size[2]-endResliceZ);
        //image_size[2] = (image_max[2] - image_min[2]) + 1;
        image_size[2] = (endResliceZ - startResliceZ) + 1;
    }

    std::stringstream ss;
    ss << "image_size (px) : ";
    for(int i = 0; i < dimensions; ++i)
    {
        maxPixels*=image_size[i];

        ss << image_size[i] << "  ";
        qDebug()<<"set Array image_size[i]: " << image_size[i];
    }
    sdm->writeToLogFile(QString::fromStdString(ss.str()));


    std::vector<Coordinates>::iterator it;

    int counter=0;
    for( it = xyz.begin(); it != xyz.end(); ++it )
    {
        //qDebug() << "hier";
        Coordinates c = *it;
        for(int i = 0; i < dimensions; ++i)
        {
            it->set(i, c.get(i)-image_min[i]);
            if(i==2 && ReslizeZ)
            {
                it->set(i, it->get(i)-startResliceZ);
                int val = it->get(i);
                if(val < 0 || val >= (endResliceZ-startResliceZ))
                {
                    counter++;
                    xyz.erase(it);
                    --it;
                }

            }
            else
            {
                if(it->get(i)>=image_size[i])
                {
                    //qDebug() << "deleted";
                    counter++;
                    xyz.erase(it);
                    --it;
                }

            }
        }

    }
    qDebug() << "deleted " << counter << " coordinates";
    qDebug() << "maxPixels: " << maxPixels;
    boost::iostreams::mapped_file file;
    boost::iostreams::mapped_file_params params;

    if(file_exists(tiff_temp_file.toLocal8Bit().data()))
        remove(tiff_temp_file.toLocal8Bit().data());
    params.path = tiff_temp_file.toLocal8Bit().data();

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


void Reconstructor::CreateGaussianKernel(Kernel &krn)
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
    else    // kernel[2k+1][2k+1]
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

uint64_t Reconstructor::linearIndex(Coordinates c)
{


    uint64_t val = (uint64_t) c.get(0) +
            (uint64_t) (image_size[0])*( c.get(1) + (uint64_t) (image_size[1])*c.get(2));
    //qDebug() << "x: " << c.get(0) << " y: " << c.get(1) << " z: " << c.get(2) << " ind: " << val;
    return val;
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
      if(max_val[i] != 0)
      {
          image_max[i] = round(max_val[i]/binning[i]);
          image_min[i] = round(min_val[i]/binning[i]);
          image_size[i] = (image_max[i]-image_min[i]) + 1;
          qDebug() << "Reconstructor setMinMax min: " << image_min[i] << "  " << image_max[i] << " size: " <<image_size[i];
          minMaxDefined = true;
      }
    }

    qDebug() << "Reconstructor::setMinMax(sdmixer::min_max m) ready";

}


void Reconstructor::getSettingsFromGUI()
{

    binning[0] = sdm->getReconstructor_xyBinning();
    binning[1] = sdm->getReconstructor_xyBinning();
    binning[2] = sdm->getReconstructor_zBinning();

    runConvolution=sdm->getRunConvolution();

    NM_PER_PX = sdm->getPixelSize();
    force2D = sdm->getForce2D();
    if(force2D)
    {
        dimensions = 2;
        image_size[2]=0;
    }

    perform_hist_eq = sdm->getNonLinearHistEq();
    sqrtCum = sdm->getSqrtCummulation();
    hist_correct_value = sdm->getHisteqCoefficient();
    hist_threshold = sdm->getThreshold();
    oneConvolutionKernel = sdm->getOneKernelForAllChannels();

    if(oneConvolutionKernel)
    {
        Kernel krn;
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
        CreateGaussianKernel(krn);
        all_kernels.push_back(krn);

        qDebug() << "One kernel, krn.sigma_xy " << krn.sigma_xy ;
    }
    else
    {
        vecKernel = sdm->getConvolutionKernel();
        for(auto i : vecKernel)
        {
            Kernel krn;
            sdmixer::gaussian_kernel gk = i;
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
            CreateGaussianKernel(krn);
            all_kernels.push_back(krn);
            qDebug() << "krn.sigma_xy " << krn.sigma_xy ;
        }
    }

    performNNStatistic = sdm->getPerformNNStatistic();

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

    sdm->getHeader(line, columns, min_maxValues, INPUT_FILE);
    if(!force2D)
        dimensions = columns.dimensions;
    setMinMax(min_maxValues);

    if( INPUT_FILE != sdmixer::XYZ_FILE )
    {
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
            }

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

            if(INPUT_FILE == sdmixer::FILTER_FILE)
            {
                loc.filter = dbl_vec[columns.filter];
                if(loc.filter>max_filter)
                    max_filter=loc.filter;
            }

            sdm->pushBackLocalization(loc);
        }
        input_data = sdm->getPfOutput();
        qDebug() << "loaded ";

        //qDebug() << " is this a filter_file?: " << FilterInput;
        if(INPUT_FILE == sdmixer::FILTER_FILE)
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
            strLine.erase(strLine.begin(),
                          std::find_if(strLine.begin(), strLine.end(),
                                       std::bind1st(std::not_equal_to<char>(), ' ')));
            boost::split(v, strLine, boost::is_any_of("\t "));

            int index = 0;
            for (auto i : v)
            {
                index++;
                dbl_vec.push_back(strtod(i.c_str(), NULL));
            }
            //qDebug() << "v.size() " << v.size();


            Coordinates c;
            Point<2> point2D;
            Point<3> point3D;
            for(int i = 0; i < dimensions; ++i)
            {
                int ind = columns.getXYZCol(i);
                c.set(i, round(dbl_vec[ind]/binning[i]));

                if(performNNStatistic)
                {
                    if(dimensions == 2 || force2D)
                    {
                        for(int i = 0; i < dimensions; ++i)
                            point2D.x[i] = c.get(i);
                    }
                    else
                    {
                        for(int i = 0; i < dimensions; ++i)
                            point3D.x[i] =  c.get(i);
                    }
                }
            }
            xyz.push_back(c);

            if(performNNStatistic)
            {
                if(dimensions == 2 || force2D)
                    pts2D.push_back(point2D);
                else
                    pts3D.push_back(point3D);
            }
        }
        // get File Name
        input_file = sdm->getCurrentFile();

        QFile qf(input_file);
        QFileInfo fi(qf);

        if(!sdm->getOutputDirectory().isEmpty())
            output_dir = sdm->getOutputDirectory();
        else
            output_dir=fi.path();

        input_base_name = fi.completeBaseName();

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
    tiff_temp_file = sdm->getTiffTempFile();
    convolved_image = sdm->getConvImgTempFile();


    qDebug() << "tiff_temp_file: " << tiff_temp_file;
    qDebug() << "convolved_image: " << convolved_image;
    /*QFile f(xyz_file);
    f.open(QIODevice::ReadOnly| QIODevice::Text);

    QTextStream in(&f);
    QString line = in.readLine();

    sdm->getHeader(line, columns, min_maxValues, INPUT_FILE);
    dimensions = columns.dimensions;
    f.close();*/
    getSettingsFromGUI();

    doWorkLater=true;

}
void Reconstructor::initData(int current_filter)
{
    qDebug() << "Rconstructor init Data";

    setMinMax(min_maxValues);

    int limage_min[3]={0};
    int limage_max[3]={0};
    int limage_size[3]={0};

    int index = 0;
    std::vector<sdmixer::Localization>::iterator it;

    for( it = input_data->begin(); it != input_data->end(); ++it )
    {
        sdmixer::Localization loc = *it;
        Coordinates c;
        Point<2> point2D;
        Point<3> point3D;
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
                //for NN Analysis
                if(performNNStatistic)
                {
                    if(dimensions == 2 || force2D)
                    {
                        for(int i = 0; i < dimensions; ++i)
                            point2D.x[i] = loc.getShortDim(i);
                    }
                    else
                    {
                        for(int i = 0; i < dimensions; ++i)
                            point3D.x[i] = loc.getShortDim(i);
                    }
                }

                if(limage_min[i]>short_val)
                        limage_min[i]=short_val;
                if(limage_max[i]<short_val)
                    limage_max[i]=short_val;

                ++index;
            }
            xyz.push_back(c);
            if(performNNStatistic)
            {
                if(dimensions == 2 || force2D)
                    pts2D.push_back(point2D);
                else
                    pts3D.push_back(point3D);
            }
            for(int i = 0; i < dimensions; ++i)
            {
                int long_val = round(loc.getLongDim(i)/binning[i]);
                if(limage_min[i]>long_val)
                        limage_min[i]=long_val;
                if(limage_max[i]<long_val)
                    limage_max[i]=long_val;

                c.set(i, long_val);
                //for NN Analysis
                if(performNNStatistic)
                {
                    if(dimensions == 2 || force2D)
                    {
                        for(int i = 0; i < dimensions; ++i)
                            point2D.x[i] = loc.getLongDim(i);
                    }
                    else
                    {
                        for(int i = 0; i < dimensions; ++i)
                            point3D.x[i] = loc.getLongDim(i);
                    }
                }
            }
            xyz.push_back(c);

            if(performNNStatistic)
            {
                if(dimensions == 2 || force2D)
                    pts2D.push_back(point2D);
                else
                    pts3D.push_back(point3D);
            }
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
    curr_filter = current_filter;
    qDebug() << "Reconstruction Constructor, curr_filter: " << curr_filter;

    this->sdm = s;
    tiff_temp_file = sdm->getTiffTempFile();
    convolved_image = sdm->getConvImgTempFile();

    //tiff_temp_file = tf.toLocal8Bit().data();
    //convolved_image = cv.toLocal8Bit().data();
    qDebug() << "tiff_temp_file: " << tiff_temp_file;
    qDebug() << "convolved_image: " << convolved_image;

    input_data = data;
    input_file = sdm->getCurrentFile();

    QFile f(input_file);
    f.open(QIODevice::ReadOnly| QIODevice::Text);

    QTextStream in(&f);
    QString line = in.readLine();

    sdm->getHeader(line, columns, min_maxValues, INPUT_FILE);
    dimensions = columns.dimensions;
    f.close();

    QFile qf(input_file);
    QFileInfo fi(qf);

    if(!sdm->getOutputDirectory().isEmpty())
        output_dir = sdm->getOutputDirectory();
    else
        output_dir=fi.path();

    input_base_name = fi.completeBaseName();
    if(input_base_name.contains("_pairs_out"))
        input_base_name.replace("_pairs_out", "");
    if(input_base_name.contains("_filter_out"))
        input_base_name.replace("_filter_out", "");
    if(input_base_name.contains("grouped_out"))
        input_base_name.replace("grouped_out", "");

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
        params.path = convolved_image.toLocal8Bit().data();
    else
        params.path = tiff_temp_file.toLocal8Bit().data();

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
