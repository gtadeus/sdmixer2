#include "reconstructor.h"


void Reconstructor::setArray()
{
    qDebug() << "start set array";

    for(int i = 0; i < dimensions; ++i)
    {
        qDebug()<<"set Array image_max[i]: " <<image_max[i];
        maxPixels*=(image_max[i]+1);
    }
    qDebug()<<maxPixels;
    boost::iostreams::mapped_file file;
    boost::iostreams::mapped_file_params params;

    remove(tiff_temp_file);
    params.path = tiff_temp_file;

    params.mode = (std::ios::out | std::ios::in);
    params.new_file_size =  maxPixels*sizeof(uint8);

    file.open(params);

    uint8 *array =static_cast<uint8*>((void*)file.data());

    //data[0] = 12;
    //array = new uint8[maxPixels];
    if (file.is_open())
    {
    qDebug() << "initializing... max_pixel = " << maxPixels;
       for ( int i = 0; i < maxPixels; ++i)
       {
           array[i]=0;
       }


       //qDebug() << array[maxPixels];
       qDebug() << "populating...";
        uint64_t temp_max =0;
       for(std::vector<int>::size_type i = 0; i != xyz.size(); i++)
       {
           //qDebug() << linearIndex(xyz[i]) << " : " << xyz[i].x << " " << xyz[i].y << " " << xyz[i].z;
           uint64_t lin_index = linearIndex(xyz[i]);
           if (i == 0)
               temp_max = lin_index;
           else
               if(lin_index > temp_max)
                   temp_max = lin_index;
            array[lin_index]+=1;
       }
       qDebug() << temp_max;
        qDebug() << "ready!";

       /*for( auto i : xyz)
       {
           //qDebug() << linearIndex(i) << " : " << i.x << " " << i.y << " " << i.z;
           array+=linearIndex(i);
           *array = 1;
           array=start;
       }*/
    }
    file.close();


    qDebug() << "finish array";

}

int Reconstructor::pow2roundup (int x)
{
    /*if (x < 0)
        return 0;
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return x+1;*/
    return x+1;
}
void Reconstructor::Convolution()
{
    int img_sizeX = image_max[0] + 1;
    int img_sizeY = image_max[1] + 1;
    int img_sizeZ = image_max[2] + 1;
    int NPixel=img_sizeX*img_sizeY*img_sizeZ;

    //double *image_out = new double[img_sizeX*img_sizeY*img_sizeZ];

    //for (int i = 0; i < NPixel; ++i)
        //image[i]=0;

    //image[13]=1;

    int M=img_sizeX;
    int N=img_sizeY;
    int Z=img_sizeZ;


    double scale = 1.0 / (M * N * Z);

    int w_fftw = pow2roundup(M + krn.k);
    int h_fftw = pow2roundup(N + krn.k);

    int z_fftw;
    if(krn.make3D)
        z_fftw = pow2roundup(Z + krn.k);
    else
        z_fftw = 1;


    int NPixelPadded = h_fftw * w_fftw * z_fftw;
    int NPixelFFT = h_fftw * w_fftw * ceil((double)z_fftw/2.0 +1.0);

    qDebug() << w_fftw << "  " << h_fftw << "  " << z_fftw;
    qDebug() << M << "  " << N << "  " << Z;

    boost::iostreams::mapped_file image_in;
    boost::iostreams::mapped_file_params params;
    params.mode = (std::ios::out | std::ios::in);
    params.path = tiff_temp_file;
    image_in.open(params);
    uint8 *img_in =static_cast<uint8*>((void*)image_in.data());

    boost::iostreams::mapped_file fft_src, fft_kernel;
    boost::iostreams::mapped_file_params params_scr, params_kernel;
    remove(fft_src_file);
    remove(fft_krnl_file);
    params_scr.path = fft_src_file;
    params_kernel.path = fft_krnl_file;

    params_scr.mode = (std::ios::out | std::ios::in);
    params_kernel.mode = (std::ios::out | std::ios::in);

    params_scr.new_file_size =  NPixelFFT*sizeof(fftwf_complex);
    params_kernel.new_file_size =  NPixelFFT*sizeof(fftwf_complex);

    fft_src.open(params_scr);
    fft_kernel.open(params_kernel);

    fftwf_plan p_forw_src, p_forw_kernel, pinv;

    //double *out_src = (double*) fftwf_malloc(sizeof(fftwf_complex) *  NPixelFFT);
    //double *out_kernel = (double*) fftwf_malloc(sizeof(fftwf_complex) * NPixelFFT);

    qDebug() << "fft tempfiles created";
    float *out_src = static_cast<float*>((void*)fft_src.data());
    float *out_kernel = static_cast<float*>((void*)fft_kernel.data());

    for (int i = 0; i < 2*NPixelFFT; ++i)
    {

        out_src[i]=0;
        out_kernel[i]=0;
    }
    qDebug() << "conv init ready";

     double tempsum = 0;
     double tempsum2 = 0;

     for (int i=0;i<NPixel; ++i)
         tempsum += (double) img_in[i];

    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j)
            for ( int k = 0; k < Z; ++k)
            {
                uint64_t indTIFF = linearIndexTIFF(i, j, k, img_sizeX, img_sizeY);
                uint64_t indFFT = linearIndex3DFFT(i, j, k, z_fftw, h_fftw);
                //int ind = i+img_sizeX*( j+ img_sizeY*k );
                //int indFFT = k + (z_fftw+2)*(j+h_fftw*i);
                int d = img_in[indTIFF];

                out_src[indFFT] = (float)d;

            }
    for (int i = 0; i < krn.size; ++i)
        for (int j = 0; j < krn.size; ++j)
            for(int k = 0; k < krn.size_z; ++k)
            {
                uint64_t indFFT = linearIndex3DFFT(i, j, k, z_fftw, h_fftw);
                //int indFFT = k + (z_fftw+2)*(j+h_fftw*i);
                uint64_t indTIFF = linearIndexTIFF(i, j, k, krn.size, krn.size);
                out_kernel[indFFT] = krn.data[indTIFF];
            }
    for (int i=0;i<NPixelFFT; ++i)
        tempsum2+=out_src[i];

    qDebug() << "tempsum1: " << tempsum << "tempsum2: " << tempsum2;


    //CONVN : SAME

    /* create plan for FFT; */
    int n[3];
    n[0]=w_fftw;
    n[1]=h_fftw;
    n[2]=z_fftw;

    fftw_init_threads();
    fftw_plan_with_nthreads(8);

    p_forw_src = fftwf_plan_dft_r2c(3, n, out_src, (fftwf_complex*) out_src, FFTW_ESTIMATE);
    p_forw_kernel = fftwf_plan_dft_r2c(3, n, out_kernel, (fftwf_complex*) out_kernel, FFTW_ESTIMATE);

    // iFFT
    pinv = fftwf_plan_dft_c2r(3, n, (fftwf_complex*) out_kernel, out_kernel, FFTW_ESTIMATE);

    qDebug()<< "created plan";
    fftwf_execute(p_forw_src);
    qDebug()<< "executed plan1";
    fftwf_execute(p_forw_kernel);
    qDebug()<< "executed plan2";


    for (int i = 0; i < 2*(NPixelFFT); i+=2)
    {
        float real_s, real_k, img_s, img_k;
        real_s = out_src[i];
        real_k = out_kernel[i];
        img_s = out_src[i+1];
        img_k = out_kernel[i+1];
        //qDebug() << real_s << "  " << img_s << "  "  << real_k << "  " << img_k;
        out_kernel[i] = (real_s * real_k - img_s * img_k);
        out_kernel[i+1] = (real_s * img_k + img_s * real_k);
        //qDebug() << out_kernel[i] << "  " << out_kernel[i+1];
    }
    //for(int i=0; i< NPixelFFT; i++)
         //qDebug() << out_kernel[i];

    fftwf_execute(pinv);
    qDebug() << "executed plan3";

    uint64_t index = 0;

    for (int i = floor(krn.size/2); i < (M+floor(krn.size/2)); ++i)
        for (int j = floor(krn.size/2); j < (N+floor(krn.size/2)); ++j)
            for(int k = floor(krn.size_z/2); k < (Z+floor(krn.size_z/2)); ++k)
            {

                //int indFFT = i + (w_fftw)*(j+h_fftw*k);
                uint64_t indTIFF = linearIndexTIFF(i-floor(krn.size/2), j-floor(krn.size/2), k-floor(krn.size/2), img_sizeX, img_sizeY);
                uint64_t indFFT = linearIndex3DFFT(i, j, k, z_fftw, h_fftw);
                /*if(indTIFF > 510087299)
                    qDebug() << "indTiff: " << indTIFF;
                if(indFFT > 527106047)
                    qDebug() << "indFFT: " << indFFT;*/

                //int indFFT = k + (z_fftw+2)*(j+h_fftw*i);
                //int index = (i-floor(krn.size/2)) + img_sizeX*((j-floor(krn.size/2))+img_sizeY*(k-floor(krn.size_z/2)));
                int val = (int) out_kernel[indFFT] * scale * 255.0;
                img_in[indTIFF] = val;
                int dd = img_in[indTIFF];

                if ( index == 0)
                {
                    dbl_image_min=dd;
                    dbl_image_max=dd;
                }
                if (dbl_image_min>dd)
                    dbl_image_min=dd;
                if(dbl_image_max<dd)
                    dbl_image_max=dd;

                //qDebug() << out_kernel[k+(h_fftw+2)*j+z_fftw*(h_fftw+2)*i];
                ++index;
            }
    qDebug() << "image_max_charval: " << dbl_image_max;
    qDebug() << "image_min_charval: " << dbl_image_min;

    image_in.close();
    fft_src.close();
    fft_kernel.close();

    fftw_cleanup_threads();

    qDebug() << "convolution ready!";

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

/*
void Reconstructor::getMinMax()
{
    int counter = 0;
    for ( auto coord : xyz)
    {
        for (int i = 0; i < dimensions; ++i)
        {
            if(counter == 0)
            {
                image_min[i] = coord.get(i);
                image_max[i] = coord.get(i);
            }
            if (coord.get(i) > image_max[i])
                image_max[i] = coord.get(i);
            if (coord.get(i) < image_min[i])
                image_min[i] = coord.get(i);
        }
        ++counter;
    }
    for (int i = 0; i < dimensions; ++i)
    {
        qDebug() << image_max[i];
    }
}
*/
uint64_t Reconstructor::linearIndex(Coordinates c)
{

    return (uint64_t) c.get(0) + (uint64_t) (image_max[0]+1)*( c.get(1) + (uint64_t) (image_max[1]+1)*c.get(2));
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


void Reconstructor::setMinMax(double min_x, double max_x,
                              double min_y, double max_y,
                              double min_z, double max_z)
{
    // min max, convert from m to nm
    double min_val[3]={0};
    double max_val[3]={0};

    min_val[0]=min_x * 1e9;
    max_val[0]=max_x * 1e9;

    min_val[1]=min_y * 1e9;
    max_val[1]=max_y * 1e9;

    min_val[2]=min_z * 1e9;
    max_val[2]=max_z * 1e9;

    for (int i = 0; i < dimensions; ++i)
    {
      if(round(max_val[i]) != 0)
      {
          image_max[i] = round(max_val[i]/binning[i]);
          image_size[i] = image_max[i]+1;
      }
    }
    //image_max[i] = round(max_val[i])/xy_binning;
    //image_max[2] = round(max_val[2])/z_binning;
}

void Reconstructor::XYZfromFilter()
{/*
    for (auto i : data)
    {

            bool ShortOK = false;
            bool LongOK = false;
            for ( int d = 0; d < dimensions; ++d)
            {
                if ( i.getShortDim(d) >= min_val[d])
                    if(i.getShortDim(d) <= max_val[d] && max_val[d] != 0)
                        ShortOK=true;
                if ( i.getLongDim(d) >= min_val[d])
                    if( i.getLongDim(d) <= max_val[d] && max_val[d] != 0)
                        LongOK=true;


            }
            if (ShortOK && LongOK)
            {
                Coordinates c;
                c.x = round(i.xShort/xy_binning);
                c.y = round(i.yShort/xy_binning);
                c.z = round(i.zShort/z_binning);
                xyz.push_back(c);
                c.x = round(i.xLong/xy_binning);
                c.y = round(i.yLong/xy_binning);
                c.z = round(i.zLong/z_binning);
                xyz.push_back(c);
            }
    }*/
    //qDebug() << "ready, first (x,y,z) tupel:" << xyz[0].x << "  " << xyz[0].y<< "  "  << xyz[0].z;
}
/*
void Reconstructor::getIndexFromXYZ()
{
    // N1xN2xN3 Array
    int Nl[3]={round(max_x), round(max_y), round(max_z)};
    int index = 0;
    std::vector<int> image;
    for ( auto i : xyz)
    {
        for (int k = 0; k < dimensions; ++k)
        {
            for(int l = k; k < dimensions; ++l)
            {
                index*=Nl[k];
            }
            if(k==0)
                index*=i.x;
            if(k==1)
                index*=i.y;
            if(k==2)
                index*=i.z;

            image.push_back(int(index));

        }
    }
}
*/
static bool pred( const std::string &s ) {
  // ...
}
Reconstructor::Reconstructor(sdmixer *s, std::vector<PairFinder::Localization>& data)
{
    this->sdm = s;
    binning[0] = s->getReconstructor_xyBinning();
    binning[1] = s->getReconstructor_xyBinning();
    binning[2] = s->getReconstructor_zBinning();
    int index = 0;

    for (auto loc : data)
    {
        Coordinates c;

        for(int i = 0; i < dimensions; ++i)
        {
            int short_val = round(loc.getShortDim(i)/binning[i]);
            if (index == 0)
            {
                image_min[i]=short_val;
                image_max[i]=short_val;
            }
            c.set(i, short_val);
            if(image_min[i]>short_val)
                    image_min[i]=short_val;
            if(image_max[i]<short_val)
                image_max[i]=short_val;

            ++index;
        }
        xyz.push_back(c);
        for(int i = 0; i < dimensions; ++i)
        {
            int long_val = round(loc.getLongDim(i)/binning[i]);
            if(image_min[i]>long_val)
                    image_min[i]=long_val;
            if(image_max[i]<long_val)
                image_max[i]=long_val;

            c.set(i, long_val);
        }
        xyz.push_back(c);

    }
    for ( int i = 0; i < dimensions; ++i)
    {
        image_size[i]=image_max[i]+1;
        qDebug() << "Reconstructor (image size, from input data): " << image_size[i];
    }


}


void Reconstructor::outputTIFF()
{
    uint16 spp, bpp, photo, res_unit;

    spp = 1; /* Samples per pixel */
    //3 == rgb, 4 == alpha channel
    bpp = 8; /* Bits per sample */
    photo = PHOTOMETRIC_MINISBLACK;

    qDebug() << "started OutputTIFF";

    boost::iostreams::mapped_file_sink file;
    boost::iostreams::mapped_file_params params;
    params.path = tiff_temp_file;

    file.open(params);

    if(!file.is_open())
    {
        qDebug() << "no temp file to create OutputTIFF";
        return;
    }
    uint8 *array =static_cast<uint8*>((void*)file.data());

    int sum_pages=0;
    int number_image=0;
    QString current_image;
    int current_page=0;

    out = TIFFOpen(tiff_out_file.toLocal8Bit(), "w");
    if (!out){return;
         // "Can't open %s for writing\n"
    }

    for (int page = 0; page < image_size[2]; page++)
    {
            /*qDebug() << "page: " << page;
            qDebug() << "sum_pages: " << sum_pages;*/

            uint64_t left =  sum_pages * image_size[0] * image_size[1];
            uint64_t right =  std::numeric_limits<uint32_t>::max();
            /*qDebug() << left;
            qDebug() << right;*/
        if( left > 0.8*right) // 80% of the theoritally allowed max file_size
        {
            //qDebug() << "current_page = 0 ";
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
        if(page%(image_size[2]/4) == 0)
            qDebug() << "writing page (total): " << page << "  current: " << current_page;

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
        TIFFSetField(out, TIFFTAG_PAGENUMBER, current_page, image_size[2]);

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
