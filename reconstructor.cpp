#include "reconstructor.h"


void Reconstructor::setArray()
{
    qDebug() << "start set array";

    for(int i = 0; i < dimensions; ++i)
    {
        maxPixels*=(image_max[i]+1);
    }

    boost::iostreams::mapped_file file;
    boost::iostreams::mapped_file_params params;
    remove(tmpfile) ;
    params.path = "tmp.file";

    params.mode = (std::ios::out | std::ios::in);
    params.new_file_size =  maxPixels*sizeof(uint8);

    file.open(params);

    uint8 * array =static_cast<uint8*>((void*)file.data());

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

   /* array = new uint8[maxPixels];
    qDebug() << "initializing... map_pixel = " << maxPixels;
    for ( int i = 0; i < maxPixels; ++i)
    {
        array[i]=0;
    }
    qDebug() << "populating...";
    for( auto i : xyz)
    {
        //qDebug() << linearIndex(i) << " : " << i.x << " " << i.y << " " << i.z;
        array[linearIndex(i)] += 1;
    }

    qDebug() << "finish array";
*/
    //delete[] array;

    setKernel();
    CreateGaussianKernel();
    /*for(int i = 0; i < krn.size * krn.size; ++i)
        qDebug() << krn.data[i];*/

    Convolution();
}
void Reconstructor::zeroPad()
{

}
int Reconstructor::pow2roundup (int x)
{
    if (x < 0)
        return 0;
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return x+1;
}
void Reconstructor::Convolution()
{
    boost::iostreams::mapped_file pad_img, pad_filter, fft_img, fft_filter, pad_img, out_img;
    boost::iostreams::mapped_file_params params;
    remove(tmpfile);
    params.path = "tmp.file";

    params.mode = (std::ios::out | std::ios::in);
    params.new_file_size =  maxPixels*sizeof(uint8);

    file.open(params);

    uint8 * array =static_cast<uint8*>((void*)file.data());
    int img_sizeX = 5;//image_max[0] + 1;
    int img_sizeY = 5;//image_max[1] + 1;
    int img_sizeZ = 1;//image_max[2] + 1;
    int NPixel=img_sizeX*img_sizeY*img_sizeZ;


    double *image = new double[img_sizeX*img_sizeY*img_sizeZ];
    double *image_out = new double[img_sizeX*img_sizeY*img_sizeZ];

    for (int i = 0; i < NPixel; ++i)
        image[i]=0;

    image[12]=1;
    int M=img_sizeX;
    int N=img_sizeY;
    int Z=img_sizeZ;


    fftw_plan p_forw_src, p_forw_kernel, pinv;

    double scale = 1.0 / (M * N * Z);

    int h_fftw = pow2roundup(M + krn.k);
    int w_fftw = pow2roundup(N + krn.k);

    int z_fftw;
    if(krn.make3D)
        z_fftw = pow2roundup(Z + krn.k);
    else
        z_fftw = 1;

    int NPixelPadded = h_fftw * w_fftw * z_fftw;
    int NPixelFFT = h_fftw * w_fftw * (z_fftw/2+1);

    qDebug() << h_fftw << "  " << h_fftw << "  " << z_fftw;

    double *in_src = new double[NPixelPadded];
    fftw_complex *out_src = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NPixelFFT);
    double *in_kernel = new double[NPixelPadded];
    fftw_complex *out_kernel = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NPixelFFT);
    //fftw_complex *mult = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NPixelFFT);
    double *dst_fft = new double[NPixelPadded];
    double *dst = new double[NPixel];

    for (int i = 0; i < NPixelPadded; ++i)
    {
            in_src[i] = 0;
            in_kernel[i] = 0;
    }

    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j)
            for ( int k = 0; k < Z; ++k)
            {
                in_src[i+w_fftw*j+w_fftw*h_fftw*k] = image[i+img_sizeX*j+img_sizeX*img_sizeY*k];
            }
    for (int i = 0; i < krn.size; ++i)
        for (int j = 0; j < krn.size; ++j)
            for(int k = 0; k < krn.size_z; ++k)
            {
                in_kernel[i+w_fftw*j + w_fftw*h_fftw*k] = krn.data[i+krn.size*j+krn.size*krn.size*k];
            }

   /* qDebug() << "insrc:";
    for (int i=0;i<NPixelPadded; ++i)
        qDebug()<< in_src[i];
    qDebug() << "kernel:";
    for (int i=0;i<NPixelPadded; ++i)
        qDebug()<< in_kernel[i];
*/
    //CONVN : SAME


    /* create plan for FFT; */
    // image
    //p_forw_src = fftw_plan_dft_r2c_2d(h_fftw, w_fftw, in_src, out_src, FFTW_ESTIMATE);
    //p_forw_kernel = fftw_plan_dft_r2c_2d(h_fftw, w_fftw, in_kernel, out_kernel, FFTW_ESTIMATE);
    int n[3];
    n[0]=h_fftw;
    n[1]=w_fftw;
    n[2]=z_fftw;

    p_forw_src = fftw_plan_dft_r2c(3, n, in_src, out_src, FFTW_ESTIMATE);
    p_forw_kernel = fftw_plan_dft_r2c(3, n, in_kernel, out_kernel, FFTW_ESTIMATE);

    /* create plan for IFFT; */
    //pinv = fftw_plan_dft_c2r_2d(h_fftw, w_fftw, mult, dst_fft, FFTW_ESTIMATE);
    pinv = fftw_plan_dft_c2r(3, n, out_kernel, dst_fft, FFTW_ESTIMATE);
    qDebug()<< "created plan";
    fftw_execute(p_forw_src);
    qDebug()<< "executed plan1";
    fftw_execute(p_forw_kernel);
    qDebug()<< "executed plan2";

    for (int i = 0; i < NPixelFFT; ++i)
    {
        out_kernel[i][0] = (out_src[i][0] * out_kernel[i][0] - out_src[i][1] * out_kernel[i][1]);
        out_kernel[i][1] = (out_src[i][0] * out_kernel[i][1] + out_src[i][1] * out_kernel[i][0]);
    }

    fftw_execute(pinv);

    //scale


   /* for(int i=0; i< NPixelPadded; ++i)
        qDebug() << dst_fft[i];*/

    int index = 0;

    for (int i = floor(krn.size/2); i < (M+floor(krn.size/2)); ++i)
        for (int j = floor(krn.size/2); j < (N+floor(krn.size/2)); ++j)
            for(int k = floor(krn.size_z/2); k < (Z+floor(krn.size_z/2)); ++k)
            {
                dst[index] = dst_fft[i+w_fftw*j+w_fftw*h_fftw*k];
                ++index;
            }

    for (int i=0;i<NPixel; ++i)
            dst[i]*=scale;

   /*for(int i=0; i< NPixel; ++i)
        qDebug() << dst[i];*/

    qDebug() << "convolution ready!";

}
void Reconstructor::CreateGaussianKernel()
{
    float f;
    int index = 0;
    float sum = 0;
    float e = 2.71828;

    // kernel[2k+1][2k+1][2k+1]
    ///////////////////////////////
    //e ^ -(x*x + y*y)/2 * sigma * sigma
    ////////////////////////////////

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
    krn.make3D=false ;
    krn.data = new float[krn.size*krn.size*krn.size];
}


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

uint64_t Reconstructor::linearIndex(Coordinates c)
{
    int index = 0;
    return c.get(0) + (uint64_t) (image_max[0]+1)*c.get(1) + (uint64_t) (image_max[1]+1)*(image_max[0]+1)*c.get(2);
    //return image_max[0]*image_max[1]*c.get(0)+image_max[1]*c.get(1)+c.get(2);
    for ( int i = 0; i < dimensions; ++i)
    {
        int prod=1;
        for(int j = (i+1); j < dimensions; ++j)
        {
            prod *= (image_max[j]-1);
        }
        index += (c.get(i) * prod);
    }

    return index;
}



void Reconstructor::setMinMax(double min_x, double max_x,
                              double min_y, double max_y,
                                double min_z, double max_z)
{
    // min max, convert from m to nm

    min_val[0]=min_x * 1e9;
    max_val[0]=max_x * 1e9;

    min_val[1]=min_y * 1e9;
    max_val[1]=max_y * 1e9;

    min_val[2]=min_z * 1e9;
    max_val[2]=max_z * 1e9;
}

void Reconstructor::XYZfromFilter(std::vector<PairFinder::Localization>& data)
{
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
    }
    qDebug() << "ready " << xyz[0].x << "  " << xyz[0].y<< "  "  << xyz[0].z;
}

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

Reconstructor::Reconstructor(sdmixer *s)
{
    this->sdm = s;
    //load filter output
    //convert xyz coordinates to
    //round vectors to full nm
    // get min max nm

    //horzcat input

    // CONVERT DATA TO XYZ
    // dann kann für jeden kanal einzelm geblurrt werden


    // image pixel

    // extract each

    // blur for each channel

}


void Reconstructor::run()
{
    // extract cols

    //output TIFF



    outputTIFF("out.tif");

}
void Reconstructor::convertTo2D()
{

}
void Reconstructor::createKernel()
{

}
void Reconstructor::convolve(double *kernel)
{

}
void Reconstructor::hist_correct()
{

}
void Reconstructor::outputTIFF(QString path)
{
    qDebug() << "writing TIFF";
    uint16 spp, bpp, photo, res_unit;

    spp = 1; /* Samples per pixel */
    //3 == rgb, 4 == alpha channel
    bpp = 8; /* Bits per sample */
    photo = PHOTOMETRIC_MINISBLACK;
    uint64_t image_width = image_max[0] + 1;
    uint64_t image_height = image_max[1] + 1;
    int image_z = image_max[2] + 1;

    boost::iostreams::mapped_file_sink file;
    boost::iostreams::mapped_file_params params;
    params.path = "tmp.file";
    //params.mode = std::ios::in;
    file.open(params);

    uint8 *array =static_cast<uint8*>((void*)file.data());

    uint64_t sum_pages=0;
    int number_image=0;
    QString current_image;
    int current_page=0;

    out = TIFFOpen(path.toLocal8Bit(), "w");
    if (!out){return;
         // "Can't open %s for writing\n"
    }

    for (int page = 0; page < image_z; page++)
    {
            qDebug() << "page: " << page;
            qDebug() << "sum_pages: " << sum_pages;

            uint64_t left =  sum_pages * image_width * image_height;
            uint64_t right =  std::numeric_limits<uint32_t>::max();
            qDebug() << left;
            qDebug() << right;
        if( left > 0.8*right) // 80% of the theoritally allowed max file_size
        {
            qDebug() << "current_page = 0 ";
            current_page=0;
            sum_pages=0;
            number_image++;
            //TIFFWriteDirectory(out);
            TIFFClose(out);

            current_image = path ;
            current_image.replace(".tif", "");
            current_image += (QString::number(number_image) + ".tif");

            qDebug()<<"TIFF " << number_image << " ready!";

            out = TIFFOpen(current_image.toLocal8Bit(), "w");
            if (!out){return;
                 // "Can't open %s for writing\n"
            }

        }
        if(page%10 == 0)
            qDebug() << "writing page (total): " << page << "  current: " << current_page;
        TIFFSetField(out, TIFFTAG_IMAGEWIDTH, image_width / spp);
        TIFFSetField(out, TIFFTAG_IMAGELENGTH, image_height);
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
        TIFFSetField(out, TIFFTAG_PAGENUMBER, current_page, image_z);

        for (int j = 0; j < image_height; j++)
            TIFFWriteScanline(out, &array[j * image_width + image_width*image_height*page], j, 0);

        TIFFWriteDirectory(out);
        sum_pages++;
        current_page++;
    }
    TIFFClose(out);
    file.close();
    qDebug()<<"TIFF out ready!!";

}
