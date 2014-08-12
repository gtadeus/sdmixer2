#include "reconstructor.h"
#include "sdmixer.h"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

void Reconstructor::Convolution2()
{
    qDebug() << " started convolution 2";
    uint64_t img_sizeX = image_size[0];
    uint64_t img_sizeY = image_size[1];
    uint64_t img_sizeZ = image_size[2];

    uint64_t M=img_sizeX;
    uint64_t N=img_sizeY;
    uint64_t Z=img_sizeZ;

    uint64_t NPixel=img_sizeX*img_sizeY*img_sizeZ;
    //double *image = new double[NPixel];

   /* for (int i = 0; i < NPixel; ++i)
    {
        image[i] = 0;
    }
    //image[4]=1;
    image[13] = 1;*/

    //CreateGaussianKernel();
    //qDebug() << krn.sigma_xy;
    //qDebug() << krn.size;

    boost::iostreams::mapped_file mmap_image_in;
    boost::iostreams::mapped_file_params params_in;
    params_in.mode = (std::ios::out | std::ios::in);
    params_in.path = tiff_temp_file;
    mmap_image_in.open(params_in);
    uint8 *image_in =static_cast<uint8*>((void*)mmap_image_in.data());

    boost::iostreams::mapped_file conv_out_file;
    boost::iostreams::mapped_file_params params_out;
    remove(convolved_image);
    params_out.path = convolved_image;
    params_out.mode = (std::ios::out | std::ios::in);
    params_out.new_file_size =  NPixel*sizeof(uint8);
    conv_out_file.open(params_out);
    uint8 *image_out = static_cast<uint8*>((void*)conv_out_file.data());



   /* double scale = 1.0 / (M * N * Z);


    int w_fft = (M + 2*krn.k);
    int h_fft = (N + 2*krn.k);

    int z_fft = 1;//(Z + 2*krn.k);

    int NPixelPadded = w_fft*h_fft*z_fft;*/

    //double *padded_in = new double[NPixelPadded];
    //double *padded_out = new double[NPixelPadded];
    //double *image_out = new double[NPixel];

    /*for (int i = 0; i < NPixelPadded; ++i)
    {
        padded_in[i] = 0;
    }*/

/*
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j)
            //for ( int k = 0; k < Z; ++k)
            {
            int k = 0;
                uint64_t indTIFF = img_sizeY * img_sizeX * k + img_sizeX * j + i;
                uint64_t indFFT = h_fft * w_fft * k + w_fft * j + i;

                int d = image[indTIFF];

                padded_in[indFFT] = (float)d;

            }*/

    qDebug() << "kernel_:size: " << krn.size;
    qDebug() << "kernel_3D: " << krn.make3D;

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
                        //
                        if(rowIndex >=0 && rowIndex < N &&
                                colIndex >= 0 && colIndex < M &&
                                zIndex >= 0 && zIndex < Z)
                        {
                            //sum+=image[M*rowIndex+colIndex]*
                            //        krn.data[krn.size*mm + nn];
                            sum+=(double)image_in[M*N*zIndex+M*rowIndex+colIndex]*
                                    krn.data[krn.size*krn.size*pp+krn.size*mm + nn];
                            //sum+=padded_in[M * (j + n) + (i + m)] *
                            //krn.data[krn.size*(krn.size-1-n) + (krn.size-1-m)];
                            //sum+=padded_in[M * (j + m) + (i + n)] *
                                                                //            krn.data[krn.size*(krn.size-1-m) + (krn.size-1-n)];
                        }
                        }
                    }
                }

                //qDebug() << i << " " << j << " " << k << " sum: " << sum;
                image_out[M*N*k+M*j+i]=round(sum*255.0);
                //if(image_out[M*N*k+M*j+i] != 0)
                    //qDebug() << i << " " << j << " " << k << " sum: " << sum;
                //padded_out[w_fft*j+i] = sum;
            }
        }
    }
    /*
    float knl[3] = {0.3302,0.3395,0.3302};
    int len = 3;
    uint64_t counter=0;
    for(uint64_t i = 0; i < M; ++i)
    {
        for ( uint64_t j = 0 ; j < N ; ++j)
        {
            for ( uint64_t k = 0; k < Z; ++k)
            {
                counter++;
                float sum = 0;
                uint64_t rowIndex = j - krn.k;
                uint64_t colIndex = i - krn.k;
                uint64_t zIndex = k - krn.k;

                for( uint64_t m = 0; m < len; ++m)
                {
                    uint64_t mm = krn.size - 1 - m;
                    rowIndex = j + m - krn.k;

                    if(rowIndex >=0     && rowIndex < N &&
                       colIndex >= 0    && colIndex < M &&
                       zIndex >= 0      && zIndex < Z)
                    {
                        sum+=(double)image_in[M*N*zIndex+M*rowIndex+colIndex]*
                                knl[mm];

                    }
                }
                for( uint64_t n = 0; n < len; ++n)
                {
                    uint64_t nn = krn.size - 1 - n;
                    colIndex = i + n - krn.k;

                    if(rowIndex >=0     && rowIndex < N &&
                       colIndex >= 0    && colIndex < M &&
                       zIndex >= 0      && zIndex < Z)
                    {
                        sum+=(double)image_in[M*N*zIndex+M*rowIndex+colIndex]*
                                knl[nn];

                    }
                }
                for( uint64_t p = 0; p < len; ++p)
                {
                    uint64_t pp = krn.size - 1 - p;
                    zIndex = k + p - krn.k;

                    if(rowIndex >=0     && rowIndex < N &&
                       colIndex >= 0    && colIndex < M &&
                       zIndex >= 0      && zIndex < Z)
                    {
                        sum+=(double)image_in[M*N*zIndex+M*rowIndex+colIndex]*
                                knl[pp];

                    }
                }
                image_out[M*N*k+M*j+i]=round(sum*255.0);
                //if(image_out[M*N*k+M*j+i] != 0)
                    //qDebug() << i << " " << j << " " << k << " sum: " << sum;
                if(counter % NPixel/10 == 0)
                    qDebug() << "10%";

            }
        }
    }*/
    mmap_image_in.close();
    conv_out_file.close();
    qDebug() << "convolution ready";


/*    int len = 3;
    int i = 0;
    int l = 0 ;

    //for ( int i = 0 ; i < z_fft; ++i)       // Z
    {
        for( int j = 0; j < h_fft; ++j)     // Y
        {
            for ( int k = 0; k < w_fft; ++k)  // X
            {
                double img_sum = 0;
                //for ( int l = 0; l < len; ++l)
                {
                    for(int m = 0; m < len; ++m)
                    {
                        for( int n = 0; n < len; ++n)
                        {
                            int s = k - krn.k;
                            int t = j - krn.k;
                            //int r = k - krn.k;
                            int ind1 = (M) * ( N) * (i + l) +
                                    (M) * (j + m) + ( k + n );
                            int ind2 = ( len* len) * l + (len) * m + n;

                            img_sum += padded_in[ind1] * krn.data[ind2];

                            //qDebug() << ind1 << "  " << ind2 << " " <<img_sum;

                        }
                    }
                }
                //qDebug() << "img_sum" << img_sum;
                padded_out[h_fft * w_fft * i + w_fft * j + k] = (double) img_sum / (double)(len*len*len);
            }
        }
    }
*/

    /*for(int i = 0; i < NPixelPadded; ++i)
    {
        //qDebug() << padded_out[i];
    }
    for(int i = 0; i < NPixel; ++i)
    {
        //qDebug() << image_out[i];
    }*/
}

void Reconstructor::doWork()
{
    qDebug() <<" started Reconstructor worker";
    bool repeat =true;

    while(repeat && curr_filter<=max_filter)
    {
        repeat=false;
        setArray();
        if(runConvolution)
            Convolution2();
            //Convolution();
        map8bit();
        outputTIFF();
        if(FilterInput)
            repeat=true;
        curr_filter++;
        initData(curr_filter);
    }


    emit finished();
}
void Reconstructor::init(bool getKernelVector, int current_filter)
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
        if(getKernelVector)
        {
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
        }
    }

    qDebug()<< "reconstructor init ready";
}

void Reconstructor::setArray()
{
    qDebug() << "start set array";

    for(int i = 0; i < dimensions; ++i)
    {
        maxPixels*=image_size[i];
        qDebug()<<"set Array image_size[i]: " << image_size[i];
        //maxPixels*=(image_max[i]+1);
    }
    std::vector<Coordinates>::iterator it;

    for( it = xyz.begin(); it != xyz.end(); ++it )
    {
        Coordinates c = *it;
        for(int i = 0; i < dimensions; ++i)
        {
            //qDebug() << c.get(i);
            c.set(i, c.get(i)-image_min[i]);
            //qDebug() << c.get(i);
        }

    }
    qDebug() << "maxPixels: " << maxPixels;
    boost::iostreams::mapped_file file;
    boost::iostreams::mapped_file_params params;

    remove(tiff_temp_file);
    params.path = tiff_temp_file;

    params.mode = (std::ios::out | std::ios::in);
    params.new_file_size =  maxPixels*sizeof(uint8);

    file.open(params);

    uint8 *array =static_cast<uint8*>((void*)file.data());

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
       qDebug()<< "found " << xyz.size() << " coordinates";

       for( it = xyz.begin(); it < xyz.end(); ++it)
       {
           Coordinates c = *it;
           uint64_t lin_index = linearIndex(c);

           array[lin_index]+=1;
           //qDebug() << lin_index << "  " << array[lin_index];
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
       qDebug() << " array init ready!";

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
    sdmixer::log(sdm, "started convolution");
    int img_sizeX = image_size[0];//image_max[0] + 1;
    int img_sizeY = image_size[0];//image_max[1] + 1;
    int img_sizeZ = image_size[0];//image_max[2] + 1;
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
    sdmixer::log(sdm, "convolution ready!");

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
    //qDebug() << "linINdex of: " << c.get(0) << "  " << c.get(1) << " " << c.get(2);
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
void Reconstructor::getHeader(QString header)
{
    header = header.replace("#", "");

    QDomDocument qd;

    qd.setContent(header);
    qDebug() << header;


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
    qDebug() << "dimensions: " << dimensions;
    qDebug() << "max Values from config";
    qDebug() << min_maxValues.min_x << "  " << min_maxValues.max_x;
    qDebug() << min_maxValues.min_y << "  " << min_maxValues.max_y;
    qDebug() << min_maxValues.min_z << "  " << min_maxValues.max_z;

    minMaxDefined=true;

}

Reconstructor::Reconstructor(sdmixer *s, QString xyz_file)
{
    this->sdm=s;
    qDebug() << "Reconstructor loading file " << xyz_file;


    QFile f(xyz_file);
    f.open(QIODevice::ReadOnly| QIODevice::Text);

    QTextStream in(&f);
    QString line = in.readLine();

    getHeader(line);

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
        qDebug() << "v.size() " << v.size();

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

    qDebug() << "this is a filter_file: " << FilterInput;
    if(FilterInput)
        initData(1);
    else
        initData(0);
    /*init(false);
    std::vector<double> input;
    int rawDataCols=0;
    int rawDataRows=0;

    std::ifstream ifs(xyz_file.toStdString());
    std::string firstLine, secondLine;
    getline (ifs, firstLine);
    getline (ifs, secondLine);
    ifs.close();
    //determine column number from second line
    std::stringstream countColsStream(secondLine);
    double dd;
    while (countColsStream >> dd)
    {
        ++rawDataCols;
    }


    QFile f(xyz_file);
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

    std::vector<double>::iterator it;

    for( it = input.begin(); it != input.end(); ++it)
    {
        double d = *it;
        Coordinates c;

        //for(int i = 0; i < )
    }*/

}
void Reconstructor::initData(int current_filter)
{
    int limage_min[3]={0};
    int limage_max[3]={0};
    int limage_size[3]={0};

    maxPixels=1;
    binning[0] = sdm->getReconstructor_xyBinning();
    binning[1] = sdm->getReconstructor_xyBinning();
    binning[2] = sdm->getReconstructor_zBinning();
    if(!dimDefined)
    {
        dimensions = sdm->getCurrentDimensions();
        qDebug()<<"dimensions_ " << dimensions;
    }


    init(true, current_filter);


    tiff_out_file = output_dir.append("/");
    tiff_out_file = tiff_out_file.append(input_base_name);
    if(current_filter != 0)
    {
        tiff_out_file = tiff_out_file.append("_ch");
        tiff_out_file = tiff_out_file.append(QString::number(current_filter));
    }
    tiff_out_file.append(".tif");

    qDebug() << output_dir;


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

    initData(current_filter);

}


void Reconstructor::outputTIFF()
{
    uint16 spp, bpp, photo, res_unit;

    spp = 1; /* Samples per pixel */
    //3 == rgb, 4 == alpha channel
    bpp = 8; /* Bits per sample */
    photo = PHOTOMETRIC_MINISBLACK;

    qDebug() << "started OutputTIFF";
    sdmixer::log(sdm, "started OutputTIFF");


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
    sdmixer::log(sdm, "TIFF out ready!!");

}
