#include "FLT_FilterFile.h"
#include <iostream>

FLT_FilterFile::FLT_FilterFile()
{
}

FLT_FilterFile::~FLT_FilterFile()
{
}

//int FLT_FilterFile::filtrate(double* in, int length, double* &out, bool tails) {
//    if (length <= 0) {
//        error_code = FILTER_ERROR_LENGTH;
//        return 0;
//    }
//    
//    frame1.setData(in, 0, length);
//    fft_filtrate(frame1);
//
//    std::cout << "frame1.add " << frame1.add << std::endl;
//    std::cout << "frame1.add_min " << frame1.add_min << std::endl;
//    std::cout << "frame1.add_min_left " << frame1.add_min_left << std::endl;
//    std::cout << "frame1.add_other " << frame1.add_other << std::endl;
//    std::cout << "frame1.add_other_left " << frame1.add_other_left << std::endl;
//    std::cout << "frame1.add_other_right " << frame1.add_other_right << std::endl;
//    std::cout << "frame1.data_size " << frame1.data_size << std::endl;
//    std::cout << "fft_size " << fft_size << std::endl;
//    
//    int len = 0;
//    if (tails) {
//        len = length + add_min;
//        out = new double[len];
//        for (int i = 0; i < len; i++)
//        {
//            out[i] = frame1.data[frame1.add_other_left + i];
//            //printf("%f", frame1.data[frame1.add_other_left + i])
//        }
//    }
//    else {
//        len = length;
//        out = new double[len];
//        for (int i = 0; i < len; i++)
//        {
//            out[i] = frame1.data[frame1.add_other_left + frame1.add_min_left + i];
//        }
//    }
//    return len;
//}

int FLT_FilterFile::filtrateBlock(double* in, int length, double* &out, bool tails) {
    if (length <= 0) {
        error_code = FILTER_ERROR_LENGTH;
        return 0;
    }



    // Если сигнал меньше стандартной выборки
    if (length < fft_size) {
        sample_size = length;
    }
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Обработать

    if (conv_frames != nullptr)     delete[] conv_frames;
    conv_frames = new double[sample_size * 2 + frame1.add];
}

//bool FLT_FilterFile::fft_filtrate(Frame& frame) {
//    pFrameData = frame.data;
//    pFrameDataFFT = frame.data_fft;
//    fftw_execute_dft_r2c(forward_signal_fft, frame.data, frame.data_fft);
//
//    for (int i = 0; i < fft_size / 2 + 1; i++) {
//        double& a = frame.data_fft[i][0];
//        double& b = frame.data_fft[i][1];
//        double& c = h_fft[i][0];
//        double& d = h_fft[i][1];
//
//        mul_frames_fft[i][0] = ((a * c) - (b * d)) / (fft_size);
//        mul_frames_fft[i][1] = ((a * d) + (b * c)) / (fft_size);
//    }
//
//    fftw_execute_dft_c2r(backward_signalF, mul_frames_fft, frame.data);
//    return false;
//}

//bool FLT_FilterFile::local_init(int N, double fd, int accurancy, int window) {
//    // ---------------- data_size & fft_size SET 
//    sample_size = N * 45;
//    fft_size = 1;
//    while (fft_size < sample_size)
//    {
//        fft_size = fft_size << 1;
//    }
//    fft_size *= accurancy;
//    add_min = N - 1;
//
//    // ---------------- Array size set
//    if (h != nullptr)               delete[] h;
//    if (h_magnitude != nullptr)     delete[] h_magnitude;
//    if (h_phase != nullptr)         delete[] h_phase;
//    if (h_attenuation != nullptr)   delete[] h_attenuation;
//    if (freq_match != nullptr)      delete[] freq_match;
//
//    if (h_fft != nullptr)           fftw_free(h_fft);
//    if (mul_frames_fft != nullptr)  fftw_free(mul_frames_fft);
//
//    h = new double[fft_size];
//    h_magnitude = new double[fft_size / 2 + 1];
//    h_phase = new double[fft_size / 2 + 1];
//    h_attenuation = new double[fft_size / 2 + 1];
//    freq_match = new double[fft_size];
//
//    frame1.init(N, sample_size, fft_size);
//    frame2.init(N, sample_size, fft_size);
//
//    forward_signal_fft = fftw_plan_dft_r2c_1d(fft_size, pFrameData, pFrameDataFFT, FFTW_ESTIMATE);
//    backward_signalF = fftw_plan_dft_c2r_1d(fft_size, mul_frames_fft, pFrameData, FFTW_ESTIMATE);
//
//    h_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (fft_size / 2 + 1));
//    mul_frames_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (fft_size / 2 + 1) * 2);
//
//    // Расчет окна
//    if (window) {
//        if (w != nullptr)   delete[] w;
//        w = new double[fft_size];
//        calc_window();
//    }
//
//    // ---------------- Array Calc
//
//    // Расчет массива соответствующих частот
//    for (int i = 0; i < fft_size; i++) {
//        freq_match[i] = double(i) / fft_size * fd;
//    }
//
//    return 0;
//}

//bool FLT_FilterFile::setIrLowpassR1B1(int N, double fd, int accurancy, double band, int window) {
//    // ---------------- Check parameters 
//    if (
//        !check_N(N) ||
//        !check_fd(fd) ||
//        !check_accurancy(accurancy) ||
//        !check_window(window) ||
//        !check_band_LowpassR1B1(band, fd)
//        )
//        return 0;
//    this->N = N;
//    this->fd = fd;
//    this->accurancy = accurancy;
//    this->window = window;
//    bands.resize(1);
//    bands.at(0) = band;
//
//    // ---------------- sample_size & fft_size SET 
//    local_init(N, fd, accurancy, window);
//
//    calc_IrLowpassR1B1();
//    calc_h_fft_mag_ph_att();
//
//    return 1;
//}

int FLT_FilterFile::get_sample_size()
{
    return sample_size;
}