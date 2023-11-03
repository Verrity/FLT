#include "FLT_BaseFilter.h"
#include <iostream>

FLT_BaseFilter::FLT_BaseFilter()
{

}

FLT_BaseFilter::~FLT_BaseFilter()
{
    if (h != nullptr)               delete[] h;
    if (h_magnitude != nullptr)     delete[] h_magnitude;
    if (h_phase != nullptr)         delete[] h_phase;
    if (h_attenuation != nullptr)   delete[] h_attenuation;
    if (freq_match != nullptr)      delete[] freq_match;
    if (w != nullptr)               delete[] w;
    if (conv_frames != nullptr)     delete[] conv_frames;
    if (h_fft != nullptr)           fftw_free(h_fft);
    if (mul_frames_fft != nullptr)  fftw_free(mul_frames_fft);
    fftw_destroy_plan(forward_signal_fft);
    fftw_destroy_plan(backward_signalF);
}

FLT_BaseFilter::Frame::Frame()
{
}

FLT_BaseFilter::Frame::~Frame() {
    if (data != nullptr)       delete[] data;
    if (data_fft != nullptr)   fftw_free(data_fft);
}

bool FLT_BaseFilter::fft_filtrate(Frame& frame) {
    pFrameData = frame.data;
    pFrameDataFFT = frame.data_fft;
    fftw_execute_dft_r2c(forward_signal_fft, frame.data, frame.data_fft);

    for (int i = 0; i < fft_size / 2 + 1; i++) {
        double& a = frame.data_fft[i][0];
        double& b = frame.data_fft[i][1];
        double& c = h_fft[i][0];
        double& d = h_fft[i][1];

        mul_frames_fft[i][0] = ((a * c) - (b * d)) / (fft_size);
        mul_frames_fft[i][1] = ((a * d) + (b * c)) / (fft_size);
        //printf("%-5d)[%d] %f\n", i, 0, mul_frames_fft[i][0]);
        //printf("%-5d)[%d] %f\n\n", i, 1, mul_frames_fft[i][1]);
    }

    fftw_execute_dft_c2r(backward_signalF, mul_frames_fft, frame.data);

    return false;
}

int FLT_BaseFilter::filtrate(double* in, int length, double*& out, int accurancy, bool tails) {
    // Check Parameters ---------------------------------------------
    if (length <= 0) {
        error_code = FILTER_ERROR_LENGTH;
        return 0;
    }
    if (!check_accurancy(accurancy))
        return 0;
    // --------------------------------------------------------------

    // Calc Parameters ----------------------------------------------
    if ( // Если массивы под эти параметры не были рассчитаны рассчитать новые
        (filtrate_length != length) ||
        (filtrate_accurancy != accurancy)
        ) 
    { 
        filtrate_length = length;
        filtrate_accurancy = accurancy;

        this->accurancy = accurancy;
        add_min = N - 1;
        fft_size = 1;
        while (fft_size < (length + add_min)) // fft_size < (длина сигнала + минимум добавочных элементов)
        {
            fft_size = fft_size << 1;
        }
        fft_size << accurancy;

        if (h != nullptr)               delete[] h;
        if (freq_match != nullptr)      delete[] freq_match;
        if (h_fft != nullptr)           fftw_free(h_fft);
        if (mul_frames_fft != nullptr)  fftw_free(mul_frames_fft);
        fftw_destroy_plan(forward_signal_fft);
        fftw_destroy_plan(backward_signalF);

        h = new double[fft_size];
        h_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (fft_size / 2 + 1));
        mul_frames_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (fft_size / 2 + 1) * 2);
        freq_match = new double[fft_size];

        frame1.init(N, length, fft_size);

        forward_signal_fft = fftw_plan_dft_r2c_1d(fft_size, pFrameData, pFrameDataFFT, FFTW_ESTIMATE);
        backward_signalF = fftw_plan_dft_c2r_1d(fft_size, mul_frames_fft, pFrameData, FFTW_ESTIMATE);

        if (window) {
            if (w != nullptr)   delete[] w;
            w = new double[fft_size];
            calc_window();
        }
        calc_h();
        calc_h_fft();
    }
    // --------------------------------------------------------------

    frame1.setData(in, 0, length);
    fft_filtrate(frame1);

    std::cout << "frame1.add " << frame1.add << std::endl;
    std::cout << "frame1.add_min " << frame1.add_min << std::endl;
    std::cout << "frame1.add_min_left " << frame1.add_min_left << std::endl;
    std::cout << "frame1.add_other " << frame1.add_other << std::endl;
    std::cout << "frame1.add_other_left " << frame1.add_other_left << std::endl;
    std::cout << "frame1.add_other_right " << frame1.add_other_right << std::endl;
    std::cout << "frame1.data_size " << frame1.data_size << std::endl;
    std::cout << "fft_size " << fft_size << std::endl;

    int len = 0;
    if (tails) {
        len = length + add_min;
        out = new double[len];
        for (int i = 0; i < len; i++)
        {
            //out[i] = frame1.data[frame1.add_other_left + i];
            out[i] = frame1.data[i];
            //printf("%d) %f\n", i, frame1.data[frame1.add_other_left + i]);
            //printf("%f\n", frame1.data[frame1.add_other_left + i]);
        }
    }
    else {
        len = length;
        out = new double[len];
        for (int i = 0; i < len; i++)
        {
            out[i] = frame1.data[frame1.add_min_left + i];
        }
    }
    return len;
}

bool FLT_BaseFilter::setIrLowpassR1B1(int N, double fd, double band, int window) {
    // ---------------- Check parameters 
    if (
        !check_N(N) ||
        !check_fd(fd) ||
        !check_window(window) ||
        !check_band_LowpassR1B1(band, fd)
        )
        return 0;
    this->N = N;
    this->fd = fd;
    this->window = window;
    bands.resize(1);
    bands.at(0) = band;
    
    type = 1;
    method = 1;
    bands_count = 1;

    return true;
}

void FLT_BaseFilter::Frame::init(int N, int data_size, int fft_size) {
    this->N = N;
    this->data_size = data_size;
    this->fft_size = fft_size;

    if (data != nullptr)       delete[] data;
    if (data_fft != nullptr)   fftw_free(data_fft);

    data = new double[fft_size];
    data_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (fft_size / 2 + 1));

    add = fft_size - data_size;
    add_min = N - 1;
    add_min_left = add_min / 2;
    add_other = fft_size - data_size - add_min;
    add_other_left = add_other / 2;
    add_other_right = add_other - add_other_left;
}

bool FLT_BaseFilter::convolFull(Frame& frame1, Frame& frame2) {
    //if (frame1.data_size != frame2.data_size) {
    //    printf("Convolution ERROR (convolFull) array sizes are different\n");
    //    return 0;
    //}
    int size = frame1.data_size;
    for (int i = 0; i < size; i++)
        conv_frames[i] = frame1.data[i];
    for (int i = 0; i < add_min; i++)
        conv_frames[size + i] = frame1.data[size + i] + frame2.data[i];
    for (int i = 0; i < size; i++)
        conv_frames[size + add_min + i] = frame2.data[add_min + i];
    return true;
}

void FLT_BaseFilter::convolDifferent(Frame& frame1, Frame& frame2) {

    int pos = 0;

    int frame1_pre1 = frame1.add_other_left;
    int frame1_len1 = frame1.add_min_left;
    int frame1_pre2 = frame1_pre1 + frame1_len1;
    int frame1_len2 = abs(frame1.data_size - frame1.add_min_left);
    int frame1_pre3 = frame1_pre2 + frame1_len2;
    int len3 = frame1.add_min_left * 2;

    int frame2_pre3 = frame2.add_other_left;
    int frame2_pre4 = frame2_pre3 + len3;
    int frame2_len4 = abs(frame2.data_size - frame2.add_min_left);
    int frame2_pre5 = frame2_pre4 + frame2_len4;
    int frame2_len5 = frame2.add_min_left;

    // -- 1 --
    for (int i = 0; i < frame1_len1; i++)
        conv_frames[i] = frame1.data[frame1_pre1 + i];
    pos = frame1_len1;

    // -- 2 --
    for (int i = 0; i < frame1_len2; i++)
        conv_frames[pos + i] = frame1.data[frame1_pre2 + i];
    pos += frame1_len2;

    // -- 3 --
    for (int i = 0; i < len3; i++)
        conv_frames[pos + i] = frame1.data[frame1_pre3 + i] + frame2.data[frame2_pre3 + i];
    pos += len3;

    // -- 4 --
    for (int i = 0; i < frame2_len4; i++)
        conv_frames[pos + i] = frame2.data[frame2_pre4 + i];
    pos += frame2_len4;

    // -- 5 --
    for (int i = 0; i < frame2_len5; i++)
        conv_frames[pos + i] = frame2.data[frame2_pre5 + i];
}

void FLT_BaseFilter::Frame::setData(double* arr, int begin, int size)
{
    //if (size != this->data_size) {
    //    if (data != nullptr) {
    //        delete[] data;
    //    }
    //    data = new double[fft_size];
    //    this->data_size = size;
    //}

    for (int i = 0; i < fft_size; i++) {
        if (i < size)  // Сигнал
            data[i] = arr[begin + i];
        else        // Дополнить нулями для БПФ
            data[i] = 0;
    }

    add_other = fft_size - size - (N - 1);
    add_other_left = add_other / 2;
    add_other_right = add_other - add_other_left;
}

void FLT_BaseFilter::Frame::setData(Frame& frame)
{
    this->setData(frame.data, 0, frame.data_size);
}

void FLT_BaseFilter::Frame::switchData(Frame& toFrame, Frame& fromFrame)
{
    double* temp = toFrame.data;
    toFrame.data = fromFrame.data;
    fromFrame.data = temp;

}

int FLT_BaseFilter::calc_IrLowpassR1B1()
{
    double& B1 = bands.at(0);
    for (int i = 0; i < fft_size; i++)
        if (i < (N - 1) / 2) {
            h[i] = sin(2 * FILTER_PI * B1 / fd * (i - (N - 1) / 2)) / (FILTER_PI * (i - (N - 1) / 2));
            h[N - 1 - i] = h[i];
        }
        else
            if (i >= N)
                h[i] = 0;
    h[(N - 1) / 2] = 2 * (B1 / fd); // central

    // if window
    if (window)
        for (int i = 0; i < N; i++)
            h[i] = h[i] * w[i];
    return 0;
}



void FLT_BaseFilter::calc_h()
{
    switch (type) // Lowpass, Bandpas...
    {
    case 1: // Lowpass
        switch (method) // Weighting, Разложения
        {
        case 1: // Weighting
            switch (bands_count) // Количество полос частот
            {
            case 1: // Одна полоса
                calc_IrLowpassR1B1(); // Lowpass, Метод взвешивания, 1 полоса
            default:
                break;
            }
            
        default:
            break;
        }
    default:
        break;
    }
}

void FLT_BaseFilter::calc_h_fft()
{
    fftw_plan forward_h = fftw_plan_dft_r2c_1d(fft_size, h, h_fft, FFTW_ESTIMATE);
    fftw_execute(forward_h);
    fftw_destroy_plan(forward_h);
}

void FLT_BaseFilter::calc_h_fft_mag_ph_att() {
    // Расчет и выполнение БПФ для Имп. хар-ки
    fftw_plan forward_h = fftw_plan_dft_r2c_1d(fft_size, h, h_fft, FFTW_ESTIMATE);
    fftw_execute(forward_h);
    // АЧХ, ФЧХ, Затухание
    for (int i = 0; i < fft_size / 2 + 1; i++) {
        h_magnitude[i] = calc_magnitude(h_fft[i][0], h_fft[i][1]);

        h_phase[i] = calc_phase(h_fft[i][0], h_fft[i][1]);
        h_attenuation[i] = 20 * log10(1 / h_magnitude[i]);
    };
    fftw_destroy_plan(forward_h);
}

void FLT_BaseFilter::calc_window() {
    switch (window)
    {
    case 1:     // Hamming
        for (int i = 0; i < fft_size; i++)
            if (i < (N - 1) / 2) {
                w[i] = 0.54 + 0.46 * cos(2 * FILTER_PI / N * (i - (N - 1) / 2));
                w[N - 1 - i] = w[i];
            }
            else
                if (i >= N)
                    w[i] = 0;
        w[(N - 1) / 2] = 1;
        break;
    case 2:     // Blackman
        for (int i = 0; i < fft_size; i++)
            if (i < (N - 1) / 2) {
                w[i] = 0.42 + 0.5 * cos(2 * FILTER_PI / N * (i - (N - 1) / 2)) + 0.08 * cos(4 * FILTER_PI / N * (i - (N - 1) / 2));
                w[N - 1 - i] = w[i];
            }
            else
                if (i >= N)
                    w[i] = 0;
        w[(N - 1) / 2] = 1;
        break;
    case 3:     // Hann
        for (int i = 0; i < fft_size; i++)
            if (i < (N - 1) / 2) {
                w[i] = double(1) / 2 * (1 + cos(2 * FILTER_PI / N * (i - (N - 1) / 2)));
                w[N - 1 - i] = w[i];
            }
            else
                if (i >= N)
                    w[i] = 0;
        w[(N - 1) / 2] = 1;
    }
}

bool FLT_BaseFilter::check_N(int N)
{
    if ((N < 17) || ((N % 2) == 0)) {
        error_code = FILTER_ERROR_N;
        return false;
    }
    return true;
}

bool FLT_BaseFilter::check_fd(double fd)
{
    if (fd <= 16) {
        error_code = FILTER_ERROR_FD;
        return false;
    }
    return true;
}

bool FLT_BaseFilter::check_accurancy(int accurancy) {
    if (accurancy < 1) {
        error_code = FILTER_ERROR_ACCURANCY;
        return false;
    }
    return true;
}


bool FLT_BaseFilter::check_window(int window)
{
    if ((window < 0) || (window > 3)) {
        error_code = FILTER_ERROR_WINDOW;
        return false;
    }
    return true;
}

bool FLT_BaseFilter::check_band_LowpassR1B1(double band, double fd)
{
    if ((band >= fd / 2) || (band <= 0)) {
        error_code = FILTER_ERROR_BAND;
        return false;
    }
    return true;
}

double FLT_BaseFilter::calc_magnitude(double& real_value, double& complex_value) {
    return sqrt(pow(real_value, 2) + pow(complex_value, 2));
}
double FLT_BaseFilter::calc_phase(double& real_value, double& complex_value) {
    return atan2(complex_value, real_value);
}

int FLT_BaseFilter::get_fft_size()
{
    return fft_size;
}

int FLT_BaseFilter::get_N()
{
    return N;
}

int FLT_BaseFilter::get_window()
{
    return window;
}

int FLT_BaseFilter::get_h(double* &h)
{
    h = new double[N];
    for (int i = 0; i < N; i++)
    {
        h[i] = this->h[i];
    }
    return N;
}

int FLT_BaseFilter::get_h_magnitude(double* &magnitude, bool symmetrical)
{
    int fft_short = fft_size / 2 + 1;
    if (symmetrical) {
        magnitude = new double[fft_size];

        for (int i = 0; i < fft_short; i++)
            magnitude[i] = h_magnitude[i];
        for (int i = 0; i < fft_short - 2; i++)
            magnitude[fft_short + i] = h_magnitude[fft_short - 2 - i];

        return fft_size;
    }
    else {
        magnitude = new double[fft_short];

        for (int i = 0; i < fft_short; i++) {
            magnitude[i] = h_magnitude[i];
        }
        return fft_short;
    }
}

int FLT_BaseFilter::get_h_phase(double* &phase, bool symmetrical)
{
    int fft_short = fft_size / 2 + 1;
    if (symmetrical) {
        phase = new double[fft_size];

        for (int i = 0; i < fft_short; i++)
            phase[i] = h_phase[i];
        for (int i = 0; i < fft_short - 2; i++)
            phase[fft_short + i] = h_phase[fft_short - 2 - i];
        return fft_size;
    }
    else {
        phase = new double[fft_short];

        for (int i = 0; i < fft_short; i++) {
            phase[i] = h_phase[i];
        }
        return fft_short;
    }
}

int FLT_BaseFilter::get_h_attenuation(double* &attenuation, bool symmetrical)
{
    int fft_short = fft_size / 2 + 1;
    if (symmetrical) {
        attenuation = new double[fft_size];

        for (int i = 0; i < fft_short; i++)
            attenuation[i] = h_attenuation[i];
        for (int i = 0; i < fft_short - 2; i++)
            attenuation[fft_short + i] = h_attenuation[fft_short - 2 - i];
        return fft_size;
    }
    else {
        attenuation = new double[fft_short];

        for (int i = 0; i < fft_size / 2 + 1; i++) {
            attenuation[i] = h_attenuation[i];
        }
        return fft_short;
    }
}

int FLT_BaseFilter::get_freq_match(double* &freq_match, bool symmetrical)
{
    int fft_short = fft_size / 2 + 1;
    if (symmetrical) {
        freq_match = new double[fft_size];

        for (int i = 0; i < fft_size; i++)
            freq_match[i] = this->freq_match[i];
        return fft_size;
    }
    else {
        freq_match = new double[fft_short];

        for (int i = 0; i < fft_short; i++) {
            freq_match[i] = this->freq_match[i];
        }
        return fft_short;
    }
}

int FLT_BaseFilter::get_w(double* &w)
{
    for (int i = 0; i < N; i++)
        w[i] = this->w[i];
    return N;
}

int FLT_BaseFilter::get_error_code()
{
    return error_code;
}
