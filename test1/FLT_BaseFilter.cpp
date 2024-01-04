#include "FLT_BaseFilter.h"

FLT_BaseFilter::FLT_BaseFilter()
{

}

FLT_BaseFilter::~FLT_BaseFilter()
{
    if (isMinAllocated)
        free_min();
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
    }

    fftw_execute_dft_c2r(backward_signalF, mul_frames_fft, frame.data);

    return false;
}

void FLT_BaseFilter::init_min(int N, int fd, int fft_size, int frame_size)
{
    h = new double[fft_size];
    h_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (fft_size / 2 + 1));
    mul_frames_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (fft_size / 2 + 1) * 2);

    frame1.init(N, frame_size, fft_size);
    frame2.init(N, frame_size, fft_size);

    forward_signal_fft = fftw_plan_dft_r2c_1d(fft_size, pFrameData, pFrameDataFFT, FFTW_ESTIMATE);
    backward_signalF = fftw_plan_dft_c2r_1d(fft_size, mul_frames_fft, pFrameData, FFTW_ESTIMATE);

    if (window) {
        if (w != nullptr)   delete[] w;
        w = new double[fft_size];
        calc_window();
    }
    calc_h();
    calc_h_fft();
    isMinAllocated = true;
}

void FLT_BaseFilter::free_min()
{
    if (h != nullptr)               delete[] h;                 h = nullptr;
    if (h_fft != nullptr)           fftw_free(h_fft);           h_fft = nullptr;
    if (w != nullptr)               fftw_free(h_fft);           h_fft = nullptr;
    if (mul_frames_fft != nullptr)  fftw_free(mul_frames_fft);  mul_frames_fft = nullptr;
    fftw_destroy_plan(forward_signal_fft);
    fftw_destroy_plan(backward_signalF);
    isMinAllocated = false;
   
}

bool FLT_BaseFilter::filtrate(double* const in, int length, int accurancy) {
    // Check Parameters ---------------------------------------------
    if (length <= 0) {
        error_code = FILTER_ERROR_LENGTH;
        return false;
    }
    if (!check_accurancy(accurancy))
        return false;
    // --------------------------------------------------------------

    if (// Если массивы под эти параметры не были рассчитаны рассчитать новые
        (signal_size != length) ||
        (this->accurancy != accurancy)
        )
    {
        signal_size = length;
        this->accurancy = accurancy;

        add_min = N - 1;
        fft_size = 1;
        while (fft_size < (length + add_min)) // fft_size < (длина сигнала + минимум добавочных элементов)
        {
            fft_size = fft_size << 1;
        }
        
        fft_size = fft_size << accurancy;
        if (isMinAllocated)
            free_min();
        init_min(this->N, this->fd, this->fft_size, length);
    }

    if (in != nullptr) {
        frame1.setData(in, 0, length);
        fft_filtrate(frame1);

        for (int i = 0; i < length; i++)
            in[i] = frame1.data[add_min2 + i];
    }
    return true;
}

int FLT_BaseFilter::filtrateT(double* const in, int length, double*& output, int accurancy) {
    // Check Parameters ---------------------------------------------
    if (length <= 0) {
        error_code = FILTER_ERROR_LENGTH;
        return 0;
    }
    if (!check_accurancy(accurancy))
        return 0;
    // --------------------------------------------------------------

    if (// Если массивы под эти параметры не были рассчитаны рассчитать новые
        (signal_size != length) ||
        (this->accurancy != accurancy)
        )
    {
        signal_size = length;
        this->accurancy = accurancy;

        add_min = N - 1;
        fft_size = 1;
        while (fft_size < (length + add_min)) // fft_size < (длина сигнала + минимум добавочных элементов)
        {
            fft_size = fft_size << 1;
        }
        fft_size = fft_size << accurancy;

        if ((in != nullptr) || (output != nullptr)) {
            if (isMinAllocated)
                free_min();
            init_min(this->N, this->fd, this->fft_size, length);
        }
    }

    if ((in != nullptr) || (output != nullptr)) {
        frame1.setData(in, 0, length);
        fft_filtrate(frame1);

        int len = length + add_min;
        output = new double[len];
        for (int i = 0; i < len; i++)
            output[i] = frame1.data[i];

        return len;
    }
}

bool FLT_BaseFilter::setIrBandpassR2B2(int N, double fd, double band1, double band2, double band3, double band4, int window)
{
    // ---------------- Check parameters 
    if (
        !check_N(N) ||
        !check_fd(fd) ||
        !check_window(window) ||
        !check_bands(band1, band2, band3, band4, fd)
        )
        return 0;
    this->N = N;
    this->fd = fd;
    this->add_min = N - 1;
    this->add_min2 = add_min / 2;
    this->window = window;
    bands.resize(4);
    bands.at(0) = band1;
    bands.at(1) = band2;
    bands.at(2) = band3;
    bands.at(3) = band4;

    type = 3;
    method = 2;
    bands_count = 4;

    return true;
}

int FLT_BaseFilter::measureAttenuation(double*& pointer, int length, int accurancy, double f_low, double step, double f_high, unsigned long &time_ns) 
{
    //   
    if ((f_low <= 0) || (step < f_low) || (f_high <= step) || (f_high > fd / 2)) {
        error_code = FILTER_ERROR_BAND;
        return 0;
    }
    if (!check_accurancy(accurancy))
        return 0;

    double AM_f_low = f_low;
    double AM_f_high = f_high;
    double AM_step = step;
    int AM_len = (AM_f_high - AM_f_low) / AM_step + 1;

    fft_size = 1;
    while (fft_size < (length + add_min)) // fft_size < (длина сигнала + минимум добавочных элементов)
    {
        fft_size = fft_size << 1;
    }
    fft_size << accurancy;

    if (isMinAllocated)
        free_min();
    init_min(this->N, this->fd, this->fft_size, 0);

    double* signal = new double[fft_size];
    pointer = new double[AM_len];
    fftw_complex* signalFFT = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (fft_size / 2 + 1));
    fftw_complex* signalFiltratedFFT = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (fft_size / 2 + 1) * 2);
    //   
    fftw_plan forward_h = fftw_plan_dft_r2c_1d(fft_size, signal, signalFFT, FFTW_ESTIMATE);
    fftw_plan backward_signalF = fftw_plan_dft_c2r_1d(fft_size, signalFiltratedFFT, signal, FFTW_ESTIMATE);

    double averageOUT = 0;
    double averageIN = 0.774596669241; // Опорное значение 0 db [0.775 Вольт]
    double sum = 0;
    double H = 0;
    int iterator = 0;

    for (double f = AM_f_low; f <= AM_f_high; f += AM_step) {
        averageOUT = 0;
        sum = 0;
        //     f
        for (int j = 0; j < fft_size; j++)
            signal[j] = averageIN * sqrt(2) * sin(2 * FILTER_PI * j * f / fd);
        // 
        if (f == AM_f_low)
            startTimer();
        fftw_execute_dft_r2c(forward_h, signal, signalFFT);
        // Умнодить БПФ сигнала и БПФ Имп. хар-ки. (фильтрация
        for (int i = 0; i < fft_size / 2 + 1; i++) {
            double& a = signalFFT[i][0];
            double& b = signalFFT[i][1];
            double& c = h_fft[i][0];
            double& d = h_fft[i][1];

            signalFiltratedFFT[i][0] = ((a * c) - (b * d)) / fft_size;
            signalFiltratedFFT[i][1] = ((a * d) + (b * c)) / fft_size;
        }
        // ОБПФ
        fftw_execute_dft_c2r(backward_signalF, signalFiltratedFFT, signal);
        if (f == AM_f_low) {
            stopTimer();
            time_ns = m_int_durationTime_ns / fft_size;
        }
        // Сумма модулей сигналов
        for (int j = 0; j < fft_size; j++)
            sum = sum + std::abs(signal[j]);

        averageOUT = sum / fft_size;            // Среднее значение
        H = averageOUT / averageIN;             // Коэффициент передачи фильтра
        pointer[iterator] = 20 * log10(H);      // Затухание в dBu
        iterator++;
    }
    delete[] signal;
    fftw_free(signalFFT);
    fftw_free(signalFiltratedFFT);
    fftw_destroy_plan(forward_h);
    fftw_destroy_plan(backward_signalF);
    return AM_len;
}

void FLT_BaseFilter::Frame::init(int N, int data_size, int fft_size) {
    if (data_size > 0){
        this->N = N;
        this->data_size = data_size;
        this->fft_size = fft_size;

        if (data != nullptr)       delete[] data;
        if (data_fft != nullptr)   fftw_free(data_fft);

        data = new double[fft_size];
        data_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (fft_size / 2 + 1));
    }
}

void FLT_BaseFilter::Frame::setData(double* arr, int begin, int size)
{
    this->data_size = size;

    for (int i = 0; i < fft_size; i++) {
        if (i < size)  // Сигнал
            data[i] = arr[begin + i];
        else        // Дополнить нулями для БПФ
            data[i] = 0;
    }
}

void FLT_BaseFilter::Frame::switchData(Frame& toFrame, Frame& fromFrame)
{
    double* temp = toFrame.data;
    toFrame.data = fromFrame.data;
    fromFrame.data = temp;
}

int FLT_BaseFilter::calc_IrBandpassR2B2()
{
    double& border_stop1 = bands.at(0);
    double& border_pass1 = bands.at(1);
    double& border_pass2 = bands.at(2);
    double& border_stop2 = bands.at(3);
    
    double numenator1_1 = 0, numenator1_2 = 0, denominator1 = 0;
    double numenator2_1 = 0, numenator2_2 = 0, denominator2 = 0;
    double multiplier = 0;
    int x1 = 0, x2 = 0;
    // Рассчитать Имп. хар-ку
    for (int i = 0; i < fft_size; i++)
        if (i < (N - 1) / 2) {
            numenator1_1 = cos(2 * FILTER_PI * ((N - 1) / 2 - i) * (border_stop1 / fd));
            numenator1_2 = cos(2 * FILTER_PI * ((N - 1) / 2 - i) * (border_pass1 / fd));
            denominator1 = border_stop1 / fd - border_pass1 / fd;
            x1 = (numenator1_1 - numenator1_2) / denominator1;

            numenator2_1 = cos(2 * FILTER_PI * ((N - 1) / 2 - i) * (border_pass2 / fd));
            numenator2_2 = cos(2 * FILTER_PI * ((N - 1) / 2 - i) * (border_stop2 / fd));
            //denominator2 = border_pass2 / fd - border_stop2 / fd;
            denominator2 = border_stop2 / fd - border_pass2 / fd;
            x2 = (numenator2_1 - numenator2_2) / denominator2;

            multiplier = 1 / (2 * pow(FILTER_PI, 2) * pow((N - 1) / 2 - i, 2));

            h[i] = multiplier * (x1 + x2);
            h[N - 1 - i] = h[i];
        }
        else
            if (i >= N)
                h[i] = 0;

    h[(N - 1) / 2] = (border_pass2 + border_stop2 - border_pass1 - border_stop1) / fd; // средний отсчет

    // if window
    if (window)
        for (int i = 0; i < N; i++)
            h[i] = h[i] * w[i];
    return 0;
}

void FLT_BaseFilter::calc_h()
{
    calc_IrBandpassR2B2(); // Bandpass, Метод разложения, 4 полосы
}

void FLT_BaseFilter::calc_h_fft()
{
    fftw_plan forward_h = fftw_plan_dft_r2c_1d(fft_size, h, h_fft, FFTW_ESTIMATE);
    fftw_execute(forward_h);
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

void FLT_BaseFilter::startTimer()
{
    m_startTime = clock_t::now();
}

void FLT_BaseFilter::stopTimer()
{
    m_endTime = clock_t::now();
    m_durationTime_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(m_endTime - m_startTime);
    m_int_durationTime_ns = m_durationTime_ns.count() /*- m_minimum_execution_time_timer_code*/; //       ""
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
    if (accurancy < 0) {
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

bool FLT_BaseFilter::check_bands(double band1, double band2, double band3, double band4, double fd)
{
    if ((band1 < 0) || (band1 > band2) || (band2 > band3) || (band3 > band4) || (band4 > fd/2)) {
        error_code = FILTER_ERROR_BAND;
        return false;
    }
    return true;
}

double FLT_BaseFilter::calc_magnitude(double& real_value, double& complex_value) {
    return sqrt(pow(real_value, 2) + pow(complex_value, 2));
}

double FLT_BaseFilter::calc_phase(double& real_value, double complex_value) {
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
            magnitude[i] = calc_magnitude(h_fft[i][0], h_fft[i][1]);
        for (int i = 0; i < fft_short - 2; i++)
            magnitude[fft_short + i] = calc_magnitude(h_fft[fft_short - 2 - i][0], h_fft[fft_short - 2 - i][1]);

        return fft_size;
    }
    else {
        magnitude = new double[fft_short];

        for (int i = 0; i < fft_short; i++) {
            magnitude[i] = calc_magnitude(h_fft[i][0], h_fft[i][1]);
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
            phase[i] = calc_phase(h_fft[i][0], h_fft[i][1]);

        for (int i = 0; i < fft_short - 2; i++)                 
            phase[fft_short + i] = calc_phase(h_fft[fft_short - 2 - i][0], -h_fft[fft_short - 2 - i][1]);

        return fft_size;
    }
    else {
        phase = new double[fft_short];

        for (int i = 0; i < fft_short; i++) {
            phase[i] = calc_phase(h_fft[i][0], h_fft[i][1]);
        }
        return fft_short;
    }
}

int FLT_BaseFilter::get_h_attenuation(double* &attenuation, bool symmetrical)
{
    int fft_short = fft_size / 2 + 1;
    double magnitude = 0;
    if (symmetrical) {
        attenuation = new double[fft_size];

        for (int i = 0; i < fft_short; i++) {
            magnitude = calc_magnitude(h_fft[i][0], h_fft[i][1]);
            attenuation[i] = 20 * log10(1 / magnitude);
        }
        for (int i = 0; i < fft_short - 2; i++) {
            magnitude = calc_magnitude(h_fft[fft_short - 2 - i][0], h_fft[fft_short - 2 - i][1]);
            attenuation[fft_short + i] = 20 * log10(1 / magnitude);
        }
        return fft_size;
    }
    else {
        attenuation = new double[fft_short];

        for (int i = 0; i < fft_short; i++) {
            magnitude = calc_magnitude(h_fft[i][0], h_fft[i][1]);
            attenuation[i] = 20 * log10(1 / magnitude);
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
            freq_match[i] = double(i) / fft_size * fd;
        return fft_size;
    }
    else {
        freq_match = new double[fft_short];

        for (int i = 0; i < fft_short; i++) {
            freq_match[i] = double(i) / fft_size * fd;
        }
        return fft_short;
    }
}

int FLT_BaseFilter::get_w(double* &w)
{
    if (!window) {
        error_code = FILTER_ERROR_WINDOW;
        return 0;
    }
        w = new double[N];
        for (int i = 0; i < N; i++)
            w[i] = this->w[i];
        return N;
}

int FLT_BaseFilter::get_error_code()
{
    return error_code;
}
