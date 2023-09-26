#include "FLT_filter.h"

FILTER FLT_Filter::next_descriptor = 1;
std::unordered_map<FILTER, FLT_Filter*> FLT_Filter::list;

FLT_Filter::FLT_Filter(FILTER& filter, int N, double fd) : N(N), fd(fd)
{
    if (list.empty()) {
        descriptor = 0;
        next_descriptor = 1;
    }
    descriptor = next_descriptor;
    next_descriptor++;
    list[descriptor] = this;
    filter = rand() * fd * N;
}

FLT_Filter::~FLT_Filter()
{

}

FILTER FLT_Filter::get_id() {
    return descriptor;
}

FLT_Filter& FLT_Filter::get_ref()
{
    return *this;
}


FLT_Filter* FLT_Filter::get_pointer(FILTER id)
{
    auto it = FLT_Filter::list.find(id);
    if (it != FLT_Filter::list.end()) {
        FLT_Filter* ptr = it->second;
        return ptr;
    }
    else {
        return nullptr;
    }
}

int FLT_Filter::check_params(const char* filter, int N, int fd, int B1, int B2, int B3, int B4, int window)
{
    if (errorCode)
        return FILTER_ERROR_PREV;
    if ((N < 17) || ((N % 2) == 0))
        return FILTER_ERROR_N;
    if (fd <= 16)
        return FILTER_ERROR_FD;
    if ((window < 0) || (window > 3))
        return FILTER_ERROR_Window;

    //switch (type)
    //{
    //case 1: // Lowpass
    //    if (real == 1)  // Realisation 1
    //        if ((B1 >= fd / 2) || (B1 <= 0)) {
    //            errorCode = FILTER_ERROR_BAND;
    //        }
    //    if (real == 2)
    //        if ((B1 <= 0) || (B1 >= B2) || (B2 >= fd / 2)) {
    //            errorCode = FILTER_ERROR_BAND;
    //        }
    //    break;
    //case 2:
    //    if (real == 1)
    //        if ((B1 >= fd / 2) || (B1 <= 0)) {
    //            errorCode = FILTER_ERROR_BAND;
    //        }
    //    if (real == 2)
    //        if ((B2 <= 0) || (B2 >= B1) || (B1 >= fd / 2)) {
    //            errorCode = FILTER_ERROR_BAND;
    //        }
    //    break;
    //case 3:
    //    if (real == 1)
    //        if ((B1 <= 0) || (B1 >= B2) || (B2 >= fd / 2)) {
    //            errorCode = FILTER_ERROR_BAND;
    //        }
    //    if (real == 2)
    //        if ((BS1 <= 0) || (BS1 >= BP1) || (BP1 >= BP2) || (BP2 >= BS2) || (BS2 >= fd / 2)) {
    //            errorCode = FILTER_ERROR_BAND;
    //        }
    //    break;
    //case 4:
    //    if (real == 1)
    //        if ((f_low <= 0) || (f_low >= f_high) || (f_high >= fd / 2)) {
    //            errorCode = FILTER_ERROR_BAND;
    //            return 0;
    //        }
    //    if (real == 2)
    //        if ((BP1 <= 0) || (BP1 >= BS1) || (BS1 >= BS2) || (BS2 >= BP2) || (BP1 >= fd / 2)) {
    //            errorCode = FILTER_ERROR_BAND;
    //            return 0;
    //        }
    //default:
    //    break;
    //}
    return 0;
}

void FLT_Filter::createIRLowpassR1B1(const char* filter, double BP, double BS, int window) {
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
    if (windowType)
        for (int i = 0; i < N; i++)
            h[i] = h[i] * w[i];
}

void FLT_Filter::calc_h_fft_mag_ph_att()
{
    fftw_plan forward_h = fftw_plan_dft_r2c_1d(fft_size, h, h_fft, FFTW_ESTIMATE);
    fftw_execute(forward_h);

    for (int i = 0; i < fft_size / 2 + 1; i++) { // Pulse response magnitude, phase, attenuation calculation
        double magnitude = calc_magnitude(h_fft[i][0], h_fft[i][1]);
        double phase = calc_phase(h_fft[i][0], h_fft[i][1]);
        double attenuation = 20 * log10(1 / magnitude);

        h_magnitude[i] = magnitude;

        h_phase[i] = phase;
        h_attenuation[i] = attenuation;
    };

    fftw_destroy_plan(forward_h);
}

double FLT_Filter::calc_magnitude(double real_value, double complex_value)
{
    return sqrt(pow(real_value, 2) + pow(complex_value, 2));
}

double FLT_Filter::calc_phase(double real_value, double complex_value)
{
    return atan2(complex_value, real_value);
}

int FLT_Filter::convolFull(Frame& frame1, Frame& frame2)
{
    if (frame1.size != frame2.size) {
        printf_s("Convolution ERROR (convolFull) array sizes are different\n");
        return 0;
    }
    int size = frame1.size;
    for (int i = 0; i < size; i++)
        conv_frames[i] = frame1.data[i];
    for (int i = 0; i < add; i++)
        conv_frames[size + i] = frame1.data[size + i] + frame2.data[i];
    for (int i = 0; i < size; i++)
        conv_frames[size + add + i] = frame2.data[add + i];
    return 1;
}

void FLT_Filter::convolDifferent(Frame& frame1, Frame& frame2) {
    int pos = 0;

    int frame1_pre1 = frame1.add_other_left;
    int frame1_len1 = min_N_left;
    int frame1_pre2 = frame1_pre1 + frame1_len1;
    int frame1_len2 = abs(frame1.size - min_N_left);
    int frame1_pre3 = frame1_pre2 + frame1_len2;
    int len3 = min_N_left * 2;

    int frame2_pre3 = frame2.add_other_left;
    int frame2_pre4 = frame2_pre3 + len3;
    int frame2_len4 = abs(frame2.size - min_N_left);
    int frame2_pre5 = frame2_pre4 + frame2_len4;
    int frame2_len5 = min_N_left;

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

FLT_Filter::Frame::Frame()
{
}

FLT_Filter::Frame::~Frame()
{
    if (data != nullptr) {
        delete[] data;
    }
}

void FLT_Filter::Frame::setData(Frame& frame)
{
    this->setData(frame.data, 0, frame.size);
}

void FLT_Filter::Frame::init(int N, int fft_size)
{
    this->N = N;
    this->fft_size = fft_size;
}

template<typename T>
void FLT_Filter::Frame::setData(T* arr, int begin, int size)
{
    if (size != this->size) {
        if (data != nullptr) {
            delete[] data;
        }
        data = new double[size];
        this->size = size;
    }

    for (int i = 0; i < fft_size; i++) {
        if (i < size)  // Сигнал
            data[i] = static_cast<double&>(arr[begin + i]);
        else        // Дополнить нулями для БПФ
            data[i] = 0;
    }

    add_other = fft_size - size - (N - 1);
    add_other_left = add_other / 2;
    add_other_right = add_other - add_other_left;
}
