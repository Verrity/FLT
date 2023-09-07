#pragma once
#include <unordered_map>
#include "fftw3.h"
typedef unsigned long FILTER;

#define FILTER_PI 3.1415926535897932384626433832795

#define FLT_OK              0
#define FLT_WRONG_PARAMETER 1
#define FLT_UNKNOWN_ERROR   2
#define FILTER_ERROR_N      3   // N value is wrong
#define FILTER_ERROR_LENGTH 4   // Sgnal length value is wrong
#define FILTER_ERROR_BAND   5   // One of frequency band is wrong
#define FILTER_ERROR_FD     6   // Fd is wrong
#define FILTER_ERROR_Window 7   // Window value is wrong

#define FILTER_ERROR_PREV   100

class FLT_Filter
{
public:

    FLT_Filter(FILTER& filter, int N, double fd);
    ~FLT_Filter();

    int N = 0;
    double fd = 0;
    int fft_size = 0;
    int add = 0;
    int min_N = 0;
    int min_N_left = 0;

    double B1 = 0, B2 = 0, B3 = 0, B4 = 0;

    double* h;                      // Real: fft_size               Logical: fft_size
    fftw_complex* h_fft;            // Real: fft_size/2 + 1         Logical: fft_size/2 + 1
    double* h_magnitude;            // Real: fft_size               Logical: fft_size
    double* h_phase;                // Real: fft_size               Logical: fft_size
    double* h_attenuation;          // Real: fft_size               Logical: fft_size
    double* freq_match;             // Real: fft_size               Logical: fft_size
    double* w;                      // Real: fft_size               Logical: fft_size

    double* prev_frame;             // Real: fft_size               Logical: fft_size
    double* curr_frame;             // Real: fft_size               Logical: fft_size
    double* conv_frames;            // Real: fft_size               Logical: fft_size
    fftw_complex* prev_frame_fft;   // Real: fft_size/2 + 1         Logical: fft_size/2 + 1
    fftw_complex* curr_frame_fft;   // Real: fft_size/2 + 1         Logical: fft_size/2 + 1
    fftw_complex* mul_frames_fft;   // Real: (fft_size/2 + 1)*2     Logical: fft_size/2 + 1

    fftw_plan forward_signal_fft;   // œÎ‡Ì ¡œ‘
    fftw_plan backward_signalF;     // œÎ‡Ì Œ¡œ‘

    struct Frame
    {
        Frame();
        ~Frame();

        int size = 0;
        double* data = nullptr;

        int add_other;
        int add_other_left;
        int add_other_right;

        template <typename T>
        void setData(T* arr, int begin, int size);
        void setData(Frame& frame);
        void init(int N, int fft_size);
    private:
        int fft_size = 0;
        int N = 0;
    };

    int windowType = 1;
    Frame frame1, frame2;

    int check_params(const char* filter, int N, int fd, int BP1, int BP2, int BS1, int BS2, int window);
    void createIR(const char* filter);
    void calc_h_fft_mag_ph_att();
    double calc_magnitude(double real_value, double complex_value);
    double calc_phase(double real_value, double complex_value);

    // —‚∏ÚÍ‡
    int convolFull(Frame& frame1, Frame& frame2);
    void convolDifferent(Frame& frame1, Frame& frame2);

    FILTER get_id();
    FLT_Filter& get_ref();
    static FLT_Filter* get_pointer(FILTER id);
private:
    FILTER descriptor = 0;
    int errorCode = 0;
    static FILTER next_descriptor;
public:
    static std::unordered_map<FILTER, FLT_Filter*> list;
};
