#pragma once
#include <unordered_map>
#include "fftw3.h"
typedef unsigned long FILTER;

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

    double BP1 = 0, BP2 = 0, BP3 = 0, BP4 = 0;

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

    Frame frame1, frame2;


    int convolFull(Frame& frame1, Frame& frame2);
    void convolDifferent(Frame& frame1, Frame& frame2);

    FILTER getId();
    static FLT_Filter* getInstanse(FILTER filter);
private:
    FILTER descriptor = 0;
    static FILTER nextDescriptor;
    //using objectPtr = std::shared_ptr<FLT_Filter>;
    static std::unordered_map<FILTER, FLT_Filter*> objects;
};
