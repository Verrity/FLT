#pragma once
#include "fftw3.h"
#include <math.h>
#include <vector>

#define FILTER_PI 3.1415926535897932384626433832795

#define FILTER_ERROR_OK          1
#define FILTER_ERROR_N          -1   // Значение N недопустимо
#define FILTER_ERROR_FD         -2   // Значение частоты дискретизации недопустимо
#define FILTER_ERROR_ACCURANCY  -3   // Значение точности недопустимо
#define FILTER_ERROR_WINDOW     -4   // Значение окна недопустимо
#define FILTER_ERROR_BAND       -5   // Одна из полос частот недопустима
#define FILTER_ERROR_LENGTH     -6 

class FLT_BaseFilter
{
public:
    FLT_BaseFilter();
    ~FLT_BaseFilter();
protected:
    struct Frame
    {
    public:
        Frame();
        ~Frame();
        double* data = nullptr;
        fftw_complex* data_fft = nullptr;
        int data_size = 0;
        int add = 0;
        int add_min = 0;
        int add_min_left = 0;
        int add_other = 0;
        int add_other_left = 0;
        int add_other_right = 0;

        void init(int N, int data_size, int fft_size);
        void setData(double* arr, int begin, int size);
        void setData(Frame& frame);
        static void switchData(Frame& toFrame, Frame& fromFrame);
    private:
        int fft_size = 0;
        int N = 0;
    };

    int type = 0;
    int method = 0;
    int bands_count = 0;

    int N = 0;
    double fd = 0;
    int fft_size = 0;
    int accurancy = 0;
    int add_min = 0;
    int error_code = FILTER_ERROR_OK;
    std::vector<double> bands;
    int window = 0;
    int conv_size = 0;

    double* h = nullptr;                    // Real: fft_size                   Logical: fft_size
    fftw_complex* h_fft = nullptr;          // Real: fft_size/2 + 1             Logical: fft_size/2 + 1
    double* h_magnitude = nullptr;          // Real: fft_size                   Logical: fft_size
    double* h_phase = nullptr;              // Real: fft_size                   Logical: fft_size
    double* h_attenuation = nullptr;        // Real: fft_size                   Logical: fft_size
    double* freq_match = nullptr;           // Real: fft_size                   Logical: fft_size
    double* w = nullptr;                    // Real: fft_size                   Logical: fft_size
    double* conv_frames = nullptr;          // Real: fft_size               Logical: fft_size
    fftw_complex* mul_frames_fft = nullptr; // Real: (fft_size/2 + 1)*2     Logical: fft_size/2 + 1

    fftw_plan forward_signal_fft;   // План БПФ
    fftw_plan backward_signalF;     // План ОБПФ
    double* pFrameData = nullptr;
    fftw_complex* pFrameDataFFT = nullptr;

    Frame frame1, frame2;

    // ---------------- CHECKER
    // --- PARAMETERS
    bool check_N(int N);
    bool check_fd(double fd);
    bool check_accurancy(int accurancy);
    bool check_window(int window);
    // --- BANDS
    bool check_band_LowpassR1B1(double band, double fd);

    bool convolFull(Frame& frame1, Frame& frame2);
    void convolDifferent(Frame& frame1, Frame& frame2);

    // ---------------- CALCULATION
    // --- IR
    int calc_IrLowpassR1B1();
    // --- SUB
    void calc_h();
    void calc_h_fft();
    void calc_h_fft_mag_ph_att();
    void calc_window();
    double calc_magnitude(double& real_value, double& complex_value);
    double calc_phase(double& real_value, double& complex_value);
    //bool init(int N, double fd, int accurancy, int window);
    bool fft_filtrate(Frame& frame);
private:
    int filtrate_length = 0;
    int filtrate_accurancy = 0;

public:
    virtual int filtrate(double* in, int length, double*& out, int accurancy, bool tails);
    // ---------------- IMPULSE RESPONSE
    bool setIrLowpassR1B1(int N, double fd, double band, int window);

    // ---------------- GET
    int get_fft_size();
    int get_N();
    int get_window();
    int get_h(double* &h);
    int get_h_magnitude(double* &magnitude, bool symmetrical);
    int get_h_phase(double* &phase, bool symmetrical);
    int get_h_attenuation(double* &attenuation, bool symmetrical);
    int get_freq_match(double* &freq_match, bool symmetrical);
    int get_w(double* &w);
    int get_error_code();
};
