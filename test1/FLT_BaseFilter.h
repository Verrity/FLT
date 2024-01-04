#pragma once
#include "fftw3.h"
#include <math.h>
#include <vector>
#include <cstring>
#include <iostream>
#include <chrono>

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
private:
    FLT_BaseFilter(const FLT_BaseFilter& other) = delete;
protected:
    struct Frame
    {
    public:
        Frame();
        ~Frame();
        double* data = nullptr;
        fftw_complex* data_fft = nullptr;
        int data_size = 0;

        void init(int N, int data_size, int fft_size);
        void setData(double* arr, int begin, int size);
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
    int signal_size = 0;
    int add_min = 0;
    int add_min2 = 0;
    int error_code = FILTER_ERROR_OK;
    std::vector<double> bands;
    int window = 0;
    int conv_size = 0;

    double* h = nullptr;                    // Real: fft_size                   Logical: fft_size
    fftw_complex* h_fft = nullptr;          // Real: fft_size/2 + 1             Logical: fft_size/2 + 1
    double* w = nullptr;                    // Real: fft_size                   Logical: fft_size
    fftw_complex* mul_frames_fft = nullptr; // Real: (fft_size/2 + 1)*2     Logical: fft_size/2 + 1
    
    fftw_plan forward_signal_fft;   // План БПФ
    fftw_plan backward_signalF;     // План ОБПФ
    double* pFrameData = nullptr;
    fftw_complex* pFrameDataFFT = nullptr;

    Frame frame1, frame2;


    using clock_t = std::chrono::high_resolution_clock; // Для таймера
    std::chrono::time_point<clock_t> m_startTime;   // Время начала
    std::chrono::time_point<clock_t> m_endTime;     // Время конца
    std::chrono::nanoseconds m_durationTime_ns;     // Для времени
    unsigned long int m_int_durationTime_ns = 0;    // Время [нс]
    void startTimer();
    void stopTimer();

    // ---------------- CHECKER
    // --- PARAMETERS
    bool check_N(int N);
    bool check_fd(double fd);
    bool check_accurancy(int accurancy);
    bool check_window(int window);
    // --- BANDS
    bool check_bands(double band, double fd);

    // ---------------- CALCULATION
    // --- IR
    int calc_IrLowpassR1B1();

    // --- SUB
    void calc_h();
    void calc_h_fft();
    void calc_window();
    double calc_magnitude(double& real_value, double& complex_value);
    double calc_phase(double& real_value, double complex_value);

    // ---------------- Filtration and initialisation
    bool fft_filtrate(Frame& frame);
    /* Инициализировать и рассчитать минимально необходимые для
    работы фильтра переменные и массивы
    (h, h_fft, mul_frames_ftt, fftw plans, frame1, frame2, w)
    */
    void init_min(int N, int fd, int fft_size, int frame_size);
    void free_min();
    bool isMinAllocated = false;

public:
    /* Simple filtering using a single fft transformation use this for small signals, preferably equal to length = fft - N + 1
    * initializes arrays and calculates the impulse response
    * signal - input & output data
    * length - signal length
    * accuracy - an increase in the length of the FFT by 2^accuracy relative to the minimum possible FFT
    - put nullptr to signal if you don't want to filtrate, just use get_... functions with arrays
    - put signal if you want to filtrate and use get_... functions with arrays
    returns false if error
    */
    bool filtrate(double* const signal, int length, int accurancy);

    /* Simple filtering using a single fft transformation use this for small signals, preferably equal to length = fft - N + 1
    * initializes arrays and calculates the impulse response
    * signal - input data
    * length - signal length
    * output - put empty pointer for output signal with impulse response (tails)
    * accuracy - an increase in the length of the FFT by 2^accuracy relative to the minimum possible FFT
    - put nullptr to signal if you don't want to filtrate, just use get_... functions with arrays
    - put signal if you want to filtrate and use get_... functions with arrays
    returns output signal length, 0 if error 
    */
    int filtrateT(double* const signal, int length, double* &output, int accurancy);

    // ---------------- IMPULSE RESPONSE

    /*Lowpass, Method: weighting, bands: 1
    * N - impulse response length, must be odd, min value 17
    * fd - sampling rate in Hz
    * band - band in Hz
    * window - window type (0 - none, 1 - Hamming, 2 - Blackman, 3 - Hann)
    */
    bool setIrLowpassR1B1(int N, double fd, double band, int window);

    /*Measures attenuation from f_low [Hz] to f_high [Hz] in increments of step in Hz
    output - put empty pointer for filter attenuation [dBu]
    * To accurately measure attenuation at low frequencies, reduce the sampling rate.
    * To accurately measure attenuation at higher frequencies, increase the sampling rate.
    * output - put empty pointer for measured attenuation
    * length - the length of your input signal
    * accurancy - an indicator of how many times the involved size
    of the fft will differ from the minimum size of the fft
    make it bigger for more attenuation, reduces the speed, the minimum value is 0
    returns output length, 0 if error, check error_code
    */
    virtual int measureAttenuation(double*& output, int length, int accurancy, double f_low, double step, double f_high, unsigned long &time_ns);

    // ---------------- GET
    int get_fft_size();
    int get_N();
    int get_window();

    /*h - put empty pointer for output
    returns output length, 0 if error 
    before that, the memory must be initialized, see filtrate functions
    */
    int get_h(double* &h);

    /*Get impulse response 
    magnitude - put empty pointer for output
    returns output length, 0 if error
    use symmetrical  true - if you want to get from 0 Hz to Fd
    false - if you want to get from 0 Hz to Fn
    before that, the memory must be initialized, see filtrate functions
    */
    int get_h_magnitude(double* &magnitude, bool symmetrical);

    /*Get phase of the approximating filter function
    phase - put empty pointer for output
    returns output length, 0 if error
    use symmetrical  true - if you want to get from 0 Hz to Fd
    false - if you want to get from 0 Hz to Fn
    before that, the memory must be initialized, see filtrate functions
    */
    int get_h_phase(double* &phase, bool symmetrical);

    /*Get attenuation of the approximating filter function
    attenuation - put empty pointer for output
    returns output length, 0 if error
    use symmetrical  true - if you want to get from 0 Hz to Fd
    false - if you want to get from 0 Hz to Fn
    before that, the memory must be initialized, see filtrate functions
    */
    int get_h_attenuation(double* &attenuation, bool symmetrical);

    /*Get array of frequencies corresponding to the harmonics of the 
    parameters of the approximating filter function
    freq_match - put empty pointer for output
    returns output length, 0 if error
    use symmetrical  true - if you want to get from 0 Hz to Fd
    false - if you want to get from 0 Hz to Fn
    before that, the memory must be initialized, see filtrate functions
    */
    int get_freq_match(double* &freq_match, bool symmetrical);

    /*Get widow
    freq_match - put empty pointer for output
    returns output length, 0 if error
    before that, the memory must be initialized, see filtrate functions
    */
    int get_w(double* &w);

    /*Returns error code*/
    int get_error_code();
};
