#pragma once
#include "FLT_BaseFilter.h"

#define FILTER_ERROR_LENGTH_PARAMETER -7;

class FLT_FilterFile : public FLT_BaseFilter {
protected:
    int frame_size = 0;
    double* left_tail = nullptr;
    double* right_tail = nullptr;
    bool check_length_parameter(int length_parameter);
    void _filtrateBlock(double* const signal, int length, double* left_tail, double* right_tail);
public:
    FLT_FilterFile();
    ~FLT_FilterFile();

    /* Splits the signal into blocks and filters. It is used on signals of any magnitude, it is effective in processing large sequences
    * initializes arrays and calculates the impulse response
    * signal - input & output data
    * length_parameter - increasing the length of the processing frame and hence the FFT by 2^length_parameter
    * accuracy - an increase in the length of the FFT by 2^accuracy relative to the minimum possible FFT
    - put nullptr to signal if you don't want to filtrate, just use get_... functions with arrays
    - put signal if you want to filtrate and use get_... functions with arrays
    returns false if error
    */
    bool filtrateBlock(double* const signal, int length, int length_parameter, int accurancy);
    
    /* Splits the signal into blocks and filters. It is used on signals of any magnitude, it is effective in processing large sequences, returns signal with impulse response (tails)
    * function use block filtrating if length >= N * 45
    * initializes arrays and calculates the impulse response
    * signal - input data
    * output - put empty pointer for output signal with impulse response (tails)
    * length_parameter - increasing the length of the processing frame and hence the FFT by 2^length_parameter
    * accuracy - an increase in the length of the FFT by 2^accuracy relative to the minimum possible FFT
    - put nullptr to signal if you don't want to filtrate, just use get_... functions with arrays
    - put signal if you want to filtrate and use get_... functions with arrays
    returns output signal length, 0 if error
    */
    int filtrateBlockT(double* const signal, double* &output, int length, int length_parameter, int accurancy);

    /*Measures attenuation from f_low [Hz] to f_high [Hz] in increments of step in Hz
    output - put empty pointer for filter attenuation [dBu]
    * To accurately measure attenuation at low frequencies, reduce the sampling rate.
    * To accurately measure attenuation at higher frequencies, increase the sampling rate.
    * output - put empty pointer for measured attenuation
    * length - the length of your input signal
    * accuracy - an increase in the length of the FFT by 2^accuracy relative to the minimum possible FFT
    of the fft will differ from the minimum size of the fft
    make it bigger for more attenuation, reduces the speed, the minimum value is 0
    returns output length, 0 if error, check error_code
    */
    virtual int measureAttenuation(double*& output, int length, int accurancy, double f_low, double step, double f_high, unsigned long& time_ns) override;
    
    // ---------------- GET
    int get_frame_size();
private:
    int filtrateBlock_length = 0;
    int filtrateBlock_accurancy = 0;

    //bool local_init(int N, double fd, int accurancy, int window);
    //bool fft_filtrate(Frame& frame);

};
