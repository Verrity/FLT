#pragma once
#include "FLT_FilterFile.h"

#define FILTER_FIRST_PKT -8
#define FILTER_ERROR_FUNCTION -9

class FLT_FilterPkt :
    public FLT_FilterFile
{
protected:
    using FLT_BaseFilter::filtrate;
private:
    //int value = 1;
    int min_fft = 0;

    int packet_size = 0;
    int packet_index = 0;

    double* ptrToAllocatedData1 = nullptr;
    double* ptrToAllocatedData2 = nullptr;
    double* packet1 = nullptr;
    double* packet2 = nullptr;
    double* right_tail0 = nullptr;

    using FLT_FilterFile::filtrateBlock;
    using FLT_FilterFile::filtrateBlockT;
    using FLT_FilterFile::measureAttenuation;
public:
    /*Don't splits the signal into blocks. It is used on signals of any size, it is effective in processing small sequences
    * initializes arrays and calculates the impulse response and starts simple filtering
    * packet_size - const size of your packet
    * accuracy - an increase in the length of the FFT by 2^accuracy relative to the minimum possible FFT
    * use filtratePkt to filtrate packets
    */
    bool startTransfer(int packet_size, int accurancy);

    /* First packet will be loaded, function return false and error code
    FILTER_FIRST_PKT, the iteration must be skipped, 
    * when the next package is loaded, the previous package
    will be returned, function return true
    * packet - input & output data
    use getLatestPkt() to get the latest packet
    */
    bool filtratePkt(double* packet);

    /*Returns pointer to latest packet size of packet_size
    * use stopTransferBlock to stop stansfer
    */
    double* getLatestPkt();

    /* Stops transfer */
    bool stopTransfer();

    // -------------------------------------------------------------------------
    
    /*Splits the signal into blocks and filters. It is used on signals of any size, it is effective in processing large sequences
    * initializes arrays and calculates the impulse response and starts block filtering
    * packet_size - const size of your packet
    * length_parameter - increasing the length of the processing frame and hence the FFT by 2^length_parameter
    * accuracy - an increase in the length of the FFT by 2^accuracy relative to the minimum possible FFT
    * the length of the packet should fit 4 blocks, will return 0
    and the code FILTER_ERROR_VALUE if this condition is not met, increase accuracy or packet_size.
    * use filtratePktBlock to filtrate packets
    */
    bool startTransferBlock(int packet_size, int length_parameter, int accurancy);

    /* First packet will be loaded, function return false and error code
    FILTER_FIRST_PKT, the iteration must be skipped
    * when the next package is loaded, the previous package
    will be returned, function return true
    * packet - input & output data
    use getLastPktBlock1() to get the latest packet
    */
    bool filtratePktBlock(double* const packet);

    /*Returns pointer to latest packet size of packet_size
    * use stopTransferBlock to stop stansfer
    */
    double* getLatestPktBlock();

    /* Stops transfer */
    bool stopTransferBlock();

    /*Measures attenuation from f_low [Hz] to f_high [Hz] in increments of step in Hz
    output - put empty pointer for filter attenuation [dBu]
    * To accurately measure attenuation at low frequencies, reduce the sampling rate.
    * To accurately measure attenuation at higher frequencies, increase the sampling rate.
    * output - put empty pointer for measured attenuation
    * length - the length of your input signal
    * accuracy - an increase in the length of the FFT by 2^accuracy relative to the minimum possible FFT
    returns output length, 0 if error, check error_code
    */
    int measureAttenuation(double*& output, int packet_size, int accurancy, double f_low, double step, double f_high, unsigned long &time_ns) override;

};

