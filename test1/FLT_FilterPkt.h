#pragma once
#include "FLT_BaseFilter.h"

#define FILTER_FIRST_PKT 123
#define FILTER_ERROR_VALUE 125
#define FILTER_ERROR_FFT 12516

class FLT_FilterPkt :
    public FLT_BaseFilter
{
protected:
    using FLT_BaseFilter::filtrate;
private:
    int value = 1;
    int min_fft = 0;

    int packet_size = 0;
    int frame_size = 0;
    int packet_index = 0;

    double* ptrToAllocatedData1 = nullptr;
    double* ptrToAllocatedData2 = nullptr;
    double* packet1 = nullptr;
    double* right_tail0 = nullptr;

    double* packet2 = nullptr;
    double* left_tail = nullptr;
    double* right_tail = nullptr;

    void filtratePacketBlock(double* packet, double* left_tail, double* right_tail);
public:

    // -----------------------------------------------
    bool startTransfer(int length, int accurancy);
    // С хвостами, в тот же массив
    bool filtratePkt1(double* packet);
    int stopTransfer(double*& lastPacket);
    // -----------------------------------------------
    
    /*Allocates memory and means ready for use filtratePktBlock*
    packet_size - size of your packet (const)
    fft_size 
    */
    bool startTransferBlock(int packet_size, unsigned int value);
    /*Puts data in your pointer (slower than filtratePktBlock2)
    packet - your input and output data
    returns false, if packet is first (you must skip iteration)
    returns true, if the packet is intermediate (returns the previous package at the current iteration)
    [use getLastPktBlock1() to get the latest packet]
    */
    bool filtratePktBlock1(double* const packet);
    // Возвращает указатель на данные (быстрее, чем filtratePktBlock1)
    /*Puts data in your pointer
    returns a pointer to an array of data
    you free up the memory yourself
    */
    double* const filtratePktBlock2(double* const packet);
    // Возвращает указатель на массив double
    double* getLatestPktBlock1();
    void stopTransferBlock();

};

