#pragma once
#include "FLT_FilterFile.h"

#define FILTER_FIRST_PKT 123
#define FILTER_ERROR_VALUE 125
#define FILTER_ERROR_FFT 12516
#define FILTER_ERROR_FUNCTION 1597

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
public:
    // -----------------------------------------------
    // use this if you want to init arrays for get_magnitude or other
    bool startTransfer(int length, int accurancy);
    // С хвостами, в тот же массив
    bool filtratePkt1(double* packet);
    double* getLatestPkt();
    bool stopTransfer();
    // -----------------------------------------------
    
    /*Allocates memory and means ready for use filtratePktBlock*
    packet_size - size of your packet (const)
    use this if you want to init arrays for get_magnitude or other
    */
    bool startTransferBlock(int packet_size, int accurancy);
    /*Puts data in your pointer (slower than filtratePktBlock2)
    packet - your input and output data
    returns false, if packet is first (you must skip iteration)
    returns true, if the packet is intermediate (returns the previous package at the current iteration)
    [use getLastPktBlock1() to get the latest packet]
    */
    bool filtratePktBlock(double* const packet);
    // Возвращает указатель на данные (быстрее, чем filtratePktBlock1)
    /*Puts data in your pointer
    returns a pointer to an array of data
    you free up the memory yourself
    */
    // Возвращает указатель на массив double
    double* getLatestPktBlock();

    bool stopTransferBlock();

};

