#pragma once
#include "FLT_BaseFilter.h"

#define FILTER_FIRST_PKT 123
#define FILTER_INPUT 125
#define FILTER_ERROR_FFT 12516

class FLT_FilterPkt :
    public FLT_BaseFilter
{
protected:

private:
    int packet_size = 0;
    int frame_size = 0;
    int packet_index = 0;
    // add_min2 = add_min / 2;
    int add_min2 = 0;

    double* ptrToAllocatedData1 = NULL;
    double* ptrToAllocatedData2 = NULL;
    double* packet1 = NULL;
    double* packet1_add_right = nullptr;

    double* packet2 = NULL;
    double* packet2_add_left = nullptr;
    double* packet2_add_right = nullptr;

    void filtrateFirstPacket(double* packet);
    void filtrateIntermediatePacket(double* packet);
    bool check_fft_size(int fft_size);
public:
    //bool setParameters(int N, int accurancy, int length);
    
    // -----------------------------------------------
    bool startTransfer(int length, int accurancy);
    // С хвостами, в тот же массив
    bool filtratePkt1(double* packet);
    int stopTransfer(double*& lastPacket);
    // -----------------------------------------------
    
    bool startTransferBlock(int packet_size, int fft_size);
    // Копирует данные в казатель (медленнее, чем filtratePktBlock2)
    bool filtratePktBlock1(double* const packet);
    // Возвращает указатель на данные (быстрее, чем filtratePktBlock1)
    double* const filtratePktBlock2(double* const packet);
    // Возвращает указатель на массив double
    double* getLastPktBlock1();
    void stopTransferBlock();

};

