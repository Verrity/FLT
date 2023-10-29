#pragma once
#include "FLT_BaseFilter.h"

#define FILTER_FIRST_PKT 123
#define FILTER_INPUT 125

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

    double* packet1 = nullptr;
    double* packet1_add_right = nullptr;

    double* packet2 = nullptr;
    double* packet2_add_left = nullptr;
    double* packet2_add_right = nullptr;

    int filtrateFirstPacket(double* packet);
public:
    //bool setParameters(int N, int accurancy, int length);
    
    // -----------------------------------------------
    bool startTransfer(int length, int accurancy);
    // С хвостами, в тот же массив
    bool filtratePkt1(double* packet);
    int stopTransfer(double*& lastPacket);
    // -----------------------------------------------
    
    bool startTransfer10(int length, int fft_size);
    bool filtratePkt10(double* packet);
    int stopTransfer10(double* packet);

};

