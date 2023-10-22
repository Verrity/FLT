#pragma once
#include "FLT_BaseFilter.h"

#define FILTER_FIRST_PKT 123
#define FILTER_INPUT 125

class FLT_FilterPkt :
    public FLT_BaseFilter
{
protected:
    int packet_size = 0;
private:
    int packet_counter = 0;
public:
    //bool setParameters(int N, int accurancy, int length);
    bool startTransfer(int length, int accurancy);
    // С хвостами, в тот же массив
    bool filtratePkt1(double* packet);

    int stopTransfer(double*& lastPacket);
};

