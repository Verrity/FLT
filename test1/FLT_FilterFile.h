#pragma once
#include "FLT_BaseFilter.h"

class FLT_FilterFile : public FLT_BaseFilter {
protected:
    int sample_size = 0;
    int length = 0;
public:
    FLT_FilterFile();
    ~FLT_FilterFile();

    // INPUT (N, fd, accurancy, window)
    // return error_code
    bool setIrLowpassR1B1(int N, double fd, int accurancy, double band, int window);

    int filtrate(double* in, int length, double* out, bool tails);
    // Returns output array size
    int filtrateBlock(double* in, int length, double* out, bool tails);

    // ---------------- GET
    int get_sample_size();
private:
    bool local_init(int N, double fd, int accurancy, int window);
    bool fft_filtrate(Frame& frame);
};
