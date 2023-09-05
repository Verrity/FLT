#include "FLT_filter.h"

FILTER FLT_Filter::nextDescriptor = 1;
std::unordered_map<FILTER, FLT_Filter*> FLT_Filter::objects;

FLT_Filter::FLT_Filter(FILTER& filter, int N, double fd) : N(N), fd(fd)
{
    descriptor = nextDescriptor;
    nextDescriptor++;
    objects[descriptor] = this;
}

FLT_Filter::~FLT_Filter()
{

}

FILTER FLT_Filter::getId() {
    return descriptor;
}

FLT_Filter* FLT_Filter::getInstanse(FILTER filter) {
    auto it = objects.find(filter);
    if (it != objects.end()) {
        FLT_Filter* ptr = it->second;
        return ptr;
    }
    return nullptr;
}

int FLT_Filter::convolFull(Frame& frame1, Frame& frame2)
{
    if (frame1.size != frame2.size) {
        printf_s("Convolution ERROR (convolFull) array sizes are different\n");
        return 0;
    }
    int size = frame1.size;
    for (int i = 0; i < size; i++)
        conv_frames[i] = frame1.data[i];
    for (int i = 0; i < add; i++)
        conv_frames[size + i] = frame1.data[size + i] + frame2.data[i];
    for (int i = 0; i < size; i++)
        conv_frames[size + add + i] = frame2.data[add + i];
    return 1;
}

void FLT_Filter::convolDifferent(Frame& frame1, Frame& frame2) {
    int pos = 0;

    int frame1_pre1 = frame1.add_other_left;
    int frame1_len1 = min_N_left;
    int frame1_pre2 = frame1_pre1 + frame1_len1;
    int frame1_len2 = abs(frame1.size - min_N_left);
    int frame1_pre3 = frame1_pre2 + frame1_len2;
    int len3 = min_N_left * 2;

    int frame2_pre3 = frame2.add_other_left;
    int frame2_pre4 = frame2_pre3 + len3;
    int frame2_len4 = abs(frame2.size - min_N_left);
    int frame2_pre5 = frame2_pre4 + frame2_len4;
    int frame2_len5 = min_N_left;

    // -- 1 --
    for (int i = 0; i < frame1_len1; i++)
        conv_frames[i] = frame1.data[frame1_pre1 + i];
    pos = frame1_len1;

    // -- 2 --
    for (int i = 0; i < frame1_len2; i++)
        conv_frames[pos + i] = frame1.data[frame1_pre2 + i];
    pos += frame1_len2;

    // -- 3 --
    for (int i = 0; i < len3; i++)
        conv_frames[pos + i] = frame1.data[frame1_pre3 + i] + frame2.data[frame2_pre3 + i];
    pos += len3;

    // -- 4 --
    for (int i = 0; i < frame2_len4; i++)
        conv_frames[pos + i] = frame2.data[frame2_pre4 + i];
    pos += frame2_len4;

    // -- 5 --
    for (int i = 0; i < frame2_len5; i++)
        conv_frames[pos + i] = frame2.data[frame2_pre5 + i];
}

FLT_Filter::Frame::Frame()
{
}

FLT_Filter::Frame::~Frame()
{
    if (data != nullptr) {
        delete[] data;
    }
}

void FLT_Filter::Frame::setData(Frame& frame)
{
    this->setData(frame.data, 0, frame.size);
}

void FLT_Filter::Frame::init(int N, int fft_size)
{
    this->N = N;
    this->fft_size = fft_size;
}

template<typename T>
void FLT_Filter::Frame::setData(T* arr, int begin, int size)
{
    if (size != this->size) {
        if (data != nullptr) {
            delete[] data;
        }
        data = new double[size];
        this->size = size;
    }

    for (int i = 0; i < fft_size; i++) {
        if (i < size)  // Сигнал
            data[i] = static_cast<double&>(arr[begin + i]);
        else        // Дополнить нулями для БПФ
            data[i] = 0;
    }

    add_other = fft_size - size - (N - 1);
    add_other_left = add_other / 2;
    add_other_right = add_other - add_other_left;
}
