#pragma once
#include "FLT_BaseFilter.h"

class FLT_FilterFile : public FLT_BaseFilter {
protected:
    int frame_size = 0;
    double* left_tail = nullptr;
    double* right_tail = nullptr;
    void _filtrateBlock(double* const signal, int length, double* left_tail, double* right_tail);
    //int length = 0;
public:
    FLT_FilterFile();
    ~FLT_FilterFile();

    /* ‘ильтрует массив поблочно, без хвостов
    * возвращает результат выполнени€
    * accurancy - начальное и рекомендуемое значение 0
    * put nullptr to signal if you want to init arrays for get_magnitude or other
    */
    bool filtrateBlock(double* const signal, int length, int accurancy);
    /* ‘ильтрует массив поблочно, с хвостами
    * возвращает указытель на массив размером length + N - 1
    * accurancy - начальное и рекомендуемое значение 0
    * put nullptr to signal if you want to init arrays for get_magnitude or other
    */
    double* filtrateBlockT(double* const signal, int length, int accurancy);
    // ---------------- GET
    int get_frame_size();
private:
    int filtrateBlock_length = 0;
    int filtrateBlock_accurancy = 0;
    //bool local_init(int N, double fd, int accurancy, int window);
    //bool fft_filtrate(Frame& frame);

};
