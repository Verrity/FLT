#include "FLT_FilterFile.h"
#include <iostream>

FLT_FilterFile::FLT_FilterFile()
{
}

FLT_FilterFile::~FLT_FilterFile()
{
    if (left_tail != nullptr)   delete[] left_tail;
    if (right_tail != nullptr)  delete[] right_tail;
}

//int FLT_FilterFile::filtrate(double* in, int length, double* &out, bool tails) {
//    if (length <= 0) {
//        error_code = FILTER_ERROR_LENGTH;
//        return 0;
//    }
//    
//    frame1.setData(in, 0, length);
//    fft_filtrate(frame1);
//
//    std::cout << "frame1.add " << frame1.add << std::endl;
//    std::cout << "frame1.add_min " << frame1.add_min << std::endl;
//    std::cout << "frame1.add_min_left " << frame1.add_min_left << std::endl;
//    std::cout << "frame1.add_other " << frame1.add_other << std::endl;
//    std::cout << "frame1.add_other_left " << frame1.add_other_left << std::endl;
//    std::cout << "frame1.add_other_right " << frame1.add_other_right << std::endl;
//    std::cout << "frame1.data_size " << frame1.data_size << std::endl;
//    std::cout << "fft_size " << fft_size << std::endl;
//    
//    int len = 0;
//    if (tails) {
//        len = length + add_min;
//        out = new double[len];
//        for (int i = 0; i < len; i++)
//        {
//            out[i] = frame1.data[frame1.add_other_left + i];
//            //printf("%f", frame1.data[frame1.add_other_left + i])
//        }
//    }
//    else {
//        len = length;
//        out = new double[len];
//        for (int i = 0; i < len; i++)
//        {
//            out[i] = frame1.data[frame1.add_other_left + frame1.add_min_left + i];
//        }
//    }
//    return len;
//}

bool FLT_FilterFile::filtrateBlock(double* const signal, int length, int accurancy)
{
    if (    // Если массивы под эти параметры еще не рассчитывались, рассчитать
        (signal_size != length) ||
        (this->accurancy != accurancy)
        )
    {
        signal_size = length;
        frame_size = N * 45;
        fft_size = 1;
        while (fft_size < (frame_size + add_min)) // fft_size < (длина сигнала + минимум добавочных элементов)
        {
            fft_size = fft_size << 1;
        }
        fft_size << accurancy - 1;
        frame_size = fft_size - N + 1;

        // Если в длину сигнала укладывается 3 или меньше кадров
        // Использовать обычную фильтрацию, иначе блочную
        if (length <= N * 45 * 3) {
            fft_size = 1;
            while (fft_size < (length + add_min)) // fft_size < (длина сигнала + минимум добавочных элементов)
            {
                fft_size = fft_size << 1;
            }
            fft_size << accurancy;
            if (isMinAllocated)
                free_min();
            init_min(N, fd, fft_size, length);
        }
        else {
            if (isMinAllocated)
                free_min();
            init_min(N, fd, fft_size, frame_size);
        }
    }

    // Если в длину сигнала укладывается 3 или меньше кадров
    // Использовать обычную фильтрацию, иначе блочную
    // 3 - потому что frame_size берётся 3 раза, потом берутся оставшиеся элементы
    if (length <= N * 45 * 3) {
        filtrate(signal, length, accurancy);
    }
    else {
        _filtrateBlock(signal, length, nullptr, nullptr);
    }
    return true;
}

double* FLT_FilterFile::filtrateBlockT(double* const signal, int length, int accurancy)
{
    if (    // Если массивы под эти параметры еще не рассчитывались, рассчитать
        (signal_size != length) ||
        (this->accurancy != accurancy)
        )
    {
        signal_size = length;
        this->accurancy = accurancy;
        frame_size = N * 45;
        fft_size = 1;
        while (fft_size < (frame_size + add_min)) // fft_size < (длина сигнала + минимум добавочных элементов)
        {
            fft_size = fft_size << 1;
        }
        fft_size << accurancy - 1;
        frame_size = fft_size - N + 1;

        left_tail = new double[add_min2];
        right_tail = new double[add_min2];
    }

    // Если в длину сигнала укладывается 3 или меньше кадров
    // Использовать обычную фильтрацию, иначе блочную
    // 3 - потому что frame_size берётся 3 раза, потом берутся оставшиеся элементы
    if (length <= frame_size * 3) {
        fft_size = 1;
        while (fft_size < (length + add_min)) // fft_size < (длина сигнала + минимум добавочных элементов)
            fft_size = fft_size << 1;
        fft_size << accurancy;

        if (isMinAllocated)
            free_min();
        init_min(N, fd, fft_size, length);

        return filtrateT(signal, length, accurancy);
    }
    else {
        if (isMinAllocated)
            free_min();
        init_min(N, fd, fft_size, frame_size);

        double* out = new double[length + add_min];
        double* ptrToOutBeginWriteSignal = out + add_min2;
        memcpy(ptrToOutBeginWriteSignal, signal, sizeof(double) * length);

        _filtrateBlock(ptrToOutBeginWriteSignal, length, left_tail, right_tail);

        for (int i = 0; i < add_min2; i++)
            out[i] = left_tail[i];
        for (int i = 0; i < add_min2; i++)
            out[add_min2 + length + i] = right_tail[i];
        return out;
    }
}



void FLT_FilterFile::_filtrateBlock(double* const signal, int length, double* left_tail, double* right_tail) {
    int frames_count = length % frame_size != 0 ? floor(double(length) / frame_size) + 1 : floor(double(length) / frame_size);
    int begFrame = 0;	// Текущая позиция в frame1
    int begSignal = 0;	// Текущая позиция в записывании Signal

    for (int j = 0; j < frames_count; j++) {
        // ------------------ Свёртка и фильтрация кадров
        switch (j)
        {
        case 0:	// Если кадр первый
            // Скопировать данные в frame1.data и фильтровать
            frame1.setData(signal, 0, frame_size);
            fft_filtrate(frame1);
            /*
                Свёртка frame1 и frame2, и запись кадра frame1 без левого хвоста
                в Signal, frame2 выйдет со свёрнутым левым хвостом,
                далее его логическое начало должно быть сдвинуто на add_min2.
            */
            if (left_tail != nullptr)
            {
                // Записываем левый хвост в left_tail
                for (int i = 0; i < add_min2; i++)
                    left_tail[i] = frame1.data[i];
            }
            begFrame += add_min2;

            // Переписываем от конца левого хвоста (обрезали) до участка суммирования
            for (int i = 0; i < frame_size - add_min2; i++)
                signal[i] = frame1.data[begFrame + i];
            begFrame += frame_size - add_min2;
            begSignal += frame_size - add_min2;
            break;

        case 1:	// Если кадр второй, то frame2 требуется загрузить и закончить свёртку с frame1 (i == 1) и frame2 (i == 2)
            // Скопировать данные в frame2.data и фильтровать
            frame2.setData(signal, frame_size, frame_size);
            fft_filtrate(frame2);
            // Получили frame2, можно закончить свёртку frame1 (i == 1) и frame2 (i == 2)
            // ---------- Свёртка от начала участка суммирования до конца участка суммирования.
            // --- на frame1
            for (int i = 0; i < add_min2; i++)
                signal[begSignal + i] = frame1.data[begFrame + i] + frame2.data[i];
            begFrame += add_min2;
            begSignal += add_min2;
            // --- на frame2 
            for (int i = 0; i < add_min2; i++)
                // Записываем данные frame2 от левого хвоста (со сдвигом на add_min2), свёрнутые с хвостом frame1
                frame2.data[add_min2 + i] = frame1.data[begFrame + i] + frame2.data[add_min2 + i];
            break;

        default:
            // Если кадр промежуточный (не первый, не второй и не последний)
            if (j != (frames_count - 1))
            {
                // Сменить полусвёрнутый frame2 на frame1
                Frame::switchData(frame1, frame2);
                // Скопировать данные в frame2.data и фильтровать
                frame2.setData(signal, frame_size * j, frame_size);
                fft_filtrate(frame2);

                /*
                Свёртка frame1 и frame2, и запись кадра frame1 в Signal,
                frame2 выйдет со свёрнутым левым хвостом, далее его логическое
                начало должно быть сдвинуто на add_min2.
                */

                // Переписываем frame1 (со сдвигом начала на add_min2) до начала свёртки
                begFrame = add_min2;
                for (int i = 0; i < frame_size - add_min2; i++)
                    signal[begSignal + i] = frame1.data[begFrame + i];
                begFrame += frame_size - add_min2;
                begSignal += frame_size - add_min2;
                // ---------- Свёртка от начала участка суммирования до конца участка суммирования.
                // --- на frame1
                for (int i = 0; i < add_min2; i++)
                    signal[begSignal + i] = frame1.data[begFrame + i] + frame2.data[i];
                begFrame += add_min2;
                begSignal += add_min2;
                // --- на frame2
                for (int i = 0; i < add_min2; i++)
                    frame2.data[add_min2 + i] = frame1.data[begFrame + i] + frame2.data[add_min2 + i];
            }
            // Если кадр последний
            else {
                // Сменить полусвёрнутый frame2 на frame1
                Frame::switchData(frame1, frame2);
                // Текущая позиция в исходном signal
                int beginOrigin = frame_size * j;
                // Сколько элементов осталось записать в Signal из frame2 в последнем кадре
                int res = length - beginOrigin;
                printf("res: %d\n", res);

                frame2.setData(signal, beginOrigin, res);
                fft_filtrate(frame2);

                /*
                Свёртка frame1 и frame2, и запись кадра frame1 в Signal,
                frame2 выйдет со свёрнутым левым хвостом, далее его логическое
                начало должно быть сдвинуто на add_min2.
                */

                // Переписываем frame1 (со сдвигом начала на add_min2) до начала свёртки
                begFrame = add_min2;
                for (int i = 0; i < frame_size - add_min2; i++)
                    signal[begSignal + i] = frame1.data[begFrame + i];
                begFrame += frame_size - add_min2;
                begSignal += frame_size - add_min2;

                // ---------- Свёртка от начала участка суммирования до конца участка суммирования.
                // --- на frame1
                int begFrame2 = 0;

                for (int i = 0; i < add_min2; i++)
                    signal[begSignal + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
                begFrame += add_min2;
                begFrame2 += add_min2;
                begSignal += add_min2;

                // --- на frame2

                // Если оставшееся количество элементов меньше длины хвоста
                if (res >= add_min2) {
                    // Свертка сигнала frame2 с хвостом frame1
                    for (int i = 0; i < add_min2; i++)
                        signal[begSignal + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
                    begFrame2 += add_min2;
                    begSignal += add_min2;

                    // Запись сигнала frame2 без хвоста
                    for (int i = 0; i < abs(res - add_min2); i++)
                        signal[begSignal + i] = frame2.data[begFrame2 + i];
                    begFrame2 += abs(res - add_min2);

                    if (right_tail != nullptr)
                    {
                        // Запись правого хвоста (add_min2) в right_tail
                        for (int i = 0; i < add_min2; i++)
                            right_tail[i] = frame2.data[begFrame2 + i];
                    }
                }
                else {
                    // Свертка сигнала frame2 с хвостом frame1
                    for (int i = 0; i < res; i++)
                        signal[begSignal + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
                    begFrame += res;
                    begFrame2 += res;

                    if (right_tail != nullptr)
                    {
                        // Запись правого хвоста (add_min2) в right_tail
                        for (int i = 0; i < add_min2; i++)
                            if (i < add_min2 - res)
                                right_tail[i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
                            else
                                right_tail[i] = frame2.data[begFrame2 + i];
                    }
                }
            }
            break;
        }
    }
}

int FLT_FilterFile::get_frame_size()
{
    return frame_size;
}