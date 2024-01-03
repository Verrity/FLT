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

bool FLT_FilterFile::filtrateBlock(double* const signal, int length, int length_parameter, int accurancy)
{
    if (    // Если массивы под эти параметры еще не рассчитывались, рассчитать
        (signal_size != length) ||
        (this->accurancy != accurancy)
        )
    {
        if (!check_accurancy(accurancy))
            return false;
        if (!check_length_parameter(length_parameter))
            return false;
        signal_size = length;
        frame_size = N;
        fft_size = 1;
        while (fft_size < (frame_size + add_min)) // fft_size < (длина сигнала + минимум добавочных элементов)
        {
            fft_size = fft_size << 1;
        }
        fft_size = fft_size << length_parameter;
        frame_size = fft_size - N + 1;
        fft_size = fft_size << accurancy;

        // Если в длину сигнала укладывается 3 или меньше кадров
        // Использовать обычную фильтрацию, иначе блочную
        if (length <= frame_size * 3) {
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

    if (signal != nullptr) {
        // Если в длину сигнала укладывается 3 или меньше кадров
        // Использовать обычную фильтрацию, иначе блочную
        // 3 - потому что frame_size берётся 3 раза, потом берутся оставшиеся элементы
        if (length <= frame_size * 3) {
            filtrate(signal, length, accurancy);
        }
        else {
            _filtrateBlock(signal, length, nullptr, nullptr);
        }
    }
    return true;
}

int FLT_FilterFile::filtrateBlockT(double* const signal, double*& output, int length, int length_parameter, int accurancy)
{
    if (    // Если массивы под эти параметры еще не рассчитывались, рассчитать
        (signal_size != length) ||
        (this->accurancy != accurancy)
        )
    {
        if (!check_accurancy(accurancy))
            return false;
        if (!check_length_parameter(length_parameter))
            return false;

        signal_size = length;
        this->accurancy = accurancy;
        //frame_size = N * 45;
        //fft_size = 1;
        //while (fft_size < (frame_size + add_min)) // fft_size < (длина сигнала + минимум добавочных элементов)
        //{
        //    fft_size = fft_size << 1;
        //}
        //fft_size << accurancy - 1;
        //frame_size = fft_size - N + 1;

        frame_size = N;
        fft_size = 1;
        while (fft_size < (frame_size + add_min)) // fft_size < (длина сигнала + минимум добавочных элементов)
        {
            fft_size = fft_size << 1;
        }
        fft_size = fft_size << length_parameter;
        frame_size = fft_size - N + 1;
        fft_size = fft_size << accurancy;

        left_tail = new double[add_min2];
        right_tail = new double[add_min2];
    }
    int size = 0;
    if (signal != nullptr) {
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

            size = filtrateT(signal, length, output, accurancy);

            return size;
        }
        else {
            if (isMinAllocated)
                free_min();
            init_min(N, fd, fft_size, frame_size);

            output = new double[length + add_min];
            double* ptrToOutBeginWriteSignal = output + add_min2;
            memcpy(ptrToOutBeginWriteSignal, signal, sizeof(double) * length);

            _filtrateBlock(ptrToOutBeginWriteSignal, length, left_tail, right_tail);

            for (int i = 0; i < add_min2; i++)
                output[i] = left_tail[i];
            for (int i = 0; i < add_min2; i++)
                output[add_min2 + length + i] = right_tail[i];
            return size;
        }
    }
}



bool FLT_FilterFile::check_length_parameter(int length_parameter)
{
    if (length_parameter < 0) {
        error_code = FILTER_ERROR_LENGTH_PARAMETER;
        return false;
    }
    return true;
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

int FLT_FilterFile::measureAttenuation(double*& pointer, int length, int accurancy, double f_low, double step, double f_high, unsigned long& time_ns)
{
    // Проверки границ частоты
    if ((f_low <= 0) || (step < f_low) || (f_high <= step) || (f_high > fd / 2)) {
        error_code = FILTER_ERROR_BAND;
        return 0;
    }
    if (!check_accurancy(accurancy))
        return 0;

    double AM_f_low = f_low;
    double AM_f_high = f_high;
    double AM_step = step;
    int AM_len = (AM_f_high - AM_f_low) / AM_step + 1;

    frame_size = N;
    fft_size = 1;
    while (fft_size < (frame_size + add_min)) // fft_size < (длина сигнала + минимум добавочных элементов)
    {
        fft_size = fft_size << 1;
    }
    fft_size = fft_size << accurancy;
    frame_size = fft_size - N + 1;

    if (isMinAllocated)
        free_min();
    init_min(this->N, this->fd, this->fft_size, 0);

    double* signal = new double[fft_size];
    pointer = new double[AM_len];
    fftw_complex* signalFFT = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (fft_size / 2 + 1));
    fftw_complex* signalFiltratedFFT = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (fft_size / 2 + 1) * 2);
    // Расчет планов БПФ
    fftw_plan forward_h = fftw_plan_dft_r2c_1d(fft_size, signal, signalFFT, FFTW_ESTIMATE);
    fftw_plan backward_signalF = fftw_plan_dft_c2r_1d(fft_size, signalFiltratedFFT, signal, FFTW_ESTIMATE);

    double averageOUT = 0;
    double averageIN = 0.774596669241; // Опорное значение 0 db [0.775 Вольт]
    double sum = 0;
    double H = 0;
    int iterator = 0;

    for (double f = AM_f_low; f <= AM_f_high; f += AM_step) {
        averageOUT = 0;
        sum = 0;
        // Создать синус для частоты f
        for (int j = 0; j < fft_size; j++)
            signal[j] = averageIN * sqrt(2) * sin(2 * FILTER_PI * j * f / fd);
        if (f == AM_f_low)
            startTimer();
        // БПФ
        fftw_execute_dft_r2c(forward_h, signal, signalFFT);
        // Умнодить БПФ сигнала и БПФ Имп. хар-ки. (фильтрация)
        for (int i = 0; i < fft_size / 2 + 1; i++) {
            double a = signalFFT[i][0];
            double b = signalFFT[i][1];
            double c = h_fft[i][0];
            double d = h_fft[i][1];

            signalFiltratedFFT[i][0] = ((a * c) - (b * d)) / fft_size;
            signalFiltratedFFT[i][1] = ((a * d) + (b * c)) / fft_size;
        }
        // ОБПФ
        fftw_execute_dft_c2r(backward_signalF, signalFiltratedFFT, signal);
        if (f == AM_f_low) {
            stopTimer();
            time_ns = m_int_durationTime_ns / fft_size;
        }
        // Сумма модулей сигналов
        for (int j = 0; j < fft_size; j++)
            sum = sum + std::abs(signal[j]);

        averageOUT = sum / fft_size;            // Среднее значение
        H = averageOUT / averageIN;             // Коэффициент передачи фильтра
        pointer[iterator] = 20 * log10(H);      // Затухание в dBu
        iterator++;
    }
    delete[] signal;
    fftw_free(signalFFT);
    fftw_free(signalFiltratedFFT);
    fftw_destroy_plan(forward_h);
    fftw_destroy_plan(backward_signalF);
    return AM_len;
}

int FLT_FilterFile::get_frame_size()
{
    return frame_size;
}