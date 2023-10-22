#include "FLT_FilterFile.h"
#include <iostream>

FLT_FilterFile::FLT_FilterFile()
{
}

FLT_FilterFile::~FLT_FilterFile()
{
    if (h != nullptr)               delete[] h;
    if (freq_match != nullptr)      delete[] freq_match;
    if (h_fft != nullptr)           fftw_free(h_fft);
    if (mul_frames_fft != nullptr)  fftw_free(mul_frames_fft);
    if (conv_frames != nullptr)     delete[] conv_frames;
    fftw_destroy_plan(forward_signal_fft);
    fftw_destroy_plan(backward_signalF);
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

int FLT_FilterFile::filtrateBlock(double* in, int length, double* &out, int accurancy, bool tails) {
    // Check Parameters ---------------------------------------------
    if (length <= 0) {
        error_code = FILTER_ERROR_LENGTH;
        return 0;
    }
    if (!check_accurancy(accurancy))
        return 0;
    // --------------------------------------------------------------
    sample_size = N * 45;
    fft_size = 1;
    while (fft_size < (sample_size + add_min)) // fft_size < (длина сигнала + минимум добавочных элементов)
    {
        fft_size = fft_size << 1;
    }
    fft_size << accurancy - 1;
    int len = 0;
    sample_size = fft_size - N + 1;

    // Если сигнал меньше стандартной выборки
    if ((length + add_min) < fft_size) {
        return filtrate(in, length, out, accurancy, tails);
    }
    else {
        // Calc Parameters ----------------------------------------------
        if ( // Если массивы под эти параметры не были рассчитаны рассчитать новые
            (filtrateBlock_length != length) ||
            (filtrateBlock_accurancy != accurancy)
            )
        {
            filtrateBlock_length = length;
            filtrateBlock_accurancy = accurancy;

            this->accurancy = accurancy;
            add_min = N - 1;

            conv_size = sample_size * 2 + add_min;

            if (h != nullptr)               delete[] h;
            if (freq_match != nullptr)      delete[] freq_match;
            if (h_fft != nullptr)           fftw_free(h_fft);
            if (mul_frames_fft != nullptr)  fftw_free(mul_frames_fft);
            if (conv_frames != nullptr)     delete[] conv_frames;
            fftw_destroy_plan(forward_signal_fft);
            fftw_destroy_plan(backward_signalF);

            h = new double[fft_size];
            conv_frames = new double[conv_size];
            h_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (fft_size / 2 + 1));
            mul_frames_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (fft_size / 2 + 1) * 2);
            freq_match = new double[fft_size];

            frame1.init(N, sample_size, fft_size);
            frame2.init(N, sample_size, fft_size);

            forward_signal_fft = fftw_plan_dft_r2c_1d(fft_size, pFrameData, pFrameDataFFT, FFTW_ESTIMATE);
            backward_signalF = fftw_plan_dft_c2r_1d(fft_size, mul_frames_fft, pFrameData, FFTW_ESTIMATE);

            if (window) {
                if (w != nullptr)   delete[] w;
                w = new double[fft_size];
                calc_window();
            }
            calc_h();
            calc_h_fft();
        }
        // --------------------------------------------------------------

        int frames_count = length % sample_size != 0 ? floor(double(length) / sample_size) + 1 : floor(double(length) / sample_size);

        int prev_begin = -sample_size;
        int curr_begin = 0;
        
        if (tails) {
            len = length + add_min;
        }
        else {
            len = length;
        }
        out = new double[len];

        for (int i = 1; i < frames_count; i++) {
            prev_begin = prev_begin + sample_size;
            curr_begin = curr_begin + sample_size;

            // Загрузить, обработать первый кадр
            if (i == 1) {
                frame1.setData(in, prev_begin, sample_size);
                fft_filtrate(frame1);
            }
            else {   // Если кадр не первый
                // M_real+L+1 = length(conv_frames) = fft_size
                for (int i = 0; i < fft_size; i++)
                    frame1.data[i] = conv_frames[sample_size + i];
            }

            // Обработать текущий кадр (CURRENT_FRAME)
            if (i == frames_count - 1) {// Если он последний
                int res = length - curr_begin;
                frame2.setData(in, curr_begin, res);
            }
            else {                      // Если он не последний
                frame2.setData(in, curr_begin, sample_size);
            }

            fft_filtrate(frame2);
            // Соеденить два кадра по алгоритму в массив conv_frames
            //printf_s("[%d] 1: %d\t2: %d\n", i, frame1.data_size, frame2.data_size);
            convolFull(frame1, frame2);

            // Записать в выходной массив
            if (i == 1) {
                if (tails) {
                    for (int i = 0; i < add_min + frame1.data_size; i++) {
                        out[prev_begin + i] = conv_frames[i];
                    }
                }
                else {
                    for (int i = 0; i < frame1.data_size; i++) {
                        out[prev_begin + i] = conv_frames[add_min + i];
                    }
                }
            }
            else {
                for (int i = 0; i < sample_size; i++) {
                    out[prev_begin + i] = conv_frames[i];
                }
            }

            // Записать последний кадр в выходной массив
            if (i == frames_count - 1) {
                int res = 0;
                if (tails) {
                    // сколько элементов требуется записать
                    res = length - curr_begin;
                    for (int i = 0; i < res + add_min; i++) {
                        out[curr_begin + i] = conv_frames[sample_size + i];
                    }
                }
                else {
                    res = length - curr_begin;
                    for (int i = 0; i < res; i++) {
                        out[curr_begin + i] = conv_frames[sample_size + i];
                    }
                }
            }
        }
    }

    return len;
}

int FLT_FilterFile::get_sample_size()
{
    return sample_size;
}