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
    if (    // ���� ������� ��� ��� ��������� ��� �� ��������������, ����������
        (signal_size != length) ||
        (this->accurancy != accurancy)
        )
    {
        signal_size = length;
        frame_size = N * 45;
        fft_size = 1;
        while (fft_size < (frame_size + add_min)) // fft_size < (����� ������� + ������� ���������� ���������)
        {
            fft_size = fft_size << 1;
        }
        fft_size << accurancy - 1;
        frame_size = fft_size - N + 1;

        // ���� � ����� ������� ������������ 3 ��� ������ ������
        // ������������ ������� ����������, ����� �������
        if (length <= N * 45 * 3) {
            fft_size = 1;
            while (fft_size < (length + add_min)) // fft_size < (����� ������� + ������� ���������� ���������)
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

    // ���� � ����� ������� ������������ 3 ��� ������ ������
    // ������������ ������� ����������, ����� �������
    // 3 - ������ ��� frame_size ������ 3 ����, ����� ������� ���������� ��������
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
    if (    // ���� ������� ��� ��� ��������� ��� �� ��������������, ����������
        (signal_size != length) ||
        (this->accurancy != accurancy)
        )
    {
        signal_size = length;
        this->accurancy = accurancy;
        frame_size = N * 45;
        fft_size = 1;
        while (fft_size < (frame_size + add_min)) // fft_size < (����� ������� + ������� ���������� ���������)
        {
            fft_size = fft_size << 1;
        }
        fft_size << accurancy - 1;
        frame_size = fft_size - N + 1;

        left_tail = new double[add_min2];
        right_tail = new double[add_min2];
    }

    // ���� � ����� ������� ������������ 3 ��� ������ ������
    // ������������ ������� ����������, ����� �������
    // 3 - ������ ��� frame_size ������ 3 ����, ����� ������� ���������� ��������
    if (length <= frame_size * 3) {
        fft_size = 1;
        while (fft_size < (length + add_min)) // fft_size < (����� ������� + ������� ���������� ���������)
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
    int begFrame = 0;	// ������� ������� � frame1
    int begSignal = 0;	// ������� ������� � ����������� Signal

    for (int j = 0; j < frames_count; j++) {
        // ------------------ ������ � ���������� ������
        switch (j)
        {
        case 0:	// ���� ���� ������
            // ����������� ������ � frame1.data � �����������
            frame1.setData(signal, 0, frame_size);
            fft_filtrate(frame1);
            /*
                ������ frame1 � frame2, � ������ ����� frame1 ��� ������ ������
                � Signal, frame2 ������ �� �������� ����� �������,
                ����� ��� ���������� ������ ������ ���� �������� �� add_min2.
            */
            if (left_tail != nullptr)
            {
                // ���������� ����� ����� � left_tail
                for (int i = 0; i < add_min2; i++)
                    left_tail[i] = frame1.data[i];
            }
            begFrame += add_min2;

            // ������������ �� ����� ������ ������ (��������) �� ������� ������������
            for (int i = 0; i < frame_size - add_min2; i++)
                signal[i] = frame1.data[begFrame + i];
            begFrame += frame_size - add_min2;
            begSignal += frame_size - add_min2;
            break;

        case 1:	// ���� ���� ������, �� frame2 ��������� ��������� � ��������� ������ � frame1 (i == 1) � frame2 (i == 2)
            // ����������� ������ � frame2.data � �����������
            frame2.setData(signal, frame_size, frame_size);
            fft_filtrate(frame2);
            // �������� frame2, ����� ��������� ������ frame1 (i == 1) � frame2 (i == 2)
            // ---------- ������ �� ������ ������� ������������ �� ����� ������� ������������.
            // --- �� frame1
            for (int i = 0; i < add_min2; i++)
                signal[begSignal + i] = frame1.data[begFrame + i] + frame2.data[i];
            begFrame += add_min2;
            begSignal += add_min2;
            // --- �� frame2 
            for (int i = 0; i < add_min2; i++)
                // ���������� ������ frame2 �� ������ ������ (�� ������� �� add_min2), �������� � ������� frame1
                frame2.data[add_min2 + i] = frame1.data[begFrame + i] + frame2.data[add_min2 + i];
            break;

        default:
            // ���� ���� ������������� (�� ������, �� ������ � �� ���������)
            if (j != (frames_count - 1))
            {
                // ������� ������������ frame2 �� frame1
                Frame::switchData(frame1, frame2);
                // ����������� ������ � frame2.data � �����������
                frame2.setData(signal, frame_size * j, frame_size);
                fft_filtrate(frame2);

                /*
                ������ frame1 � frame2, � ������ ����� frame1 � Signal,
                frame2 ������ �� �������� ����� �������, ����� ��� ����������
                ������ ������ ���� �������� �� add_min2.
                */

                // ������������ frame1 (�� ������� ������ �� add_min2) �� ������ ������
                begFrame = add_min2;
                for (int i = 0; i < frame_size - add_min2; i++)
                    signal[begSignal + i] = frame1.data[begFrame + i];
                begFrame += frame_size - add_min2;
                begSignal += frame_size - add_min2;
                // ---------- ������ �� ������ ������� ������������ �� ����� ������� ������������.
                // --- �� frame1
                for (int i = 0; i < add_min2; i++)
                    signal[begSignal + i] = frame1.data[begFrame + i] + frame2.data[i];
                begFrame += add_min2;
                begSignal += add_min2;
                // --- �� frame2
                for (int i = 0; i < add_min2; i++)
                    frame2.data[add_min2 + i] = frame1.data[begFrame + i] + frame2.data[add_min2 + i];
            }
            // ���� ���� ���������
            else {
                // ������� ������������ frame2 �� frame1
                Frame::switchData(frame1, frame2);
                // ������� ������� � �������� signal
                int beginOrigin = frame_size * j;
                // ������� ��������� �������� �������� � Signal �� frame2 � ��������� �����
                int res = length - beginOrigin;
                printf("res: %d\n", res);

                frame2.setData(signal, beginOrigin, res);
                fft_filtrate(frame2);

                /*
                ������ frame1 � frame2, � ������ ����� frame1 � Signal,
                frame2 ������ �� �������� ����� �������, ����� ��� ����������
                ������ ������ ���� �������� �� add_min2.
                */

                // ������������ frame1 (�� ������� ������ �� add_min2) �� ������ ������
                begFrame = add_min2;
                for (int i = 0; i < frame_size - add_min2; i++)
                    signal[begSignal + i] = frame1.data[begFrame + i];
                begFrame += frame_size - add_min2;
                begSignal += frame_size - add_min2;

                // ---------- ������ �� ������ ������� ������������ �� ����� ������� ������������.
                // --- �� frame1
                int begFrame2 = 0;

                for (int i = 0; i < add_min2; i++)
                    signal[begSignal + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
                begFrame += add_min2;
                begFrame2 += add_min2;
                begSignal += add_min2;

                // --- �� frame2

                // ���� ���������� ���������� ��������� ������ ����� ������
                if (res >= add_min2) {
                    // ������� ������� frame2 � ������� frame1
                    for (int i = 0; i < add_min2; i++)
                        signal[begSignal + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
                    begFrame2 += add_min2;
                    begSignal += add_min2;

                    // ������ ������� frame2 ��� ������
                    for (int i = 0; i < abs(res - add_min2); i++)
                        signal[begSignal + i] = frame2.data[begFrame2 + i];
                    begFrame2 += abs(res - add_min2);

                    if (right_tail != nullptr)
                    {
                        // ������ ������� ������ (add_min2) � right_tail
                        for (int i = 0; i < add_min2; i++)
                            right_tail[i] = frame2.data[begFrame2 + i];
                    }
                }
                else {
                    // ������� ������� frame2 � ������� frame1
                    for (int i = 0; i < res; i++)
                        signal[begSignal + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
                    begFrame += res;
                    begFrame2 += res;

                    if (right_tail != nullptr)
                    {
                        // ������ ������� ������ (add_min2) � right_tail
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