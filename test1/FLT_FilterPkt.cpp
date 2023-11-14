#include "FLT_FilterPkt.h"
#include <iostream>

bool FLT_FilterPkt::startTransfer(int packet_length, int accurancy)
{
	if (packet_length < 1) {
		error_code = FILTER_ERROR_LENGTH;
		return false;
	}
	if (!check_accurancy(accurancy))
		return false;

	this->N = N;
	this->accurancy = accurancy;
	this->packet_size = packet_length;
	add_min = N - 1;

	fft_size = 1;
	while (fft_size < (packet_size + add_min)) // fft_size < (длина сигнала + минимум добавочных элементов)
	{
		fft_size = fft_size << 1;
	}
	fft_size << accurancy - 1;

	conv_size = packet_size * 2 + add_min;

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

	frame1.init(N, packet_size, fft_size);
	frame2.init(N, packet_size, fft_size);

	forward_signal_fft = fftw_plan_dft_r2c_1d(fft_size, pFrameData, pFrameDataFFT, FFTW_ESTIMATE);
	backward_signalF = fftw_plan_dft_c2r_1d(fft_size, mul_frames_fft, pFrameData, FFTW_ESTIMATE);

	if (window) {
		if (w != nullptr)   delete[] w;
		w = new double[fft_size];
		calc_window();
	}
	calc_h();
	calc_h_fft();

	return 1;
}

bool FLT_FilterPkt::filtratePkt1(double* packet)
{
	packet_index++;
	// Если пакет первый
	if (packet_index == 1) {
		frame1.setData(packet, 0, packet_size);
		fft_filtrate(frame1);
		error_code = FILTER_FIRST_PKT;
		return false;
	}
	else {
		frame2.setData(packet, 0, packet_size);
		fft_filtrate(frame2);
		convolFull(frame1, frame2);

		// Записать в выходной массив frame1.data
		for (int i = 0; i < packet_size; i++)
			packet[i] = conv_frames[i];

		// Переместить соединенный текущий массив в предыдущий массив
		for (int i = 0; i < fft_size; i++) // add_min + packet_size + 1 = packet_length(conv_frames) = fft_size
			frame1.data[i] = conv_frames[packet_size + i];

		return true;
	}
}

int FLT_FilterPkt::stopTransfer(double* &lastPacket)
{
	lastPacket = new double[packet_size + N - 1];

	for (int i = 0; i < packet_size + N - 1; i++)
		lastPacket[i] = conv_frames[packet_size + i];
	return packet_size + N - 1;

}

void FLT_FilterPkt::filtratePacketBlock(double* packet, double* left_tail, double* right_tail)
{
	int frames_count = packet_size % frame_size != 0 ? floor(double(packet_size) / frame_size) + 1 : floor(double(packet_size) / frame_size);
	int begFrame = 0;	// Текущая позиция в frame1
	int begPacket = 0;	// Текущая позиция в записывании Packet1

	for (int j = 0; j < frames_count; j++) {
		// ------------------ Свёртка и фильтрация кадров
		switch (j)
		{
		case 0:	// Если кадр первый
			// Скопировать данные в frame1.data и фильтровать
			frame1.setData(packet, 0, frame_size);
			fft_filtrate(frame1);
			/*
				Свёртка frame1 и frame2, и запись кадра frame1 без левого хвоста
				в Packet1, frame2 выйдет со свёрнутым левым хвостом,
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
				packet[i] = frame1.data[begFrame + i];
			begFrame += frame_size - add_min2;
			begPacket = frame_size - add_min2;
			break;

		case 1:	// Если кадр второй, то frame2 требуется загрузить и закончить свёртку с frame1 (i == 1) и frame2 (i == 2)
			// Скопировать данные в frame2.data и фильтровать
			frame2.setData(packet, frame_size, frame_size);
			fft_filtrate(frame2);
			// Получили frame2, можно закончить свёртку frame1 (i == 1) и frame2 (i == 2)
			// ---------- Свёртка от начала участка суммирования до конца участка суммирования.
			// --- на frame1
			for (int i = 0; i < add_min2; i++)
				packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[i];
			begFrame += add_min2;
			begPacket += add_min2;
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
				frame2.setData(packet, frame_size * j, frame_size);
				fft_filtrate(frame2);

				/*
				Свёртка frame1 и frame2, и запись кадра frame1 в Packet1,
				frame2 выйдет со свёрнутым левым хвостом, далее его логическое
				начало должно быть сдвинуто на add_min2.
				*/

				// Переписываем frame1 (со сдвигом начала на add_min2) до начала свёртки
				begFrame = add_min2;
				for (int i = 0; i < frame_size - add_min2; i++)
					packet[begPacket + i] = frame1.data[begFrame + i];
				begFrame += frame_size - add_min2;
				begPacket += frame_size - add_min2;
				// ---------- Свёртка от начала участка суммирования до конца участка суммирования.
				// --- на frame1
				for (int i = 0; i < add_min2; i++)
					packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[i];
				begFrame += add_min2;
				begPacket += add_min2;
				// --- на frame2
				for (int i = 0; i < add_min2; i++)
					frame2.data[add_min2 + i] = frame1.data[begFrame + i] + frame2.data[add_min2 + i];
			}
			// Если кадр последний
			else {
				// Сменить полусвёрнутый frame2 на frame1
				Frame::switchData(frame1, frame2);
				// Текущая позиция в исходном Packet
				int beginOrigin = frame_size * j;
				// Сколько элементов осталось записать в Packet1 из frame2 в последнем кадре
				int res = packet_size - beginOrigin;
				printf("res: %d\n", res);

				frame2.setData(packet, beginOrigin, res);
				fft_filtrate(frame2);

				/*
				Свёртка frame1 и frame2, и запись кадра frame1 в Packet1,
				frame2 выйдет со свёрнутым левым хвостом, далее его логическое
				начало должно быть сдвинуто на add_min2.
				*/

				// Переписываем frame1 (со сдвигом начала на add_min2) до начала свёртки
				begFrame = add_min2;
				for (int i = 0; i < frame_size - add_min2; i++)
					packet[begPacket + i] = frame1.data[begFrame + i];
				begFrame += frame_size - add_min2;
				begPacket += frame_size - add_min2;

				// ---------- Свёртка от начала участка суммирования до конца участка суммирования.
				// --- на frame1
				int begFrame2 = 0;

				for (int i = 0; i < add_min2; i++)
					packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
				begFrame += add_min2;
				begFrame2 += add_min2;
				begPacket += add_min2;

				// --- на frame2
				
				// Если оставшееся количество элементов меньше длины хвоста
				if (res >= add_min2) {
					// Свертка сигнала frame2 с хвостом frame1
					for (int i = 0; i < add_min2; i++)
						packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
					begFrame2 += add_min2;
					begPacket += add_min2;

					// Запись сигнала frame2 без хвоста
					for (int i = 0; i < abs(res - add_min2); i++)
						packet[begPacket + i] = frame2.data[begFrame2 + i];
					begFrame2 += abs(res - add_min2);

					// Запись правого хвоста (add_min2) в right_tail
					for (int i = 0; i < add_min2; i++)
						right_tail[i] = frame2.data[begFrame2 + i];
				}
				else {
					// Свертка сигнала frame2 с хвостом frame1
					for (int i = 0; i < res; i++)
						packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
					begFrame += res;
					begFrame2 += res;

					// Запись правого хвоста (add_min2) в right_tail
					for (int i = 0; i < add_min2; i++)
						if (i < add_min2 - res)
							right_tail[i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
						else
							right_tail[i] = frame2.data[begFrame2 + i];
				}
			}
			break;
		}
	}
}

bool FLT_FilterPkt::startTransferBlock(int packet_size, unsigned int value)
{
	this->value = 1;
	min_fft = 1;
	this->packet_size = 0;

	if (value < 0)
		return false;
	// Максимально близкое значение БПФ к числу N справа
	while (min_fft < N) {
		min_fft = min_fft << 1;
	}
	if (value > 0) 
		this->value = value;
	// Определяем размер БПФ
	fft_size = min_fft << this->value;
	add_min = N - 1;
	add_min2 = add_min / 2;

	if ((packet_size / 4) < (fft_size - add_min)){ // Если длина меньше, чем 4 кадра
		error_code = FILTER_ERROR_VALUE;
		return false;
	}

	frame_size = fft_size - add_min;
	this->packet_size = packet_size;
	packet_index = 0;
	printf("\n\n======= Local Transfer parameters =======\n\n");
	printf("add_min = %d\n", add_min);
	printf("add_min2 = %d\n", add_min2);
	printf("fft_size = %d\n", fft_size);
	printf("frame_size = %d\n", frame_size);

	// Освобождается в этом Классе ---------------
	if (right_tail0	!= nullptr)	delete[] right_tail0;
	if (left_tail	!= nullptr)	delete[] left_tail;
	if (right_tail	!= nullptr)	delete[] right_tail;
	if (ptrToAllocatedData1	!= nullptr)	delete[] ptrToAllocatedData1;
	if (ptrToAllocatedData2	!= nullptr)	delete[] ptrToAllocatedData2;

	right_tail0	= new double[add_min2];
	left_tail	= new double[add_min2];
	right_tail	= new double[add_min2];
	ptrToAllocatedData1	= new double[packet_size];
	ptrToAllocatedData2	= new double[packet_size];
	// ------------------------------------------

	// Освобождается в родительском классе ------
	if (h != nullptr)               delete[] h;
	if (h_fft != nullptr)           fftw_free(h_fft);
	if (mul_frames_fft != nullptr)  fftw_free(mul_frames_fft);
	if (conv_frames != nullptr)     delete[] conv_frames;
	fftw_destroy_plan(forward_signal_fft);
	fftw_destroy_plan(backward_signalF);

	h = new double[fft_size];
	h_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (fft_size / 2 + 1));
	mul_frames_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (fft_size / 2 + 1) * 2);
	
	frame1.init(N, frame_size, fft_size);
	frame2.init(N, frame_size, fft_size);

	forward_signal_fft = fftw_plan_dft_r2c_1d(fft_size, pFrameData, pFrameDataFFT, FFTW_ESTIMATE);
	backward_signalF = fftw_plan_dft_c2r_1d(fft_size, mul_frames_fft, pFrameData, FFTW_ESTIMATE);

	if (window) {
		if (w != nullptr)   delete[] w;
		w = new double[fft_size];
		calc_window();
	}
	calc_h();
	calc_h_fft();
	// ------------------------------------------
	return true;
}

bool FLT_FilterPkt::filtratePktBlock1(double* const packet)
{
	packet_index++;

	if (packet_index == 1)
	{
		// Первый кадр не будет возвращен, поэтому его требуеся сохранить
		packet1 = ptrToAllocatedData1;
		packet2 = ptrToAllocatedData2;

		memcpy(packet1, packet, sizeof(double) * packet_size);
		filtratePacketBlock(packet1, nullptr, right_tail0);

		return false;
	}
	else {
		memcpy(packet2, packet, sizeof(double) * packet_size);
		filtratePacketBlock(packet2, left_tail, right_tail);

		// Оба пакета фильтрованы, требуется их свернуть
		int begPacket1 = packet_size - add_min2;
		for (int i = 0; i < add_min2; i++)
		{
			// Сворачиваем правую часть packet1 с левым хвостом packet2
			packet1[begPacket1 + i] = packet1[begPacket1 + i] + left_tail[i];
			// Сворачиваем левую часть packet2 с правым хвостом packet1
			packet2[i] = packet2[i] + right_tail0[i];
		}

		memcpy(packet, packet1, sizeof(double) * packet_size);
		
		double* temp0 = packet1;
		packet1 = packet2;
		packet2 = temp0;
		// Правый хвост Пакета 2 нужно присвоить правому хвосту пакета 1, т.к. Пакет 2 --> Пакет 1
		double* temp = right_tail0;
		right_tail0 = right_tail;
		right_tail = temp;

		return true;
	}
}

double* const FLT_FilterPkt::filtratePktBlock2(double* const packet)
{
	//packet_index++;
	//int begPacket1 = packet_size - add_min;

	//switch (packet_index)
	//{
	//case 1:
	//	// Первый кадр не будет возвращен, поэтому его требуеся сохранить
	//	packet1 = ptrToAllocatedData1;
	//	packet2 = ptrToAllocatedData2;
	//	memcpy(packet1, packet, sizeof(double) * packet_size);

	//	filtrateFirstPacket(packet1);
	//	return nullptr;
	//	break;

	//case 2:
	//	memcpy(packet2, packet, sizeof(double) * packet_size);
	//	filtratePacketBlock(packet2);
	//	// Оба пакета фильтрованы, требуется их свернуть
	//	for (int i = 0; i < add_min2; i++)
	//	{
	//		// Сворачиваем правую часть packet1 с левым хвостом packet2
	//		packet1[begPacket1 + i] = packet1[begPacket1 + i] + left_tail[i];
	//		// Сворачиваем левую часть packet2 с правым хвостом packet1
	//		packet2[i] = packet2[i] + right_tail0[i];
	//	}
	//	return packet1;
	//	break;

	//default:
	//	packet1 = packet2;
	//	packet2 = ptrToAllocatedData1;

	//	memcpy(packet2, packet, sizeof(double) * packet_size);
	//	filtratePacketBlock(packet2);
	//	// Оба пакета фильтрованы, требуется их свернуть
	//	for (int i = 0; i < add_min2; i++)
	//	{
	//		// Сворачиваем правую часть packet1 с левым хвостом packet2
	//		packet1[begPacket1 + i] = packet1[begPacket1 + i] + left_tail[i];
	//		// Сворачиваем левую часть packet2 с правым хвостом packet1
	//		packet2[i] = packet2[i] + right_tail0[i];
	//	}
	//	return packet1;

	//	break;
	//}
	return 0;
}

double* FLT_FilterPkt::getLatestPktBlock1()
{
	double* lastPacket = new double[packet_size];
	memcpy(lastPacket, packet1, sizeof(double) * packet_size);
	return lastPacket;
}

void FLT_FilterPkt::stopTransferBlock()
{
	value = 1;
	min_fft = 1;
	packet_size	= 0;
	frame_size = 0;
	packet_index = 0;
	add_min2 = 0;

	if (right_tail0	!= nullptr)			delete[] right_tail0;			right_tail0	= nullptr;
	if (left_tail	!= nullptr)			delete[] left_tail;				left_tail	= nullptr;
	if (right_tail	!= nullptr)			delete[] right_tail;			right_tail	= nullptr;
	if (ptrToAllocatedData1 != nullptr)	delete[] ptrToAllocatedData1;	ptrToAllocatedData1 = nullptr;
	if (ptrToAllocatedData2 != nullptr)	delete[] ptrToAllocatedData2;	ptrToAllocatedData2 = nullptr;
	packet1 = nullptr;
	packet2 = nullptr;
}

