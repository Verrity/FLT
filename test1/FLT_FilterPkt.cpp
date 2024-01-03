#include "FLT_FilterPkt.h"
#include <iostream>

bool FLT_FilterPkt::startTransfer(int packet_length, int accurancy)
{
	this->accurancy = 1;
	if (packet_length < 1) {
		error_code = FILTER_ERROR_LENGTH;
		return false;
	}
	if (accurancy < 0) {
		error_code = FILTER_ERROR_ACCURANCY;
		return false;
	}

	this->N = N;
	this->accurancy = accurancy;
	this->packet_size = packet_length;
	add_min = N - 1;

	fft_size = 1;
	while (fft_size < (packet_size + add_min)) // fft_size < (длина сигнала + минимум добавочных элементов)
	{
		fft_size = fft_size << 1;
	}

	frame_size = packet_size;
	fft_size = fft_size << accurancy;

	if (isMinAllocated)
		free_min();
	init_min(N, fd, fft_size, packet_length);

	return true;
}

bool FLT_FilterPkt::filtratePkt(double* packet)
{
	packet_index++;
	// Если пакет первый
	switch (packet_index) {
	case 1:
		frame1.setData(packet, 0, packet_size);
		fft_filtrate(frame1);
		error_code = FILTER_FIRST_PKT;
		return false;
		break;
	default:
		frame2.setData(packet, 0, packet_size);
		fft_filtrate(frame2);

		for (int i = 0; i < frame_size - add_min2; i++) {
			packet[i] = frame1.data[add_min2 + i];
		}
		for (int i = 0; i < add_min2; i++) {
			packet[frame_size - add_min2 + i] = frame1.data[frame_size + i] + frame2.data[i];
		}
		return true;
		Frame::switchData(frame1, frame2);
		break;
	}
}

double* FLT_FilterPkt::getLatestPkt()
{
	double* lastPacket = new double[packet_size];

	for (int i = 0; i < frame_size - add_min2; i++) {
		lastPacket[i] = frame1.data[add_min2 + i];
	}
	for (int i = 0; i < add_min2; i++) {
		lastPacket[frame_size - add_min2 + i] = frame1.data[frame_size + i] + frame2.data[i];
	}
	return lastPacket;
}

bool FLT_FilterPkt::stopTransfer()
{
	accurancy = 1;
	min_fft = 1;
	packet_size = 0;
	frame_size = 0;
	packet_index = 0;
	add_min2 = 0;
	return true;
}

bool FLT_FilterPkt::startTransferBlock(int packet_size, int length_parameter, int accurancy)
{
	this->accurancy = 0;
	min_fft = 1;
	this->packet_size = 0;

	if (!check_accurancy(accurancy))
		return false;
	if (!check_length_parameter(length_parameter))
		return false;

	frame_size = N;
	fft_size = 1;
	while (fft_size < (frame_size + add_min)) // fft_size < (длина сигнала + минимум добавочных элементов)
	{
		fft_size = fft_size << 1;
	}

	fft_size = fft_size << length_parameter;
	frame_size = fft_size - N + 1;
	fft_size = fft_size << accurancy;

	add_min = N - 1;
	add_min2 = add_min / 2;

	if ((packet_size / 4) < (frame_size)) { // Если длина меньше, чем 4 кадра
		error_code = FILTER_ERROR_ACCURANCY;
		return false;
	}

	//---// Максимально близкое значение БПФ к числу N справа
	//---while (min_fft < N) {
	//---	min_fft = min_fft << 1;
	//---}

	//---this->accurancy = accurancy;
	//---// Определяем размер БПФ
	//---fft_size = min_fft << accurancy;
	//---add_min = N - 1;
	//---add_min2 = add_min / 2;

	//---if ((packet_size / 4) < (fft_size - add_min)){ // Если длина меньше, чем 4 кадра
	//---	error_code = FILTER_ERROR_ACCURANCY;
	//---	return false;
	//---}
	//---frame_size = fft_size - add_min;
	this->packet_size = packet_size;
	packet_index = 0;
	printf("\n\n======= Local Transfer parameters =======\n\n");
	printf("add_min = %d\n", add_min);
	printf("add_min2 = %d\n", add_min2);
	printf("fft_size = %d\n", fft_size);
	printf("frame_size = %d\n", frame_size);

	// Освобождается в этом Классе ---------------
	if (right_tail0 != nullptr)	delete[] right_tail0;
	if (left_tail != nullptr)	delete[] left_tail;
	if (right_tail != nullptr)	delete[] right_tail;
	if (ptrToAllocatedData1 != nullptr)	delete[] ptrToAllocatedData1;
	if (ptrToAllocatedData2 != nullptr)	delete[] ptrToAllocatedData2;

	right_tail0 = new double[add_min2];
	left_tail = new double[add_min2];
	right_tail = new double[add_min2];
	ptrToAllocatedData1 = new double[packet_size];
	ptrToAllocatedData2 = new double[packet_size];
	// ------------------------------------------

	// Освобождается в родительском классе ------
	if (isMinAllocated)
		free_min();
	init_min(N, fd, fft_size, frame_size);
	// ------------------------------------------
	return true;
}

bool FLT_FilterPkt::filtratePktBlock(double* const packet)
{
	packet_index++;

	if (packet_index == 1)
	{
		// Первый кадр не будет возвращен, поэтому его требуеся сохранить
		packet1 = ptrToAllocatedData1;
		packet2 = ptrToAllocatedData2;

		memcpy(packet1, packet, sizeof(double) * packet_size);
		_filtrateBlock(packet1, packet_size, nullptr, right_tail0);

		return false;
	}
	else {
		memcpy(packet2, packet, sizeof(double) * packet_size);
		_filtrateBlock(packet2, packet_size, left_tail, right_tail);

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

double* FLT_FilterPkt::getLatestPktBlock()
{
	if (packet_index == 0) {
		error_code = FILTER_ERROR_FUNCTION;
		return nullptr;
	}
	else {
		double* lastPacket = new double[packet_size];
		memcpy(lastPacket, packet1, sizeof(double) * packet_size);
		return lastPacket;
	}
}

bool FLT_FilterPkt::stopTransferBlock()
{
	if (packet_index == 0) {
		error_code = FILTER_ERROR_FUNCTION;
		return false;
	}

	accurancy = 1;
	min_fft = 1;
	packet_size = 0;
	frame_size = 0;
	packet_index = 0;
	add_min2 = 0;

	if (right_tail0 != nullptr)			delete[] right_tail0;			right_tail0 = nullptr;
	if (left_tail != nullptr)			delete[] left_tail;				left_tail = nullptr;
	if (right_tail != nullptr)			delete[] right_tail;			right_tail = nullptr;
	if (ptrToAllocatedData1 != nullptr)	delete[] ptrToAllocatedData1;	ptrToAllocatedData1 = nullptr;
	if (ptrToAllocatedData2 != nullptr)	delete[] ptrToAllocatedData2;	ptrToAllocatedData2 = nullptr;
	packet1 = nullptr;
	packet2 = nullptr;
	return true;
}

int FLT_FilterPkt::measureAttenuation(double*& pointer, int packet_size, int accurancy, double f_low, double step, double f_high, unsigned long& time_ns)
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
	fft_size << accurancy;
	frame_size = fft_size - N + 1;

	if ((packet_size / 4) < (fft_size - add_min)) { // Если длина меньше, чем 4 кадра
		error_code = FILTER_ERROR_ACCURANCY;
		return 0;
	}

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
		// БПФ
		if (f == AM_f_low)
			startTimer();
		fftw_execute_dft_r2c(forward_h, signal, signalFFT);
		// Умнодить БПФ сигнала и БПФ Имп. хар-ки. (фильтрация)
		for (int i = 0; i < fft_size / 2 + 1; i++) {
			double& a = signalFFT[i][0];
			double& b = signalFFT[i][1];
			double& c = h_fft[i][0];
			double& d = h_fft[i][1];

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

