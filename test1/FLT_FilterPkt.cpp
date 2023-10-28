#include "FLT_FilterPkt.h"

//bool FLT_FilterPkt::setParameters(int N, int accurancy, int length)
//{
//	if (!check_N(N) || !check_accurancy(accurancy)) {
//		return 0;
//	}
//	if (length < 1) {
//		error_code = FILTER_ERROR_LENGTH;
//		return 0;
//	}
//
//	this->N = N;
//	this->accurancy = accurancy;
//	this->packet_size = length;
//	add_min = N - 1;
//
//	fft_size = 1;
//	while (fft_size < (packet_size + add_min)) // fft_size < (длина сигнала + минимум добавочных элементов)
//	{
//		fft_size = fft_size << 1;
//	}
//	fft_size << accurancy - 1;
//
//	conv_size = packet_size * 2 + add_min;
//}



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
		{
			packet[i] = conv_frames[i];
			//printf("%f\n", packet[i]);
		}

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

// ================================================================

int FLT_FilterPkt::filtrateBlock(double* packet, bool lastPacket)
{
	int frames_count = packet_size % frame_size != 0 ? floor(double(packet_size) / frame_size) + 1 : floor(double(packet_size) / frame_size);

	int prev_begin = -frame_size;
	int curr_begin = 0;

	for (int i = 1; i < frames_count; i++) {
		prev_begin = prev_begin + packet_size;
		curr_begin = curr_begin + packet_size;

		// Загрузить, обработать первый кадр
		if (i == 1) {
			frame1.setData(packet, prev_begin, packet_size);
			fft_filtrate(frame1);
		}
		else {   // Если кадр не первый
			// M_real+L+1 = length(conv_frames) = fft_size
			for (int i = 0; i < fft_size; i++)
				frame1.data[i] = conv_frames[packet_size + i];
		}

		// Обработать текущий кадр (CURRENT_FRAME)
		if (i == frames_count - 1) {// Если он последний
			int res = packet_size - curr_begin;
			frame2.setData(packet, curr_begin, res);
		}
		else {                      // Если он не последний
			frame2.setData(packet, curr_begin, packet_size);
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
			for (int i = 0; i < packet_size; i++) {
				out[prev_begin + i] = conv_frames[i];
			}
		}

		// Если пакет первый и в первом пакете, отрезать хвост
		if ((i == 1) && (packet_index == 1)) {

		}

		// Записать последний кадр в выходной массив
		if (i == frames_count - 1) {
			int res = 0;
			if (tails) {
				// сколько элементов требуется записать
				res = length - curr_begin;
				for (int i = 0; i < res + add_min; i++) {
					out[curr_begin + i] = conv_frames[packet_size + i];
				}
			}
			else {
				res = length - curr_begin;
				for (int i = 0; i < res; i++) {
					out[curr_begin + i] = conv_frames[packet_size + i];
				}
			}
		}
	}
	return 1;
}

int FLT_FilterPkt::filtrate1Packet(double* packet)
{
	//int frames_count = packet_size % frame_size != 0 ? floor(double(packet_size) / frame_size) + 1 : floor(double(packet_size) / frame_size);

	//int prev_begin = -frame_size;
	//int curr_begin = 0;

	//int begFrame = 0;	// Текущая позиция в frame1
	//int begPacket = 0;	// Текущая позиция в Packet1

	//// Почему [prev_begin, curr_begin] и [begFrame, begPacket]
	//// Логика в начале старая, а разбираться лень

	//for (int i = 1; i < frames_count; i++) {
	//	prev_begin = prev_begin + packet_size;
	//	curr_begin = curr_begin + packet_size;

	//	// Загрузить, обработать frame1
	//	if (i == 1) {
	//		// --- Если кадр первый
	//		frame1.setData(packet, prev_begin, packet_size);
	//		fft_filtrate(frame1);
	//	}
	//	else {  
	//		// --- Если кадр не первый
	//		// Сменить полусвёрнутый frame2 на frame1
	//		Frame::switchData(frame1, frame2);
	//	}

	//	// Обработать текущий кадр (CURRENT_FRAME)
	//	if (i == frames_count - 1) {// Если он последний
	//		int res = packet_size - curr_begin;
	//		frame2.setData(packet, curr_begin, res);
	//	}
	//	else {                      // Если он не последний
	//		frame2.setData(packet, curr_begin, packet_size);
	//	}

	int frames_count = packet_size % frame_size != 0 ? floor(double(packet_size) / frame_size) + 1 : floor(double(packet_size) / frame_size);

	int begFrame = 0;	// Текущая позиция в frame1
	int begPacket = 0;	// Текущая позиция в Packet1

	for (int i = 1; i < frames_count; i++) {
		// Загрузить, обработать frame1
		if (i == 1) {
			// --- Если кадр первый
			frame1.setData(packet, 0, packet_size);
			fft_filtrate(frame1);
		}
		else {
			// --- Если кадр не первый
			// Сменить полусвёрнутый frame2 на frame1
			Frame::switchData(frame1, frame2);
		}

		// Сколько элементов осталось записать в Packet1 из frame2 в последнем кадре
		int res = packet_size - (frames_count - 1) * frame_size;
		// Сколько вообще лишних элементов в frame2 (тк. (res + add_min), может быть != fft_size)
		int other = fft_size - res - add_min;

		// Обработать текущий кадр frame2
		if (i == frames_count - 1) {// Если он последний
			frame2.setData(packet, begPacket, res);
		}
		else {                      // Если он не последний
			frame2.setData(packet, begPacket, packet_size);
		}

		fft_filtrate(frame2);
		convolFull(frame1, frame2);

		// ------------------ Свёртка кадров
		// Если кадр первый
		if (i == 1) { 
			// Скопировать данные в frame1.data и фильтровать
			frame1.setData(packet, 0, packet_size);
			fft_filtrate(frame1);

			/*
				Свёртка frame1 и frame2, и запись кадра frame1 без левого хвоста
				в Packet1, frame2 выйдет со свёрнутым левым хвостом,
				далее его логическое начало должно быть сдвинуто на add_min2.
			*/

			// Переписываем от конца левого хвоста (обрезали) до участка суммирования
			begFrame = add_min2;
			for (int i = 0; i < frame_size - add_min2; i++)
				packet[begPacket + i] = frame1.data[begFrame + i];
			begFrame += frame_size - add_min2;
			begPacket += begFrame;
			// ---------- Свёртка от начала участка суммирования до конца участка суммирования.
			// --- на frame1
			for (int i = 0; i < add_min2; i++)
				packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[i];
			begFrame += add_min2;
			begPacket += begFrame;
			// --- на frame2
			for (int i = 0; i < add_min2; i++)
				frame2.data[add_min2 + i] = frame1.data[begFrame + i] + frame2.data[add_min2 + i];
		}
		else {
			// Если кадр промежуточный (не первый и не последний)
			if ((i != 1) && (i != (frames_count - 1))) { 
				// Если кадр второй, то его требуется загрузить и закончить конец свёртки с frame1
				if (i == 2) {

				}
				// Сменить полусвёрнутый frame2 на frame1
				Frame::switchData(frame1, frame2);
				// Скопировать данные в frame2.data и фильтровать
				frame2.setData(packet, begPacket, packet_size);
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
				begPacket += begFrame;
				// ---------- Свёртка от начала участка суммирования до конца участка суммирования.
				// --- на frame1
				for (int i = 0; i < add_min2; i++)
					packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[i];
				begFrame += add_min2;
				begPacket += begFrame;
				// --- на frame2
				for (int i = 0; i < add_min2; i++)
					frame2.data[add_min2 + i] = frame1.data[begFrame + i] + frame2.data[add_min2 + i];
			}
			// Если кадр последний
			else {	
				frame2.setData(packet, begPacket, res);
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
				begPacket += begFrame;
				// ---------- Свёртка от начала участка суммирования до конца участка суммирования.
				// --- на frame1
				
				// Сколько нужно отступить слева в frame2
				int begFrame2 = other / 2; 

				for (int i = 0; i < add_min2; i++)
					packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
				begFrame += add_min2;
				begFrame2 += add_min2;
				begPacket += begFrame;

				// --- на frame2
				// Свертка сигнала frame2 с хвостом frame1
				for (int i = 0; i < add_min2; i++)
					packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
				begFrame += add_min2;
				begFrame2 += add_min2;
				begPacket += begFrame;

				// Запись сигнала frame2 без хвоста
				for (int i = 0; i < res - add_min2; i++)
					packet[begPacket + i] = frame2.data[begFrame2 + i];
			}
		}

	return 1;
}

bool FLT_FilterPkt::startTransfer10(int length, int fft_size)
{
	packet_size = length;
	if (packet1 != nullptr)		delete[] packet1;
	if (packet2 != nullptr)		delete[] packet2;

	this->fft_size = fft_size;
}

bool FLT_FilterPkt::filtratePkt10(double* packet)
{
	packet_index++;
	// load packet
	
	int frames_count = packet_size % frame_size != 0
		? floor(double(packet_size) / frame_size) + 1
		: floor(double(packet_size) / frame_size);

	if (packet_index == 1) {

	}

	return false;
}

int FLT_FilterPkt::stopTransfer10(double* packet)
{
	return 0;
}

