#include "FLT_FilterPkt.h"
#include <iostream>
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

void FLT_FilterPkt::filtrateFirstPacket(double* packet)
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

		// ==========================================================================================

		// ------------------ Свёртка и фильтрация кадров
		switch (i)
		{
		case 1:	// Если кадр первый
			// Скопировать данные в frame1.data и фильтровать
			frame1.setData(packet, 0, frame_size);
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
			begPacket += frame_size - add_min2;
			break;

		case 2:	// Если кадр второй, то frame2 требуется загрузить и закончить свёртку с frame1 (i == 1) и frame2 (i == 2)
			// Скопировать данные в frame2.data и фильтровать
			frame2.setData(packet, begPacket, frame_size);
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
				// Записываем данные frame2 от начала (со сдвигом на add_min2), свёрнутые с хвостом frame1
				frame2.data[add_min2 + i] = frame1.data[begFrame + i] + frame2.data[add_min2 + i];
			break;

		default:
			// Если кадр промежуточный (не первый, не второй и не последний)
			if ((i != 1) && (i != (frames_count - 1)))
			{
				// Сменить полусвёрнутый frame2 на frame1
				Frame::switchData(frame1, frame2);
				// Скопировать данные в frame2.data и фильтровать
				frame2.setData(packet, begPacket, frame_size);
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
				// Сколько элементов осталось записать в Packet1 из frame2 в последнем кадре
				int res = packet_size - (frames_count - 1) * frame_size;

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
				begPacket += frame_size - add_min2;
				// ---------- Свёртка от начала участка суммирования до конца участка суммирования.
				// --- на frame1

				// Сколько вообще лишних элементов в frame2 (тк. (res + add_min), может быть != fft_size)
				int other = fft_size - frame_size - add_min;
				// Сколько нужно отступить слева в frame2
				int begFrame2 = other / 2;

				for (int i = 0; i < add_min2; i++)
					packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
				begFrame += add_min2;
				begFrame2 += add_min2;
				begPacket += add_min2;

				// --- на frame2
				// Свертка сигнала frame2 с хвостом frame1
				for (int i = 0; i < add_min2; i++)
					packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
				begFrame += add_min2;
				begFrame2 += add_min2;
				begPacket += add_min2;

				// Запись сигнала frame2 без хвоста
				for (int i = 0; i < frame_size - add_min2; i++)
					packet[begPacket + i] = frame2.data[begFrame2 + i];
				begFrame2 += frame_size - add_min2;

				// Запись правого хвоста add_min2 в packet1_add_right 
				// (нужно учитывать, что add_min2 элементов слева от правого хвоста 
				// требуется свернуть со следующим пакетом)
				for (int i = 0; i < add_min2; i++)
					packet1_add_right[i] = frame2.data[begFrame2 + i];
			}
			break;
		}


		// ==========================================================================================================

		//// Загрузить, обработать frame1
		//if (i == 1) {
		//	// --- Если кадр первый
		//	frame1.setData(packet, 0, packet_size);
		//	fft_filtrate(frame1);
		//}
		//else {
		//	// --- Если кадр не первый
		//	// Сменить полусвёрнутый frame2 на frame1
		//	Frame::switchData(frame1, frame2);
		//}

		//// Сколько элементов осталось записать в Packet1 из frame2 в последнем кадре
		//int res = packet_size - (frames_count - 1) * frame_size;
		//// Сколько вообще лишних элементов в frame2 (тк. (res + add_min), может быть != fft_size)
		//int other = fft_size - packet_size - add_min;

		//// Обработать текущий кадр frame2
		//if (i == frames_count - 1) {// Если он последний
		//	frame2.setData(packet, begPacket, packet_size);
		//}
		//else {                      // Если он не последний
		//	frame2.setData(packet, begPacket, packet_size);
		//}

		//fft_filtrate(frame2);
		//convolFull(frame1, frame2);

		//// ==========================================================================================

		//// ------------------ Свёртка кадров
		//// Если кадр первый
		//if (i == 1) { 
		//	// Скопировать данные в frame1.data и фильтровать
		//	frame1.setData(packet, 0, packet_size);
		//	fft_filtrate(frame1);

		//	/*
		//		Свёртка frame1 и frame2, и запись кадра frame1 без левого хвоста
		//		в Packet1, frame2 выйдет со свёрнутым левым хвостом,
		//		далее его логическое начало должно быть сдвинуто на add_min2.
		//	*/

		//	// Переписываем от конца левого хвоста (обрезали) до участка суммирования
		//	begFrame = add_min2;
		//	for (int i = 0; i < frame_size - add_min2; i++)
		//		packet[begPacket + i] = frame1.data[begFrame + i];
		//	begFrame += frame_size - add_min2;
		//	begPacket += begFrame;
		//}
		//else {
		//	// Если кадр второй, то frame2 требуется загрузить и закончить конец свёртки с frame1
		//	if (i == 2) 
		//	{
		//		// Скопировать данные в frame2.data и фильтровать
		//		frame2.setData(packet, begPacket, packet_size);
		//		fft_filtrate(frame2);

		//		// Получили frame2, можно закончить свёртку frame1 (i == 1) и frame2 (i == 2)
		//		// ---------- Свёртка от начала участка суммирования до конца участка суммирования.
		//		// --- на frame1
		//		for (int i = 0; i < add_min2; i++)
		//			packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[i];
		//		begFrame += add_min2;
		//		begPacket += begFrame;
		//		// --- на frame2 
		//		for (int i = 0; i < add_min2; i++)
		//			frame2.data[add_min2 + i] = frame1.data[begFrame + i] + frame2.data[add_min2 + i];
		//	}
		//	else
		//	{
		//		// Если кадр промежуточный (не первый, не второй и не последний)
		//		if ((i != 1) && (i != (frames_count - 1)))
		//		{
		//			// Сменить полусвёрнутый frame2 на frame1
		//			Frame::switchData(frame1, frame2);
		//			// Скопировать данные в frame2.data и фильтровать
		//			frame2.setData(packet, begPacket, packet_size);
		//			fft_filtrate(frame2);

		//			/*
		//			Свёртка frame1 и frame2, и запись кадра frame1 в Packet1,
		//			frame2 выйдет со свёрнутым левым хвостом, далее его логическое
		//			начало должно быть сдвинуто на add_min2.
		//			*/

		//			// Переписываем frame1 (со сдвигом начала на add_min2) до начала свёртки
		//			begFrame = add_min2;
		//			for (int i = 0; i < frame_size - add_min2; i++)
		//				packet[begPacket + i] = frame1.data[begFrame + i];
		//			begFrame += frame_size - add_min2;
		//			begPacket += begFrame;
		//			// ---------- Свёртка от начала участка суммирования до конца участка суммирования.
		//			// --- на frame1
		//			for (int i = 0; i < add_min2; i++)
		//				packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[i];
		//			begFrame += add_min2;
		//			begPacket += begFrame;
		//			// --- на frame2
		//			for (int i = 0; i < add_min2; i++)
		//				frame2.data[add_min2 + i] = frame1.data[begFrame + i] + frame2.data[add_min2 + i];
		//		}
		//		// Если кадр последний
		//		else {
		//			frame2.setData(packet, begPacket, res);
		//			fft_filtrate(frame2);

		//			/*
		//			Свёртка frame1 и frame2, и запись кадра frame1 в Packet1,
		//			frame2 выйдет со свёрнутым левым хвостом, далее его логическое
		//			начало должно быть сдвинуто на add_min2.
		//			*/

		//			// Переписываем frame1 (со сдвигом начала на add_min2) до начала свёртки
		//			begFrame = add_min2;
		//			for (int i = 0; i < frame_size - add_min2; i++)
		//				packet[begPacket + i] = frame1.data[begFrame + i];
		//			begFrame += frame_size - add_min2;
		//			begPacket += begFrame;
		//			// ---------- Свёртка от начала участка суммирования до конца участка суммирования.
		//			// --- на frame1

		//			// Сколько нужно отступить слева в frame2
		//			int begFrame2 = other / 2;

		//			for (int i = 0; i < add_min2; i++)
		//				packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
		//			begFrame += add_min2;
		//			begFrame2 += add_min2;
		//			begPacket += begFrame;

		//			// --- на frame2
		//			// Свертка сигнала frame2 с хвостом frame1
		//			for (int i = 0; i < add_min2; i++)
		//				packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
		//			begFrame += add_min2;
		//			begFrame2 += add_min2;
		//			begPacket += begFrame;

		//			// Запись сигнала frame2 без хвоста
		//			for (int i = 0; i < res - add_min2; i++)
		//				packet[begPacket + i] = frame2.data[begFrame2 + i];
		//		}
		//	}
		//}
		//// ==========================================================================================
	}
}

void FLT_FilterPkt::filtrateIntermediatePacket(double* packet)
{

	int frames_count = packet_size % frame_size != 0 ? floor(double(packet_size) / frame_size) + 1 : floor(double(packet_size) / frame_size);

	int begFrame = 0;	// Текущая позиция в frame1
	int begPacket = 0;	// Текущая позиция в Packet1

	for (int i = 1; i < frames_count; i++) {
		// ------------------ Свёртка и фильтрация кадров
		switch (i)
		{
		case 1:	// Если кадр первый
			// Скопировать данные в frame1.data и фильтровать
			frame1.setData(packet, 0, frame_size);
			fft_filtrate(frame1);

			/*
				Свёртка frame1 и frame2, и запись кадра frame1 без левого хвоста
				в Packet1, frame2 выйдет со свёрнутым левым хвостом,
				далее его логическое начало должно быть сдвинуто на add_min2.
			*/

			// Записываем левый хвост в packet2_add_left
			for (int i = 0; i < add_min2; i++)
				packet2_add_left[i] = frame1.data[i];
			begFrame += add_min2;

			// Переписываем от конца левого хвоста (обрезали) до участка суммирования
			for (int i = 0; i < frame_size - add_min2; i++)
				packet[i] = frame1.data[begFrame + i];
			begFrame += frame_size - add_min2;
			begPacket = frame_size - add_min2;
			printf("\n-Первый кадр обработан-\n");
			break;

		case 2:	// Если кадр второй, то frame2 требуется загрузить и закончить свёртку с frame1 (i == 1) и frame2 (i == 2)
			// Скопировать данные в frame2.data и фильтровать
			printf("-Задаю дату-\n");
			std::cout << "data sett start" << std::endl;
			frame2.setData(packet, begPacket, frame_size);
			printf("-Дата задана-\n");
			printf("-Фильтрую-\n");
			fft_filtrate(frame2);
			printf("-Профильтровано-\n");
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
			// beginFrame дошёл до конца frame1
			printf("-Второй кадр обработан-\n");
			break;

		default:
			// Если кадр промежуточный (не первый, не второй и не последний)
			if ((i != 1) && (i != 2) && (i != (frames_count - 1)))
			{
				// Сменить полусвёрнутый frame2 на frame1
				Frame::switchData(frame1, frame2);
				// Скопировать данные в frame2.data и фильтровать
				frame2.setData(packet, begPacket, frame_size);
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
				// beginFrame дошёл до конца frame1
				printf("-Промежуточный кадр обработан-\n");
			}
			// Если кадр последний
			else {
				// Сколько элементов осталось записать в Packet1 из frame2 в последнем кадре
				int res = packet_size - (frames_count - 1) * frame_size;

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
				begPacket += frame_size - add_min2;
				// ---------- Свёртка от начала участка суммирования до конца участка суммирования.
				// --- на frame1

				// Сколько вообще лишних элементов в frame2 (тк. (res + add_min), может быть != fft_size)
				int other = fft_size - frame_size - add_min;
				// Сколько нужно отступить слева в frame2
				int begFrame2 = other / 2;

				for (int i = 0; i < add_min2; i++)
					packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
				begFrame += add_min2;
				begFrame2 += add_min2;
				begPacket += add_min2;

				// --- на frame2
				// Свертка сигнала frame2 с хвостом frame1
				for (int i = 0; i < add_min2; i++)
					packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
				begFrame += add_min2;	// beginFrame дошёл до конца frame1
				begFrame2 += add_min2;
				begPacket += add_min2;

				// Запись сигнала frame2 без хвоста
				for (int i = 0; i < frame_size - add_min2; i++)
					packet[begPacket + i] = frame2.data[begFrame2 + i];
				begFrame2 += frame_size - add_min2;

				// Запись правого хвоста (add_min2) в packet2_add_right 
				// (нужно учитывать, что add_min2 элементов слева от правого хвоста 
				// требуется свернуть со следующим пакетом)
				for (int i = 0; i < add_min2; i++)
					packet2_add_right[i] = frame2.data[begFrame2 + i];
				printf("-Последний кадр обработан-\n");
			}
			break;
		}
	}
}

bool FLT_FilterPkt::check_fft_size(int fft_size)
{
	if (fft_size <= 0) {
		return 0;  // Нулевое или отрицательное число не является степенью двойки
	}
	// 0 - степень двойки
	bool result = ((fft_size & (fft_size - 1)) == 0);
	if (!result || ((fft_size / 2) < (N - 1))) { // если не степень двойки, fft_size должен быть больше N минимум в 2 раза
		error_code = FILTER_ERROR_FFT;
		return false;
	}
	return true;
}

bool FLT_FilterPkt::startTransferBlock(int packet_size, int fft_size)
{
	add_min = N - 1;
	add_min2 = add_min / 2;

	if (
		!(check_fft_size(fft_size)) ||
		((packet_size / 4) < (fft_size - add_min)) // Если длина меньше, чем 4 кадра
		)
		return false;

	frame_size = fft_size - add_min;
	this->packet_size = packet_size;
	this->fft_size = fft_size;
	add_min = N - 1;
	add_min2 = add_min / 2;
	packet_index = 0;

	// Освобождается в этом Классе ---------------
	if (packet1_add_right	!= nullptr)	delete[] packet1_add_right;
	if (packet2_add_left	!= nullptr)	delete[] packet2_add_left;
	if (packet2_add_right	!= nullptr)	delete[] packet2_add_right;
	if (ptrToAllocatedData1	!= nullptr)	delete[] ptrToAllocatedData1;
	if (ptrToAllocatedData2	!= nullptr)	delete[] ptrToAllocatedData2;

	packet1_add_right	= new double[add_min2];
	packet2_add_left	= new double[add_min2];
	packet2_add_right	= new double[add_min2];
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
		printf("=== %-3d) Try to get first packet\t/\t", packet_index);
		// Первый кадр не будет возвращен, поэтому его требуеся сохранить
		packet1 = ptrToAllocatedData1;
		packet2 = ptrToAllocatedData2;

		memcpy(packet1, packet, sizeof(double) * packet_size);

		filtrateFirstPacket(packet1);
		printf("first packet getted\n");
		return false;
	}
	else {
		printf("=== %-3d) Try to get packet      / ", packet_index);
		printf("TRY memcpy / ");
		memcpy(packet2, packet, sizeof(double) * packet_size);
		printf("OK / ");
		printf("TRY filtrate / ");
		filtrateIntermediatePacket(packet2);
		printf("OK / ");
		printf("TRY convol / ");
		// Оба пакета фильтрованы, требуется их свернуть
		int begPacket1 = packet_size - add_min;
		for (int i = 0; i < add_min2; i++)
		{
			// Сворачиваем правую часть packet1 с левым хвостом packet2
			packet1[begPacket1 + i] = packet1[begPacket1 + i] + packet2_add_left[i];
			// Сворачиваем левую часть packet2 с правым хвостом packet1
			packet2[i] = packet2[i] + packet1_add_right[i];
		}
		printf("OK / ");
		// 
		printf("TRY memcpy / ");
		memcpy(packet, packet1, sizeof(double) * packet_size);
		printf("OK / ");
		packet1 = packet2;
		packet2 = ptrToAllocatedData1;
		printf("packet getted\n");
		return true;
	}
}

double* const FLT_FilterPkt::filtratePktBlock2(double* const packet)
{
	packet_index++;
	int begPacket1 = packet_size - add_min;

	switch (packet_index)
	{
	case 1:
		// Первый кадр не будет возвращен, поэтому его требуеся сохранить
		packet1 = ptrToAllocatedData1;
		packet2 = ptrToAllocatedData2;
		memcpy(packet1, packet, sizeof(double) * packet_size);

		filtrateFirstPacket(packet1);
		return nullptr;
		break;

	case 2:
		memcpy(packet2, packet, sizeof(double) * packet_size);
		filtrateIntermediatePacket(packet2);
		// Оба пакета фильтрованы, требуется их свернуть
		for (int i = 0; i < add_min2; i++)
		{
			// Сворачиваем правую часть packet1 с левым хвостом packet2
			packet1[begPacket1 + i] = packet1[begPacket1 + i] + packet2_add_left[i];
			// Сворачиваем левую часть packet2 с правым хвостом packet1
			packet2[i] = packet2[i] + packet1_add_right[i];
		}
		return packet1;
		break;

	default:
		packet1 = packet2;
		packet2 = ptrToAllocatedData1;

		memcpy(packet2, packet, sizeof(double) * packet_size);
		filtrateIntermediatePacket(packet2);
		// Оба пакета фильтрованы, требуется их свернуть
		for (int i = 0; i < add_min2; i++)
		{
			// Сворачиваем правую часть packet1 с левым хвостом packet2
			packet1[begPacket1 + i] = packet1[begPacket1 + i] + packet2_add_left[i];
			// Сворачиваем левую часть packet2 с правым хвостом packet1
			packet2[i] = packet2[i] + packet1_add_right[i];
		}
		return packet1;

		break;
	}
}

double* FLT_FilterPkt::getLastPktBlock1()
{
	double* lastPacket = new double[packet_size];
	memcpy(lastPacket, packet2, sizeof(double) * packet_size);
	return lastPacket;
}

void FLT_FilterPkt::stopTransferBlock()
{
	packet_size = 0;
	frame_size = 0;
	packet_index = 0;
	add_min2 = 0;

	if (packet1_add_right != nullptr)	delete[] packet1_add_right;		packet1_add_right = nullptr;
	if (packet2_add_left != nullptr)	delete[] packet2_add_left;		packet2_add_left = nullptr;
	if (packet2_add_right != nullptr)	delete[] packet2_add_right;		packet2_add_right = nullptr;
	if (ptrToAllocatedData1 != nullptr)	delete[] ptrToAllocatedData1;	ptrToAllocatedData1 = nullptr;
	if (ptrToAllocatedData2 != nullptr)	delete[] ptrToAllocatedData2;	ptrToAllocatedData2 = nullptr;
}

