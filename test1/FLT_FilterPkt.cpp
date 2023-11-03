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
//	while (fft_size < (packet_size + add_min)) // fft_size < (����� ������� + ������� ���������� ���������)
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
	while (fft_size < (packet_size + add_min)) // fft_size < (����� ������� + ������� ���������� ���������)
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
	// ���� ����� ������
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

		// �������� � �������� ������ frame1.data
		for (int i = 0; i < packet_size; i++)
		{
			packet[i] = conv_frames[i];
			//printf("%f\n", packet[i]);
		}

		// ����������� ����������� ������� ������ � ���������� ������
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

	//int begFrame = 0;	// ������� ������� � frame1
	//int begPacket = 0;	// ������� ������� � Packet1

	//// ������ [prev_begin, curr_begin] � [begFrame, begPacket]
	//// ������ � ������ ������, � ����������� ����

	//for (int i = 1; i < frames_count; i++) {
	//	prev_begin = prev_begin + packet_size;
	//	curr_begin = curr_begin + packet_size;

	//	// ���������, ���������� frame1
	//	if (i == 1) {
	//		// --- ���� ���� ������
	//		frame1.setData(packet, prev_begin, packet_size);
	//		fft_filtrate(frame1);
	//	}
	//	else {  
	//		// --- ���� ���� �� ������
	//		// ������� ������������ frame2 �� frame1
	//		Frame::switchData(frame1, frame2);
	//	}

	//	// ���������� ������� ���� (CURRENT_FRAME)
	//	if (i == frames_count - 1) {// ���� �� ���������
	//		int res = packet_size - curr_begin;
	//		frame2.setData(packet, curr_begin, res);
	//	}
	//	else {                      // ���� �� �� ���������
	//		frame2.setData(packet, curr_begin, packet_size);
	//	}

	int frames_count = packet_size % frame_size != 0 ? floor(double(packet_size) / frame_size) + 1 : floor(double(packet_size) / frame_size);

	int begFrame = 0;	// ������� ������� � frame1
	int begPacket = 0;	// ������� ������� � Packet1

	for (int i = 1; i < frames_count; i++) {

		// ==========================================================================================

		// ------------------ ������ � ���������� ������
		switch (i)
		{
		case 1:	// ���� ���� ������
			// ����������� ������ � frame1.data � �����������
			frame1.setData(packet, 0, frame_size);
			fft_filtrate(frame1);

			/*
				������ frame1 � frame2, � ������ ����� frame1 ��� ������ ������
				� Packet1, frame2 ������ �� �������� ����� �������,
				����� ��� ���������� ������ ������ ���� �������� �� add_min2.
			*/

			// ������������ �� ����� ������ ������ (��������) �� ������� ������������
			begFrame = add_min2;
			for (int i = 0; i < frame_size - add_min2; i++)
				packet[begPacket + i] = frame1.data[begFrame + i];
			begFrame += frame_size - add_min2;
			begPacket += frame_size - add_min2;
			break;

		case 2:	// ���� ���� ������, �� frame2 ��������� ��������� � ��������� ������ � frame1 (i == 1) � frame2 (i == 2)
			// ����������� ������ � frame2.data � �����������
			frame2.setData(packet, begPacket, frame_size);
			fft_filtrate(frame2);

			// �������� frame2, ����� ��������� ������ frame1 (i == 1) � frame2 (i == 2)
			// ---------- ������ �� ������ ������� ������������ �� ����� ������� ������������.
			// --- �� frame1
			for (int i = 0; i < add_min2; i++)
				packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[i];
			begFrame += add_min2;
			begPacket += add_min2;
			// --- �� frame2 
			for (int i = 0; i < add_min2; i++)
				// ���������� ������ frame2 �� ������ (�� ������� �� add_min2), �������� � ������� frame1
				frame2.data[add_min2 + i] = frame1.data[begFrame + i] + frame2.data[add_min2 + i];
			break;

		default:
			// ���� ���� ������������� (�� ������, �� ������ � �� ���������)
			if ((i != 1) && (i != (frames_count - 1)))
			{
				// ������� ������������ frame2 �� frame1
				Frame::switchData(frame1, frame2);
				// ����������� ������ � frame2.data � �����������
				frame2.setData(packet, begPacket, frame_size);
				fft_filtrate(frame2);

				/*
				������ frame1 � frame2, � ������ ����� frame1 � Packet1,
				frame2 ������ �� �������� ����� �������, ����� ��� ����������
				������ ������ ���� �������� �� add_min2.
				*/

				// ������������ frame1 (�� ������� ������ �� add_min2) �� ������ ������
				begFrame = add_min2;
				for (int i = 0; i < frame_size - add_min2; i++)
					packet[begPacket + i] = frame1.data[begFrame + i];
				begFrame += frame_size - add_min2;
				begPacket += frame_size - add_min2;
				// ---------- ������ �� ������ ������� ������������ �� ����� ������� ������������.
				// --- �� frame1
				for (int i = 0; i < add_min2; i++)
					packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[i];
				begFrame += add_min2;
				begPacket += add_min2;
				// --- �� frame2
				for (int i = 0; i < add_min2; i++)
					frame2.data[add_min2 + i] = frame1.data[begFrame + i] + frame2.data[add_min2 + i];
			}
			// ���� ���� ���������
			else {
				// ������� ��������� �������� �������� � Packet1 �� frame2 � ��������� �����
				int res = packet_size - (frames_count - 1) * frame_size;

				frame2.setData(packet, begPacket, res);
				fft_filtrate(frame2);

				/*
				������ frame1 � frame2, � ������ ����� frame1 � Packet1,
				frame2 ������ �� �������� ����� �������, ����� ��� ����������
				������ ������ ���� �������� �� add_min2.
				*/

				// ������������ frame1 (�� ������� ������ �� add_min2) �� ������ ������
				begFrame = add_min2;
				for (int i = 0; i < frame_size - add_min2; i++)
					packet[begPacket + i] = frame1.data[begFrame + i];
				begFrame += frame_size - add_min2;
				begPacket += frame_size - add_min2;
				// ---------- ������ �� ������ ������� ������������ �� ����� ������� ������������.
				// --- �� frame1

				// ������� ������ ������ ��������� � frame2 (��. (res + add_min), ����� ���� != fft_size)
				int other = fft_size - frame_size - add_min;
				// ������� ����� ��������� ����� � frame2
				int begFrame2 = other / 2;

				for (int i = 0; i < add_min2; i++)
					packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
				begFrame += add_min2;
				begFrame2 += add_min2;
				begPacket += add_min2;

				// --- �� frame2
				// ������� ������� frame2 � ������� frame1
				for (int i = 0; i < add_min2; i++)
					packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
				begFrame += add_min2;
				begFrame2 += add_min2;
				begPacket += add_min2;

				// ������ ������� frame2 ��� ������
				for (int i = 0; i < frame_size - add_min2; i++)
					packet[begPacket + i] = frame2.data[begFrame2 + i];
				begFrame2 += frame_size - add_min2;

				// ������ ������� ������ add_min2 � packet1_add_right 
				// (����� ���������, ��� add_min2 ��������� ����� �� ������� ������ 
				// ��������� �������� �� ��������� �������)
				for (int i = 0; i < add_min2; i++)
					packet1_add_right[i] = frame2.data[begFrame2 + i];
			}
			break;
		}


		// ==========================================================================================================

		//// ���������, ���������� frame1
		//if (i == 1) {
		//	// --- ���� ���� ������
		//	frame1.setData(packet, 0, packet_size);
		//	fft_filtrate(frame1);
		//}
		//else {
		//	// --- ���� ���� �� ������
		//	// ������� ������������ frame2 �� frame1
		//	Frame::switchData(frame1, frame2);
		//}

		//// ������� ��������� �������� �������� � Packet1 �� frame2 � ��������� �����
		//int res = packet_size - (frames_count - 1) * frame_size;
		//// ������� ������ ������ ��������� � frame2 (��. (res + add_min), ����� ���� != fft_size)
		//int other = fft_size - packet_size - add_min;

		//// ���������� ������� ���� frame2
		//if (i == frames_count - 1) {// ���� �� ���������
		//	frame2.setData(packet, begPacket, packet_size);
		//}
		//else {                      // ���� �� �� ���������
		//	frame2.setData(packet, begPacket, packet_size);
		//}

		//fft_filtrate(frame2);
		//convolFull(frame1, frame2);

		//// ==========================================================================================

		//// ------------------ ������ ������
		//// ���� ���� ������
		//if (i == 1) { 
		//	// ����������� ������ � frame1.data � �����������
		//	frame1.setData(packet, 0, packet_size);
		//	fft_filtrate(frame1);

		//	/*
		//		������ frame1 � frame2, � ������ ����� frame1 ��� ������ ������
		//		� Packet1, frame2 ������ �� �������� ����� �������,
		//		����� ��� ���������� ������ ������ ���� �������� �� add_min2.
		//	*/

		//	// ������������ �� ����� ������ ������ (��������) �� ������� ������������
		//	begFrame = add_min2;
		//	for (int i = 0; i < frame_size - add_min2; i++)
		//		packet[begPacket + i] = frame1.data[begFrame + i];
		//	begFrame += frame_size - add_min2;
		//	begPacket += begFrame;
		//}
		//else {
		//	// ���� ���� ������, �� frame2 ��������� ��������� � ��������� ����� ������ � frame1
		//	if (i == 2) 
		//	{
		//		// ����������� ������ � frame2.data � �����������
		//		frame2.setData(packet, begPacket, packet_size);
		//		fft_filtrate(frame2);

		//		// �������� frame2, ����� ��������� ������ frame1 (i == 1) � frame2 (i == 2)
		//		// ---------- ������ �� ������ ������� ������������ �� ����� ������� ������������.
		//		// --- �� frame1
		//		for (int i = 0; i < add_min2; i++)
		//			packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[i];
		//		begFrame += add_min2;
		//		begPacket += begFrame;
		//		// --- �� frame2 
		//		for (int i = 0; i < add_min2; i++)
		//			frame2.data[add_min2 + i] = frame1.data[begFrame + i] + frame2.data[add_min2 + i];
		//	}
		//	else
		//	{
		//		// ���� ���� ������������� (�� ������, �� ������ � �� ���������)
		//		if ((i != 1) && (i != (frames_count - 1)))
		//		{
		//			// ������� ������������ frame2 �� frame1
		//			Frame::switchData(frame1, frame2);
		//			// ����������� ������ � frame2.data � �����������
		//			frame2.setData(packet, begPacket, packet_size);
		//			fft_filtrate(frame2);

		//			/*
		//			������ frame1 � frame2, � ������ ����� frame1 � Packet1,
		//			frame2 ������ �� �������� ����� �������, ����� ��� ����������
		//			������ ������ ���� �������� �� add_min2.
		//			*/

		//			// ������������ frame1 (�� ������� ������ �� add_min2) �� ������ ������
		//			begFrame = add_min2;
		//			for (int i = 0; i < frame_size - add_min2; i++)
		//				packet[begPacket + i] = frame1.data[begFrame + i];
		//			begFrame += frame_size - add_min2;
		//			begPacket += begFrame;
		//			// ---------- ������ �� ������ ������� ������������ �� ����� ������� ������������.
		//			// --- �� frame1
		//			for (int i = 0; i < add_min2; i++)
		//				packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[i];
		//			begFrame += add_min2;
		//			begPacket += begFrame;
		//			// --- �� frame2
		//			for (int i = 0; i < add_min2; i++)
		//				frame2.data[add_min2 + i] = frame1.data[begFrame + i] + frame2.data[add_min2 + i];
		//		}
		//		// ���� ���� ���������
		//		else {
		//			frame2.setData(packet, begPacket, res);
		//			fft_filtrate(frame2);

		//			/*
		//			������ frame1 � frame2, � ������ ����� frame1 � Packet1,
		//			frame2 ������ �� �������� ����� �������, ����� ��� ����������
		//			������ ������ ���� �������� �� add_min2.
		//			*/

		//			// ������������ frame1 (�� ������� ������ �� add_min2) �� ������ ������
		//			begFrame = add_min2;
		//			for (int i = 0; i < frame_size - add_min2; i++)
		//				packet[begPacket + i] = frame1.data[begFrame + i];
		//			begFrame += frame_size - add_min2;
		//			begPacket += begFrame;
		//			// ---------- ������ �� ������ ������� ������������ �� ����� ������� ������������.
		//			// --- �� frame1

		//			// ������� ����� ��������� ����� � frame2
		//			int begFrame2 = other / 2;

		//			for (int i = 0; i < add_min2; i++)
		//				packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
		//			begFrame += add_min2;
		//			begFrame2 += add_min2;
		//			begPacket += begFrame;

		//			// --- �� frame2
		//			// ������� ������� frame2 � ������� frame1
		//			for (int i = 0; i < add_min2; i++)
		//				packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
		//			begFrame += add_min2;
		//			begFrame2 += add_min2;
		//			begPacket += begFrame;

		//			// ������ ������� frame2 ��� ������
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

	int begFrame = 0;	// ������� ������� � frame1
	int begPacket = 0;	// ������� ������� � Packet1

	for (int i = 1; i < frames_count; i++) {
		// ------------------ ������ � ���������� ������
		switch (i)
		{
		case 1:	// ���� ���� ������
			// ����������� ������ � frame1.data � �����������
			frame1.setData(packet, 0, frame_size);
			fft_filtrate(frame1);

			/*
				������ frame1 � frame2, � ������ ����� frame1 ��� ������ ������
				� Packet1, frame2 ������ �� �������� ����� �������,
				����� ��� ���������� ������ ������ ���� �������� �� add_min2.
			*/

			// ���������� ����� ����� � packet2_add_left
			for (int i = 0; i < add_min2; i++)
				packet2_add_left[i] = frame1.data[i];
			begFrame += add_min2;

			// ������������ �� ����� ������ ������ (��������) �� ������� ������������
			for (int i = 0; i < frame_size - add_min2; i++)
				packet[i] = frame1.data[begFrame + i];
			begFrame += frame_size - add_min2;
			begPacket = frame_size - add_min2;
			printf("\n-������ ���� ���������-\n");
			break;

		case 2:	// ���� ���� ������, �� frame2 ��������� ��������� � ��������� ������ � frame1 (i == 1) � frame2 (i == 2)
			// ����������� ������ � frame2.data � �����������
			printf("-����� ����-\n");
			std::cout << "data sett start" << std::endl;
			frame2.setData(packet, begPacket, frame_size);
			printf("-���� ������-\n");
			printf("-��������-\n");
			fft_filtrate(frame2);
			printf("-��������������-\n");
			// �������� frame2, ����� ��������� ������ frame1 (i == 1) � frame2 (i == 2)
			// ---------- ������ �� ������ ������� ������������ �� ����� ������� ������������.
			// --- �� frame1
			for (int i = 0; i < add_min2; i++)
				packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[i];
			begFrame += add_min2;
			begPacket += add_min2;
			// --- �� frame2 
			for (int i = 0; i < add_min2; i++)
				// ���������� ������ frame2 �� ������ ������ (�� ������� �� add_min2), �������� � ������� frame1
				frame2.data[add_min2 + i] = frame1.data[begFrame + i] + frame2.data[add_min2 + i];
			// beginFrame ����� �� ����� frame1
			printf("-������ ���� ���������-\n");
			break;

		default:
			// ���� ���� ������������� (�� ������, �� ������ � �� ���������)
			if ((i != 1) && (i != 2) && (i != (frames_count - 1)))
			{
				// ������� ������������ frame2 �� frame1
				Frame::switchData(frame1, frame2);
				// ����������� ������ � frame2.data � �����������
				frame2.setData(packet, begPacket, frame_size);
				fft_filtrate(frame2);

				/*
				������ frame1 � frame2, � ������ ����� frame1 � Packet1,
				frame2 ������ �� �������� ����� �������, ����� ��� ����������
				������ ������ ���� �������� �� add_min2.
				*/

				// ������������ frame1 (�� ������� ������ �� add_min2) �� ������ ������
				begFrame = add_min2;
				for (int i = 0; i < frame_size - add_min2; i++)
					packet[begPacket + i] = frame1.data[begFrame + i];
				begFrame += frame_size - add_min2;
				begPacket += frame_size - add_min2;
				// ---------- ������ �� ������ ������� ������������ �� ����� ������� ������������.
				// --- �� frame1
				for (int i = 0; i < add_min2; i++)
					packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[i];
				begFrame += add_min2;
				begPacket += add_min2;
				// --- �� frame2
				for (int i = 0; i < add_min2; i++)
					frame2.data[add_min2 + i] = frame1.data[begFrame + i] + frame2.data[add_min2 + i];
				// beginFrame ����� �� ����� frame1
				printf("-������������� ���� ���������-\n");
			}
			// ���� ���� ���������
			else {
				// ������� ��������� �������� �������� � Packet1 �� frame2 � ��������� �����
				int res = packet_size - (frames_count - 1) * frame_size;

				frame2.setData(packet, begPacket, res);
				fft_filtrate(frame2);

				/*
				������ frame1 � frame2, � ������ ����� frame1 � Packet1,
				frame2 ������ �� �������� ����� �������, ����� ��� ����������
				������ ������ ���� �������� �� add_min2.
				*/

				// ������������ frame1 (�� ������� ������ �� add_min2) �� ������ ������
				begFrame = add_min2;
				for (int i = 0; i < frame_size - add_min2; i++)
					packet[begPacket + i] = frame1.data[begFrame + i];
				begFrame += frame_size - add_min2;
				begPacket += frame_size - add_min2;
				// ---------- ������ �� ������ ������� ������������ �� ����� ������� ������������.
				// --- �� frame1

				// ������� ������ ������ ��������� � frame2 (��. (res + add_min), ����� ���� != fft_size)
				int other = fft_size - frame_size - add_min;
				// ������� ����� ��������� ����� � frame2
				int begFrame2 = other / 2;

				for (int i = 0; i < add_min2; i++)
					packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
				begFrame += add_min2;
				begFrame2 += add_min2;
				begPacket += add_min2;

				// --- �� frame2
				// ������� ������� frame2 � ������� frame1
				for (int i = 0; i < add_min2; i++)
					packet[begPacket + i] = frame1.data[begFrame + i] + frame2.data[begFrame2 + i];
				begFrame += add_min2;	// beginFrame ����� �� ����� frame1
				begFrame2 += add_min2;
				begPacket += add_min2;

				// ������ ������� frame2 ��� ������
				for (int i = 0; i < frame_size - add_min2; i++)
					packet[begPacket + i] = frame2.data[begFrame2 + i];
				begFrame2 += frame_size - add_min2;

				// ������ ������� ������ (add_min2) � packet2_add_right 
				// (����� ���������, ��� add_min2 ��������� ����� �� ������� ������ 
				// ��������� �������� �� ��������� �������)
				for (int i = 0; i < add_min2; i++)
					packet2_add_right[i] = frame2.data[begFrame2 + i];
				printf("-��������� ���� ���������-\n");
			}
			break;
		}
	}
}

bool FLT_FilterPkt::check_fft_size(int fft_size)
{
	if (fft_size <= 0) {
		return 0;  // ������� ��� ������������� ����� �� �������� �������� ������
	}
	// 0 - ������� ������
	bool result = ((fft_size & (fft_size - 1)) == 0);
	if (!result || ((fft_size / 2) < (N - 1))) { // ���� �� ������� ������, fft_size ������ ���� ������ N ������� � 2 ����
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
		((packet_size / 4) < (fft_size - add_min)) // ���� ����� ������, ��� 4 �����
		)
		return false;

	frame_size = fft_size - add_min;
	this->packet_size = packet_size;
	this->fft_size = fft_size;
	add_min = N - 1;
	add_min2 = add_min / 2;
	packet_index = 0;

	// ������������� � ���� ������ ---------------
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

	// ������������� � ������������ ������ ------
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
		// ������ ���� �� ����� ���������, ������� ��� �������� ���������
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
		// ��� ������ �����������, ��������� �� ��������
		int begPacket1 = packet_size - add_min;
		for (int i = 0; i < add_min2; i++)
		{
			// ����������� ������ ����� packet1 � ����� ������� packet2
			packet1[begPacket1 + i] = packet1[begPacket1 + i] + packet2_add_left[i];
			// ����������� ����� ����� packet2 � ������ ������� packet1
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
		// ������ ���� �� ����� ���������, ������� ��� �������� ���������
		packet1 = ptrToAllocatedData1;
		packet2 = ptrToAllocatedData2;
		memcpy(packet1, packet, sizeof(double) * packet_size);

		filtrateFirstPacket(packet1);
		return nullptr;
		break;

	case 2:
		memcpy(packet2, packet, sizeof(double) * packet_size);
		filtrateIntermediatePacket(packet2);
		// ��� ������ �����������, ��������� �� ��������
		for (int i = 0; i < add_min2; i++)
		{
			// ����������� ������ ����� packet1 � ����� ������� packet2
			packet1[begPacket1 + i] = packet1[begPacket1 + i] + packet2_add_left[i];
			// ����������� ����� ����� packet2 � ������ ������� packet1
			packet2[i] = packet2[i] + packet1_add_right[i];
		}
		return packet1;
		break;

	default:
		packet1 = packet2;
		packet2 = ptrToAllocatedData1;

		memcpy(packet2, packet, sizeof(double) * packet_size);
		filtrateIntermediatePacket(packet2);
		// ��� ������ �����������, ��������� �� ��������
		for (int i = 0; i < add_min2; i++)
		{
			// ����������� ������ ����� packet1 � ����� ������� packet2
			packet1[begPacket1 + i] = packet1[begPacket1 + i] + packet2_add_left[i];
			// ����������� ����� ����� packet2 � ������ ������� packet1
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

