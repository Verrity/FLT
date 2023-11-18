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
	while (fft_size < (packet_size + add_min)) // fft_size < (����� ������� + ������� ���������� ���������)
	{
		fft_size = fft_size << 1;
	}
	fft_size << accurancy;
	frame_size = packet_size;

	if (isMinAllocated)
		free_min();
	init_min(N, fd, fft_size, packet_length);

	return true;
}

bool FLT_FilterPkt::filtratePkt1(double* packet)
{
	packet_index++;
	// ���� ����� ������
	switch (packet_index){
	case 1:
		frame1.setData(packet, 0, packet_size);
		fft_filtrate(frame1);
		error_code = FILTER_FIRST_PKT;
		return false;
		break;
	default:
		frame2.setData(packet, 0, packet_size);
		fft_filtrate(frame2);

		for (int i = 0; i < frame_size - add_min2; i++){
			packet[i] = frame1.data[add_min2 + i];
		}
		for (int i = 0; i < add_min2; i++){
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

bool FLT_FilterPkt::startTransferBlock(int packet_size, int accurancy)
{
	this->accurancy = 1;
	min_fft = 1;
	this->packet_size = 0;

	if (accurancy < 0)
		return false;
	// ����������� ������� �������� ��� � ����� N ������
	while (min_fft < N) {
		min_fft = min_fft << 1;
	}
	if (accurancy > 0)
		this->accurancy = accurancy;
	// ���������� ������ ���
	fft_size = min_fft << this->accurancy;
	add_min = N - 1;
	add_min2 = add_min / 2;

	if ((packet_size / 4) < (fft_size - add_min)){ // ���� ����� ������, ��� 4 �����
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

	// ������������� � ���� ������ ---------------
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

	// ������������� � ������������ ������ ------
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
		// ������ ���� �� ����� ���������, ������� ��� �������� ���������
		packet1 = ptrToAllocatedData1;
		packet2 = ptrToAllocatedData2;

		memcpy(packet1, packet, sizeof(double) * packet_size);
		_filtrateBlock(packet1, packet_size, nullptr, right_tail0);

		return false;
	}
	else {
		memcpy(packet2, packet, sizeof(double) * packet_size);
		_filtrateBlock(packet2, packet_size, left_tail, right_tail);

		// ��� ������ �����������, ��������� �� ��������
		int begPacket1 = packet_size - add_min2;
		for (int i = 0; i < add_min2; i++)
		{
			// ����������� ������ ����� packet1 � ����� ������� packet2
			packet1[begPacket1 + i] = packet1[begPacket1 + i] + left_tail[i];
			// ����������� ����� ����� packet2 � ������ ������� packet1
			packet2[i] = packet2[i] + right_tail0[i];
		}

		memcpy(packet, packet1, sizeof(double) * packet_size);
		
		double* temp0 = packet1;
		packet1 = packet2;
		packet2 = temp0;
		// ������ ����� ������ 2 ����� ��������� ������� ������ ������ 1, �.�. ����� 2 --> ����� 1
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
	return true;
}

