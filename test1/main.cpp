#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include "fftw3.h"
#include "FLT_FilterFile.h"
#include "FLT_FilterPkt.h"
using namespace std;

double* generate_pulse(int length, double basicF, double fd, int harmonic_count) {
	double* signal = new double[length];
	for (int i = 0; i < length; i++) {
		double s = 0;
		for (int j = 1; j <= (harmonic_count + harmonic_count - 1); j += 2) {
			s = s + 1 / double(j) * sin(((2 * FILTER_PI * basicF * j) / double(fd)) * i);
		}
		signal[i] = s;
	}
	return signal;
}

void writeToFile(int number, double* arr, int length) {
	string str1 = "C:/Users/Даниил Левин/Documents/MATLAB/FLT/graph";
	string str2 = std::to_string(number);
	string str3 = ".txt";
	string strPath = str1 + str2 + str3;
	const char* chrPath = strPath.c_str();
	//printf_s("strPath: %s\n", strPath);
	fstream file;
	file.open(chrPath, std::fstream::in | std::fstream::out | std::fstream::trunc);
	if (!file.is_open()) {
		printf_s("Error opening file\n");
		exit(1);
	}
	for (int i = 0; i < length; i++) {
		file << arr[i] << std::endl;
	}
	file.close();
	//delete[] chrPath;
}

void writeToFile(int number, fftw_complex* arr, int length) {
	string str1 = "C:/Users/Даниил Левин/Documents/MATLAB/test4/graph";
	string str2 = std::to_string(number);
	string str3 = "c.txt";
	string strPath = str1 + str2 + str3;
	const char* chrPath = strPath.c_str();
	//printf_s("strPath: %s\n", strPath);
	fstream file;
	file.open(chrPath, std::fstream::in | std::fstream::out | std::fstream::trunc);
	if (!file.is_open()) {
		printf_s("Error opening file\n");
		exit(1);
	}
	for (int i = 0; i < length; i++) {
		file << arr[i][0] << ", " << arr[i][1] << std::endl;
	}
	file.close();
	//delete[] chrPath;
}

int main() {

	// ----------------------------------------------------------------------------------------------------
	//int N = 257;
	//double fd = 44100;
	//int accurancy = 1;
	//int window = 1;
	//
	//int length_in = 30'000;
	//int baseFreq = 30;
	//int harmonicsCount = 15;

	//double* signal_in = generate_pulse(length_in, baseFreq, fd, harmonicsCount);

	//double B1 = 500;

	////FLT_BaseFilter filter_base;
	////if (!filter_base.setIrLowpassR1B1(N, fd, B1, window)) {
	////	printf("Error in set type: %d", filter_base.get_error_code());
	////	exit(-1);
	////}
	//FLT_FilterFile filter_file;
	//if (!filter_file.setIrLowpassR1B1(N, fd, B1, window)) {
	//	printf("Error in set type: %d", filter_file.get_error_code());
	//	exit(-1);
	//}

	//int ft_size = 0;

	////double* h = nullptr;
	////filter_file.get_h(h);
	////writeToFile(1, h, N);

	////double* ph = nullptr;
	////ft_size = filter_file.get_h_phase(ph, 0);
	////writeToFile(2, ph, ft_size);

	////double* mag = nullptr;
	////ft_size = filter_file.get_h_magnitude(mag, 0);
	////writeToFile(3, mag, ft_size);

	////double* att = nullptr;
	////ft_size = filter_file.get_h_attenuation(att, 0);
	////writeToFile(4, att, ft_size);

	////delete[] h;
	////delete[] ph;
	////delete[] mag;
	////delete[] att;
	//int length_out = 0;
	//double* signal_out = nullptr;
	////length_out = filter_base.filtrate(signal_in, length_in, signal_out, accurancy, 0);
	////length_out = filter_file.filtrate(signal_in, length_in, signal_out, accurancy, 0);
	//length_out = filter_file.filtrateBlock(signal_in, length_in, signal_out, accurancy, 1);

	//printf("N           | %d\n", filter_file.get_N());
	//printf("accurancy   | %d\n", accurancy);
	//printf("sample_size | %d\n", filter_file.get_sample_size());
	//printf("fft_size:   | %d\n", filter_file.get_fft_size());
	//printf("fd          | %-10.3f kHz\n", fd / 1000);
	//printf("Window      | %d\n", filter_file.get_window());

	//printf("\n\n");
	//printf("Signal IN\n");
	//printf("\tLength          | %d\n", length_in);
	//printf("\tBase freq       | %d Hz:\n", baseFreq);
	//printf("\tHarmonics count | %d\n", harmonicsCount);
	//printf("\tMax frequency   | %d Hz:\n", baseFreq * harmonicsCount * harmonicsCount);

	//printf("\n");
	//printf("Signal OUT\n");
	//printf("\tLength          | %d\n", length_out);

	//writeToFile(5, signal_in, length_in);
	//writeToFile(6, signal_out, length_out);

	//delete[] signal_in;
	//delete[] signal_out;

	// ----------------------------------------------------------------------------------------------------

	int N = 257;
	double fd = 44100;
	int accurancy = 1;
	int window = 1;
	int B1 = 500;
	
	int integer_frames_count = 5;
	int frame_length = 1777;
	int additional_length = 0;
	int length_in = integer_frames_count * frame_length + additional_length;
	int frames_count = integer_frames_count + (additional_length > 0 ? 1 : 0);

	int baseFreq = 30;
	int harmonicsCount = 15;
	double* signal_in = nullptr;
	signal_in = generate_pulse(length_in, baseFreq, fd, harmonicsCount);

	FLT_FilterPkt filter_pkt;
	if (!filter_pkt.setIrLowpassR1B1(N, fd, B1, window)) {
		printf("Error in set type: %d", filter_pkt.get_error_code());
		exit(-1);
	}

	if (!filter_pkt.startTransfer(frame_length, accurancy))
		printf("ERROR StartTransfer\n");

	double* packet = new double[frame_length];
	int length_out = frame_length * frames_count;
	double* signal_out = new double[length_out];
	int index = 0;

	printf("frames_count: %d\n", frames_count);

	for (int i = 0; i < frames_count; i++)
	{
		printf("frame [%d] --------------------------- \n", i);
		for (int j = 0; j < frame_length; j++)
		{
			packet[j] = signal_in[i * frame_length + j];
		}

		bool success = filter_pkt.filtratePkt1(packet);
		if (!success) {
			if (filter_pkt.get_error_code() == FILTER_FIRST_PKT) {
				continue;
			}
			else {
				printf("First packet Error");
			}
		}
		else {
			for (int j = 0; j < frame_length; j++)
			{
				signal_out[index] = packet[j];
				index++;
				printf("%f\n", packet[j]);
			}
		}

	}

	double* lastPkt = nullptr;
	int lastLength = filter_pkt.stopTransfer(lastPkt);
	for (int i = 0; i < lastLength; i++)
	{
		signal_out[index] = lastPkt[i];
		index++;
	}

	writeToFile(5, signal_in, length_in);
	writeToFile(6, signal_out, length_out);

	delete[] lastPkt;

	delete[] signal_in;
	delete[] packet;
	delete[] signal_out;

	return 0;
}