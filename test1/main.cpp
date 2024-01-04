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
		printf("Error opening file\n");
		exit(1);
	}
	for (int i = 0; i < length; i++) {
		file << arr[i] << std::endl;
	}
	file.close();
	//delete[] chrPath;
}

int main() {
	setlocale(LC_ALL, "ru");

	//int N = 587;
	//double fd = 44100;
	//int accurancy = 1;
	//int window = 1;

	//int length_in = 100;
	//int basefreq = 30;
	//int harmonicscount = 15;

	//double* signal_in = generate_pulse(length_in, basefreq, fd, harmonicscount);
	////writeToFile(5, signal_in, length_in);

	//double b1 = 3401;

	//FLT_BaseFilter filter_base;
	//if (!filter_base.setIrLowpassR1B1(N, fd, b1, window)) {
	//	printf("error in set type: %d", filter_base.get_error_code());
	//	exit(-1);
	//}

	//int ft_size = 0;

	//int length_out = length_in + N - 1;
	////filter_base.filtrate(nullptr, length_in, 0);
	////double* signal_out = filter_base.filtrateT(signal_in, length_in, 0);
	//filter_base.filtrate(signal_in, length_in, 0);
	//filter_base.filtrate(signal_in, length_in, 0);
	//double* parameter = nullptr;
	//double* parameter2 = nullptr;
	////double* freq_match = nullptr;
	//bool symetrical = 1;
	//int size = 0;
	////filter_base.get_h_phase(parameter, symetrical);
	////size = filter_base.get_h_magnitude(parameter2, symetrical);
	////filter_base.get_freq_match(freq_match, symetrical);

	//double* mAttenuation = nullptr;
	//unsigned long time = 0;
	////int mSize = filter_base.measureAttenuation(mAttenuation, length_in, 1, 1, 1, fd / 2 - 1, time);

	////writeToFile(5, mAttenuation, mSize);
	////writeToFile(6, parameter, size);
	////writeToFile(6, parameter2, size);
	////writeToFile(6, signal_in, length_in);
	////writeToFile(6, signal_out, length_out);

	////delete[] mAttenuation;

	////delete[] parameter;
	////delete[] parameter2;
	////delete[] freq_match;
	//delete[] signal_in;
	////delete[] signal_out;

	// ----------------------------------------------------------------------------------------------------
	////int N = 257;
	//int N = 257;
	//double fd = 44100;
	//int accurancy = 1;
	//int window = 1;
	//
	//int length_in = 131072;
	//int baseFreq = 700;
	//int harmonicsCount = 4;

	//double* signal_in = generate_pulse(length_in, baseFreq, fd, harmonicsCount);
	//writeToFile(5, signal_in, length_in);

	//double B1 = 1401;

	//FLT_FilterFile filter_file;
	//if (!filter_file.setIrLowpassR1B1(N, fd, B1, window)) {
	//	printf("Error in set type: %d", filter_file.get_error_code());
	//	exit(-1);
	//}

	////FLT_FilterFile filter_file;
	////if (!filter_file.setIrBandstopR2B2(N, fd, 250, 300, 1400, 1450, window)) {
	////	printf("Error in set type: %d", filter_file.get_error_code());
	////	exit(-1);
	////}

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
	//using clock_t = std::chrono::high_resolution_clock; // пїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅ
	//std::chrono::time_point<clock_t> m_startTime;   // пїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅ
	//std::chrono::time_point<clock_t> m_endTime;     // пїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅ
	//std::chrono::nanoseconds m_durationTime_ns;     // пїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅ
	//unsigned long int m_int_durationTime_ns = 0;    // пїЅпїЅпїЅпїЅпїЅ [пїЅпїЅ]
	//m_startTime = clock_t::now();

	//writeToFile(5, signal_in, length_in);

	//int length_out = length_in + N - 1;
	//if (!filter_file.filtrateBlock(signal_in, length_in, 0, 0)) {
	//	printf("Error code: %d", filter_file.get_error_code());
	//	exit(-1);
	//}

	//m_endTime = clock_t::now();
	//m_durationTime_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(m_endTime - m_startTime);
	//m_int_durationTime_ns = m_durationTime_ns.count() /*- m_minimum_execution_time_timer_code*/; // пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅ "пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ"
	//printf("General time: %f\n", m_int_durationTime_ns / double(length_in));

	//if (!filter_file.filtrateBlock(signal_in, length_in, 0, 0)) {
	//	printf("Error code: %d", filter_file.get_error_code());
	//	exit(-1);
	//}

	//writeToFile(6, signal_in, length_in);
	////double* signal_out = filter_file.filtrateBlock(signal_in, length_in, 0);

	//printf("N           | %d\n", filter_file.get_N());
	//printf("accurancy   | %d\n", accurancy);
	////printf("sample_size | %d\n", filter_file.get_sample_size());
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

	////double* magnitude = nullptr;
	////int size = filter_file.get_h_attenuation(magnitude, 0);
	////writeToFile(6, magnitude, size);

	//double* mAttenuation = nullptr;
	//unsigned long time = 0;
	////int mSize = filter_file.measureAttenuation(mAttenuation, length_in,0, 1, 1, fd / 2, time);
	////printf("Time: %lu ns\n", time);

	////writeToFile(5, mAttenuation, mSize);
	////writeToFile(5, signal_in, length_in);
	////writeToFile(6, signal_out, length_out);

	////delete[] signal_out;
	//delete[] signal_in;
	////delete[] mAttenuation;

	// ----------------------------------------------------------------------------------------------------

	int N = 257;
	double fd = 44100;
	int window = 1;
	//int B1 = 120 + 120*3 + 120*5 + 120*7 + 120*9;
	int B1 = 1401;
	int baseFreq = 700;
	int harmonicsCount = 4;

	//int packet_count = 8;
	//int packet_length = 32768;
	//int signal_length = packet_length * packet_count;

	int signal_length = 131072*8;
	int packet_count = 4;
	int packet_length = signal_length / packet_count;

	printf("\n-------- Creating arrays             --------\n");

	double* signal_in = generate_pulse(signal_length, baseFreq, fd, harmonicsCount);
	double* signal_out = new double[signal_length];
	double** master_packet = new double* [packet_count];

	printf("\n\n======= Filter parameters          =======\n\n");
	printf("N:               %d\n", N);
	printf("Fd:              %.2f Hz\n", fd);
	printf("window:          %d\n", window);
	printf("Band:            %d Hz\n", B1);

	printf("\n\n======= Generate Signal parameters =======\n\n");
	printf("Signal length:   %d\n", signal_length);
	printf("Base frequency:  %d Hz\n", baseFreq);
	printf("Harmonics count: %d\n", harmonicsCount);
	printf("Max frequency:   %d Hz\n", baseFreq* harmonicsCount * 2);

	printf("\n\n======= Signal transfer parameters =======\n\n");
	printf("Signal length: %d\n", signal_length);
	printf("Packets count: %d\n", packet_count);
	printf("Packet length: %d\n", packet_length);

	printf("\n-------- Arrays created              --------\n");
	printf("\n-------- Setting impulse response    --------\n");

	FLT_FilterPkt filter_pkt;

	if (!filter_pkt.setIrBandstopR1B1(N, fd, B1, B1*3, window)) {
	//if (!filter_pkt.setIrHighpassR2B2(N, fd, 200, 300, window)) {
		printf("Error in set type: %d", filter_pkt.get_error_code());
		exit(-1);
	}

	printf("\n-------- Impulse response set        --------\n");
	printf("\n-------- Starting transfer           --------\n");

	if (!filter_pkt.startTransferBlock(packet_length, 0, 0))
	//if (!filter_pkt.startTransfer(packet_length, 0))
	{
		printf("Error in start transfer, error code: %d", filter_pkt.get_error_code());
		exit(-1);
	}

	printf("\n-------- Transfer started            --------\n");
	printf("\n-------- Creating packets            --------\n");
	// пїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅ
	int begin = 0;
	for (int j = 0; j < packet_count; j++)
	{
		master_packet[j] = new double[packet_length];
		for (int i = 0; i < packet_length; i++)
			master_packet[j][i] = signal_in[begin + i];
		begin += packet_length;
	}
	
	printf("\n-------- Packets created             --------\n");
	//printf("\n--------                             --------\n");

	printf("\n-------- Starting filtration         --------\n");
	// пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅ пїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅ, пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ
	begin = 0;
	double* packet;

	using clock_t = std::chrono::high_resolution_clock; // пїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅ
	std::chrono::time_point<clock_t> m_startTime;   // пїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅ
	std::chrono::time_point<clock_t> m_endTime;     // пїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅ
	std::chrono::nanoseconds m_durationTime_ns;     // пїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅ
	unsigned long int m_int_durationTime_ns = 0;    // пїЅпїЅпїЅпїЅпїЅ [пїЅпїЅ]
	m_startTime = clock_t::now();

	for (int i = 0; i < packet_count; i++)
	{
		packet = master_packet[i];
		if (!filter_pkt.filtratePktBlock(packet)) {
		//if (!filter_pkt.filtratePkt(packet)){
			// пїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅ
			continue;
		}
		else {
			// пїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅ пїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅ, пїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅ пїЅпїЅпїЅпїЅпїЅ
			//for (int i = 0; i < packet_length; i++)
			//	signal_out[begin + i] = packet[i];
			//begin += packet_length;
		}
	}

	printf("\n-------- Filtration completed        --------\n");
	printf("\n-------- Getting last packet         --------\n");

	// пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅ
	double* lastPacket = nullptr;
	lastPacket = filter_pkt.getLatestPktBlock();
	//lastPacket = filter_pkt.getLatestPkt();
	//for (int i = 0; i < packet_length; i++)
		//signal_out[begin + i] = lastPacket[i];
	
	m_endTime = clock_t::now();
	m_durationTime_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(m_endTime - m_startTime);
	m_int_durationTime_ns = m_durationTime_ns.count() /*- m_minimum_execution_time_timer_code*/; // пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅ "пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ"
	printf("General time: %f\n", m_int_durationTime_ns / double(signal_length));

	printf("\n-------- Last packet got             --------\n");
	printf("\n-------- Stopping transfer           --------\n");

	//пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ
	//filter_pkt.stopTransfer();
	filter_pkt.stopTransferBlock();
	

	printf("\n-------- Transfer stopped            --------\n");

	//printf("\n-------- Writting signals to files   --------\n");
	//// пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅ пїЅпїЅпїЅпїЅ
	//writeToFile(5, signal_in, signal_length);
	//writeToFile(6, signal_out, signal_length);
	//printf("\n-------- Signals wrotted             --------\n");

	printf("\n-------- Starting free               --------\n");

	delete[] lastPacket;

	if (!filter_pkt.startTransferBlock(packet_length, 0, 0)){
		printf("Error in start transfer, error code: %d", filter_pkt.get_error_code());
		exit(-1);
	}

	begin = 0;
	for (int i = 0; i < packet_count; i++)
	{
		packet = master_packet[i];
		if (!filter_pkt.filtratePktBlock(packet)) {
			continue;
		}
		else {
			for (int i = 0; i < packet_length; i++)
				signal_out[begin + i] = packet[i];
			begin += packet_length;
		}
	}

	lastPacket = filter_pkt.getLatestPktBlock();
	for (int i = 0; i < packet_length; i++)
		signal_out[begin + i] = lastPacket[i];

	filter_pkt.stopTransferBlock();

	printf("\n-------- Writting signals to files   --------\n");
	writeToFile(5, signal_in, signal_length);
	writeToFile(6, signal_out, signal_length);
	printf("\n-------- Signals wrotted             --------\n");

	// пїЅпїЅпїЅпїЅпїЅпїЅ
	delete[] lastPacket;

	for (int i = 0; i < packet_count; i++)
		delete[] master_packet[i];
	delete[] master_packet;

	delete[] signal_in;
	delete[] signal_out;

	printf("\n-------- Free ended                  --------\n");

	return 0;
}

// пїЅпїЅпїЅпїЅпїЅпїЅпїЅ packet filter (пїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅ, пїЅпїЅпїЅпїЅпїЅпїЅпїЅ) ----------------------------------------------
//int N = 257;
//double fd = 44100;
//int accurancy = 1;
//int window = 1;
//int B1 = 500;
//
//int integer_frames_count = 5;
//int frame_length = 1777;
//int additional_length = 0;
//int length_in = integer_frames_count * frame_length + additional_length;
//int frames_count = integer_frames_count + (additional_length > 0 ? 1 : 0);
//
//int baseFreq = 30;
//int harmonicsCount = 15;
//double* signal_in = nullptr;
//signal_in = generate_pulse(length_in, baseFreq, fd, harmonicsCount);
//
//FLT_FilterPkt filter_pkt;
//if (!filter_pkt.setIrLowpassR1B1(N, fd, B1, window)) {
//	printf("Error in set type: %d", filter_pkt.get_error_code());
//	exit(-1);
//}
//
//if (!filter_pkt.startTransfer(frame_length, accurancy))
//printf("ERROR StartTransfer\n");
//
//double* packet = new double[frame_length];
//int length_out = frame_length * frames_count + N - 1;
//double* signal_out = new double[length_out];
//int index = 0;
//
//printf("frames_count: %d\n", frames_count);
//
//for (int i = 0; i < frames_count; i++)
//{
//	printf("frame [%d] --------------------------- \n", i);
//	for (int j = 0; j < frame_length; j++)
//	{
//		packet[j] = signal_in[i * frame_length + j];
//	}
//
//	bool success = filter_pkt.filtratePkt1(packet);
//	if (!success) {
//		if (filter_pkt.get_error_code() == FILTER_FIRST_PKT) {
//			continue;
//		}
//		else {
//			printf("First packet Error");
//		}
//	}
//	else {
//		for (int j = 0; j < frame_length; j++)
//		{
//			signal_out[index] = packet[j];
//			index++;
//			printf("%f\n", packet[j]);
//		}
//	}
//
//}
//
//double* lastPkt = nullptr;
//int lastLength = filter_pkt.stopTransfer(lastPkt);
//for (int i = 0; i < lastLength; i++)
//{
//	signal_out[index] = lastPkt[i];
//	index++;
//}
//
//writeToFile(5, signal_in, length_in);
//writeToFile(6, signal_out, length_out);
//
//delete[] lastPkt;
//
//delete[] signal_in;
//delete[] packet;
//delete[] signal_out;
// --------------------------------------------------------------------------------------------





