#pragma once
#include <unordered_map>
typedef unsigned long IR;

class FLT_ImpulseResponse
{
	
	int N = 0;
	double fd = 0;
	int window = 0;
	double* h = nullptr;
	int sample_size = 0;
	int fft_size = 0;
	double accurancy = 0;
	void(*pIRFunction) = nullptr;
	int errorCode = 0;

public:
	int setIRFunctionFile(int N, double fd, double B1, double accurancy, double fft_size, int window, double* h,
		void(*pIRFunction)	(int N, double fd, double B1, double accurancy, int window, double* h));

	FLT_ImpulseResponse();
	~FLT_ImpulseResponse();

	IR get_id();
	FLT_ImpulseResponse& get_ref();
	static FLT_ImpulseResponse* get_pointer(IR id);
private:
	IR descriptor = 0;
	static IR next_descriptor;
public:
	static std::unordered_map<IR, FLT_ImpulseResponse*> list;

};

