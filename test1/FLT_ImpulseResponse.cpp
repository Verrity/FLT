#include "FLT_ImpulseResponse.h"
#include <time.h>

IR FLT_ImpulseResponse::next_descriptor = 1;
std::unordered_map<IR, FLT_ImpulseResponse*> FLT_ImpulseResponse::list;

int FLT_ImpulseResponse::setIRFunctionFile(int N, double fd, double B1, double accurancy, double fft_size, int window, double* h,
	void(*pIRFunction)	(int N, double fd, double B1, double accurancy, int window, double* h))
{
	// CHECKS --------
	this->pIRFunction = pIRFunction;
	this->h = h;
	this->fd = fd;
	this->N = N;
	this->accurancy = accurancy;
	this->window = window;
	this->fft_size = fft_size;
	this->accurancy = accurancy;

	sample_size = fft_size / accurancy;
	// TRY CATCH -------------------------------------------
	h = new double[N];
	pIRFunction(N, fd, B1, this->accurancy, window, this->h);
	if (window) {
		//useWindow --- 
	}

	// RETURN ERROR CODE
	return 0;
}
FLT_ImpulseResponse::FLT_ImpulseResponse()
{
	if (list.empty()) {
		descriptor = 0;
		next_descriptor = 1;
	}
	descriptor = next_descriptor;
	next_descriptor++;
	list[descriptor] = this;
}

FLT_ImpulseResponse::~FLT_ImpulseResponse()
{
	delete[] h;
}
IR FLT_ImpulseResponse::get_id()
{
	return descriptor;
}
FLT_ImpulseResponse& FLT_ImpulseResponse::get_ref()
{
	return *this;
}
FLT_ImpulseResponse* FLT_ImpulseResponse::get_pointer(IR id)
{
	auto it = FLT_ImpulseResponse::list.find(id);
	if (it != FLT_ImpulseResponse::list.end()) {
		FLT_ImpulseResponse* ptr = it->second;
		return ptr;
	}
	else {
		return nullptr;
	}
};