#pragma once
#include <unordered_map>
#include "fftw3.h"
#include <iostream>

typedef unsigned long TESTER;

class Test
{
public:
	int value1 = 0;
	int value2 = 0;

	typedef unsigned long DESCRIPTOR;
	Test(int _value1, int _value2);
	~Test();
	DESCRIPTOR get_id();
	Test& get_ref();
	static Test* get_pointer(DESCRIPTOR id);
private:
	DESCRIPTOR descriptor = 0;
	static DESCRIPTOR next_descriptor;
public:
	static std::unordered_map<DESCRIPTOR, Test*> list;
};

