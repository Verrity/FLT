#include "Test.h"

std::unordered_map<Test::DESCRIPTOR, Test*> Test::list;
Test::DESCRIPTOR Test::next_descriptor;

Test::Test(int _value1, int _value2) : value1(_value1), value2(_value2)
{
	if (list.empty()) {
		descriptor = 0;
		next_descriptor = 1;
	}
	descriptor = next_descriptor;
	next_descriptor++;
	list[descriptor] = this;
}

Test::~Test()
{
}

Test::DESCRIPTOR Test::get_id()
{
	return descriptor;
}

Test& Test::get_ref()
{
	return *this;
}

Test* Test::get_pointer(DESCRIPTOR id)
{
	auto it = list.find(id);
	if (it != list.end()) {
		Test* ptr = it->second;
		return ptr;
	}
	else {
		return nullptr;
	}
}


