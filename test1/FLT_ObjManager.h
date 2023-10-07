#pragma once
#include <unordered_map>

typedef unsigned long FLTSOURCE;

class FLT_ObjManager
{
private:
	FLTSOURCE descriptor = 0;
	static FLTSOURCE next_descriptor;
public:
	static std::unordered_map<FLTSOURCE, FLT_ObjManager*> list;
	FLTSOURCE get_id();
	FLT_ObjManager& get_ref();
	static FLT_ObjManager* get_pointer(FLTSOURCE id);
};

