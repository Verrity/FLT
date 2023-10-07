#include "FLT_ObjManager.h"

FLTSOURCE FLT_ObjManager::next_descriptor = 1;
std::unordered_map<FLTSOURCE, FLT_ObjManager*> FLT_ObjManager::list;

FLT_ObjManager::FLT_ObjManager()
{
	if (list.empty()) {
		descriptor = 0;
		next_descriptor = 1;
	}
	descriptor = next_descriptor;
	next_descriptor++;
	list[descriptor] = this;
}

FLTSOURCE FLT_ObjManager::get_id()
{
	return descriptor;
}
FLT_ObjManager& FLT_ObjManager::get_ref()
{
	return *this;
}
FLT_ObjManager* FLT_ObjManager::get_pointer(FLTSOURCE id)
{
	auto it = FLT_ObjManager::list.find(id);
	if (it != FLT_ObjManager::list.end()) {
		FLT_ObjManager* ptr = it->second;
		return ptr;
	}
	else {
		return nullptr;
	}
}