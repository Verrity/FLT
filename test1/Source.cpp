#include <iostream>
#include "FLT_filter.h"
#include "Test.h"

#define DEBUG_OBJ

//int FLT_CreateLowpassR1B1File(FILTER filter, int N, double fd, double BP, double BS, int window) {
//    //FLT_Filter* _filter = new FLT_Filter(filter, N, fd);
//    return 1;
//}
//
//void FLT_Free(FILTER filter) {
//    
//}

//=======================================================================================================


//int FLT_create(TESTER &filter, int N, double fd, double BP, double BS, int window) {
//    /*
//        Check parameters
//        if (...)
//        if (...)
//        else { return ERROR }
//    */
//    Test* pFilter = new Test(N, fd);
//    filter = pFilter->get_id();
//    #ifdef DEBUG_OBJ
//        std::cout << "[OBJ ACTION] id(" << pFilter->get_id() << ") address(" << pFilter << ") " << "Created" << std::endl;
//    #endif
//    return FLT_OK;
//}
//
//int FLT_Free(TESTER filter) {
//    auto it = Test::list.find(filter);
//    if (it != Test::list.end()) {
//#ifdef DEBUG_OBJ
//        Test* pFilter = it->second->get_pointer(filter);
//#endif
//        delete it->second;
//        Test::list.erase(filter);
//#ifdef DEBUG_OBJ
//        Test& ref = *pFilter;
//        std::cout << "[OBJ ACTION] id(" << filter << ") address(" << &ref << ") " << "Destroyed" << std::endl;
//#endif
//        return FLT_OK;
//    }
//    else {
//        return FLT_WRONG_PARAMETER;
//    }
//}
//
//void FLT_Free() {
//    for (const auto& pair : Test::list) {
//        Test::DESCRIPTOR descriptor = pair.first;
//        Test* pFilter = pair.second;
//        if (pFilter != nullptr) {
//#ifdef DEBUG_OBJ
//            Test& ref = *pFilter;
//            std::cout << "[OBJ ACTION] id(" << descriptor << ") address(" << &ref << ") " << "Destroyed" << std::endl;
//#endif
//            delete pair.second;
//        }
//    }
//    Test::list.clear();
//}

//=======================================================================================================

int FLT_CreateLowpassR1B1File(FILTER &filter, int N, double fd, double BP, double BS, int window) {
    /*
        Check parameters
        if (...)
        if (...)
        else { return ERROR }
    */
    FLT_Filter* pFilter = new FLT_Filter(filter, N, fd);
    filter = pFilter->get_id();
#ifdef DEBUG_OBJ
    std::cout << "[OBJ ACTION] id(" << pFilter->get_id() << ") address(" << pFilter << ") " << "Created" << std::endl;
#endif

    
    return FLT_OK;
}

int FLT_Free(FILTER filter) {
    auto it = FLT_Filter::list.find(filter);
    if (it != FLT_Filter::list.end()) {
#ifdef DEBUG_OBJ
        FLT_Filter* pFilter = it->second->get_pointer(filter);
#endif
        delete it->second;
        FLT_Filter::list.erase(filter);
#ifdef DEBUG_OBJ
        FLT_Filter& ref = *pFilter;
        std::cout << "[OBJ ACTION] id(" << filter << ") address(" << &ref << ") " << "Destroyed" << std::endl;
#endif
        return FLT_OK;
    }
    else {
        return FLT_WRONG_PARAMETER;
    }
}

void FLT_Free() {
    for (const auto& pair : FLT_Filter::list) {
        FILTER descriptor = pair.first;
        FLT_Filter* pFilter = pair.second;
        if (pFilter != nullptr) {
#ifdef DEBUG_OBJ
            FLT_Filter& ref = *pFilter;
            std::cout << "[OBJ ACTION] id(" << descriptor << ") address(" << &ref << ") " << "Destroyed" << std::endl;
#endif
            delete pair.second;
        }
    }
    FLT_Filter::list.clear();
}

int main() {
    setlocale(LC_ALL, "rus");

    FILTER F1, F2, F3;
    FLT_CreateLowpassR1B1File(F1, 512, 44100, 10'000, 15'000, 1);
    FLT_CreateLowpassR1B1File(F2, 512, 44100, 10'000, 15'000, 1);
    FLT_CreateLowpassR1B1File(F3, 512, 44100, 10'000, 15'000, 1);
    FLT_Free(F2);
    FLT_Free();

    //TESTER lool;
    //TESTER lool2;
    //TESTER lool3;
    //FLT_create(lool, 10, 20, 1, 1, 1);
    //FLT_create(lool2, 10, 20, 1, 1, 1);
    //FLT_create(lool3, 10, 20, 1, 1, 1);
    //FLT_Free(2);
    //FLT_Free();

    return 0;
}

