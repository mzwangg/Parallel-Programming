#include "config.h"
#include "groebnerBase.hpp"
#include "GroebnerPthread.hpp"
#include "GroebnerOpenMp.hpp"

int main()
{
#ifndef __ARM_NEON
    void(*funArr[])(Eliminater*&, Eliminatee*&, int) = { serialGroebner, avx2Groebner, pthreadGroebner, pthreadAVX2Groebner,openmpGroebner,openmpSIMDGroebner };
    string nameArr[] = {"serialGroebner", "avx2Groebner" ,"pthreadGroebner","pthreadAVX2Groebner","openmpGroebner","openmpSIMDGroebner"};
    timingAll(nameArr, funArr, 6);
#else
    void(*funArr[])(Eliminater*&, Eliminatee*&, int) = { serialGroebner, neonGroebner, pthreadNeonGroebner, openmpSIMDGroebner };
    string nameArr[] = { "serialGroebner", "neonGroebner" ,"pthreadNeonGroebner","openmpSIMDGroebner"};
    timingAll(nameArr, funArr, 4);
#endif
    return 0;
}