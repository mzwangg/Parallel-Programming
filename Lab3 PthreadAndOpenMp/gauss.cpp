#include "config.h"
#include "GaussBase.hpp"
#include "GaussPthread.hpp"
#include "GaussOpenMp.hpp"

int main()
{
#ifndef __ARM_NEON
 	//测试不同方法在不同问题规模下的性能
 	void(*funArr1[])(int, float*, int) = { serial, AVX2Aligned ,pthreadBlock,pthreadCyclic,pthreadSSE,pthreadCyclicSem, pthreadBlockCyclic,pthreadDynamic,\
 		pthreadNoSimd,openmpCyclicStatic,openmpBlockCyclicStatic, openmpBlockStatic, openmpDynamic, openmpBlockCyclicDynamic, openmpColumn,\
 		openmpGuided , openmpManualSIMD };
   string nameArr1[] = {"serial", "AVX2Aligned" ,"pthreadBlock","pthreadCyclic","pthreadSSE","pthreadCyclicSem","pthreadBlockCyclic","pthreadDynamic",\
 		"pthreadNoSimd","openmpCyclicStatic","openmpBlockCyclicStatic","openmpBlockStatic","openmpDynamic","openmpBlockCyclicDynamic","openmpColumn",\
 		"openmpGuided","openmpManualSIMD"};
 	timingAll(nameArr1, funArr1, 17, 8);

 	//测试不同方法在不同线程数下的性能
 	void(*funArr2[])(int, float*, int) = { pthreadCyclic,pthreadDynamic,openmpBlockCyclicStatic,openmpDynamic };
 	string nameArr2[] = { "pthreadCyclic", "pthreadDynamic" ,"openmpBlockCyclicStatic","openmpBlockCyclicDynamic" };
 	timingAllThreadNum(nameArr2, funArr2, 4);
#else
    //测试ARM平台下的性能
    void(*funArr2[])(int, float*, int) = { serial,neonAlignedMLS,pthreadNeon,openmpBlockCyclicDynamic };
    string nameArr2[] = { "serial", "neonAlignedMLS" ,"pthreadNeon","openmpBlockCyclicDynamic" };
    timingAll(nameArr2, funArr2, 4, 8);
#endif
	return 0;
}