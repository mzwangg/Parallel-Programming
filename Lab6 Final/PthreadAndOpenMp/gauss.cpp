#include "config.h"
#include "GaussBase.hpp"
#include "GaussPthread.hpp"
#include "GaussOpenMp.hpp"

int main()
{
#ifndef __ARM_NEON
 	//���Բ�ͬ�����ڲ�ͬ�����ģ�µ�����
 	void(*funArr1[])(int, float*, int) = { serial, AVX2Aligned ,pthreadBlock,pthreadCyclic,pthreadSSE,pthreadCyclicSem, pthreadBlockCyclic,pthreadDynamic,\
 		pthreadNoSimd,openmpCyclicStatic,openmpBlockCyclicStatic, openmpBlockStatic, openmpDynamic, openmpBlockCyclicDynamic, openmpColumn,\
 		openmpGuided , openmpManualSIMD };
   string nameArr1[] = {"serial", "AVX2Aligned" ,"pthreadBlock","pthreadCyclic","pthreadSSE","pthreadCyclicSem","pthreadBlockCyclic","pthreadDynamic",\
 		"pthreadNoSimd","openmpCyclicStatic","openmpBlockCyclicStatic","openmpBlockStatic","openmpDynamic","openmpBlockCyclicDynamic","openmpColumn",\
 		"openmpGuided","openmpManualSIMD"};
 	timingAll(nameArr1, funArr1, 17, 8);

 	//���Բ�ͬ�����ڲ�ͬ�߳����µ�����
 	void(*funArr2[])(int, float*, int) = { pthreadCyclic,pthreadDynamic,openmpBlockCyclicStatic,openmpDynamic };
 	string nameArr2[] = { "pthreadCyclic", "pthreadDynamic" ,"openmpBlockCyclicStatic","openmpBlockCyclicDynamic" };
 	timingAllThreadNum(nameArr2, funArr2, 4);
#else
    //����ARMƽ̨�µ�����
    void(*funArr2[])(int, float*, int) = { serial,neonAlignedMLS,pthreadNeon,openmpBlockCyclicDynamic };
    string nameArr2[] = { "serial", "neonAlignedMLS" ,"pthreadNeon","openmpBlockCyclicDynamic" };
    timingAll(nameArr2, funArr2, 4, 8);
#endif
	return 0;
}