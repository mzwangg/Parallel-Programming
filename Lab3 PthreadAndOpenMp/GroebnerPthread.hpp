#include "config.h"
#include "EliminateClass.hpp"

pthread_t threads[32];
groebnerPara thread_param_t[32];
pthread_mutex_t mutexI;
pthread_barrier_t startBarrier;
pthread_barrier_t finishBarrier;

Eliminater* eliminater;
Eliminatee* eliminatee;
int threadNum, jIndex, eeIndex;
extern int Nmat;//矩阵阶数
extern int Nee;//被消元行行数
extern int Ner;//消元子行数

void* subthreadGroebner(void* _params)
{
    groebnerPara* thread_param = (groebnerPara*)_params;
    int t_id = thread_param->t_id;
    uint* elier, * eliee;

    while (true)
    {
        if (t_id == 0)
        {
            eeIndex = 0;
            while(--jIndex>=0)
            { // 遍历消元子
                elier = eliminater->m_elierList[jIndex];
                if (elier == nullptr)
                { // 不存在对应消元子，则找出第一个被消元行升格
                    for (int i = 0; i < Nee; i++)
                    { // 遍历被消元行
                        if (eliminatee->m_isElierList[i])
                            continue;
                        eliee = eliminatee->m_elieeList[i];
                        if (getBit(eliee, jIndex))
                        {
                            eliminater->m_elierList[jIndex] = eliee;
                            eliminatee->m_isElierList[i] = true;
                            ++jIndex;
                            break;
                        }
                    }
                    continue;
                }
                break;
            }
        }

        pthread_barrier_wait(&startBarrier);
        if (jIndex < 0)
        {
            if (t_id != 0)pthread_exit(NULL);
            return NULL;
        }

        while (eeIndex < Nee)
        {
            //动态划分
            pthread_mutex_lock(&mutexI);
            while (eeIndex < Nee && (eliminatee->m_isElierList[eeIndex] || !getBit(eliminatee->m_elieeList[eeIndex], jIndex)))++eeIndex;
            int i = eeIndex++;
            pthread_mutex_unlock(&mutexI);
            if (i >= Nee) break;

            //进行消元
            eliee = eliminatee->m_elieeList[i];
            elier = eliminater->m_elierList[jIndex];
            int strLen = eliminatee->m_Ncol >> 5;
            for (int k = 0; k < strLen; ++k)
                eliee[k] ^= elier[k];
        }
        pthread_barrier_wait(&finishBarrier);
    }
}

void pthreadGroebner(Eliminater*& m_eliminater, Eliminatee*& m_eliminatee, int m_threadNum)
{
    // 初始化
    jIndex = Nmat;
    eliminater = m_eliminater;
    eliminatee = m_eliminatee;
    threadNum = m_threadNum;
    pthread_barrier_init(&startBarrier, NULL, threadNum);
    pthread_barrier_init(&finishBarrier, NULL, threadNum);
    pthread_mutex_init(&mutexI, NULL);

    for (int th = 0; th < threadNum; th++)
    {
        thread_param_t[th].t_id = th;
        if(th!=0)pthread_create(&threads[th], NULL, subthreadGroebner, (void*)&thread_param_t[th]);
    }
    subthreadGroebner((void*)&thread_param_t[0]);

    //销毁
    pthread_barrier_destroy(&startBarrier);
    pthread_barrier_destroy(&finishBarrier);
    pthread_mutex_destroy(&mutexI);
}

#ifndef __ARM_NEON
void* subthreadAVX2Groebner(void* _params)
{
    groebnerPara* thread_param = (groebnerPara*)_params;
    int t_id = thread_param->t_id;
    __m256i velier, veliee;
    uint* elier, * eliee;

    while (true)
    {
        if (t_id == 0)
        {
            eeIndex = 0;
            while (--jIndex >= 0)
            { // 遍历消元子
                elier = eliminater->m_elierList[jIndex];
                if (elier == nullptr)
                { // 不存在对应消元子，则找出第一个被消元行升格
                    for (int i = 0; i < Nee; i++)
                    { // 遍历被消元行
                        if (eliminatee->m_isElierList[i])
                            continue;
                        eliee = eliminatee->m_elieeList[i];
                        if (getBit(eliee, jIndex))
                        {
                            eliminater->m_elierList[jIndex] = eliee;
                            eliminatee->m_isElierList[i] = true;
                            ++jIndex;
                            break;
                        }
                    }
                    continue;
                }
                break;
            }
        }

        pthread_barrier_wait(&startBarrier);
        if (jIndex < 0)
        {
            if (t_id != 0)pthread_exit(NULL);
            return NULL;
        }

        while (eeIndex < Nee)
        {
            //动态划分
            pthread_mutex_lock(&mutexI);
            while (eeIndex < Nee && (eliminatee->m_isElierList[eeIndex] || !getBit(eliminatee->m_elieeList[eeIndex], jIndex)))++eeIndex;
            int i = eeIndex++;
            pthread_mutex_unlock(&mutexI);
            if (i >= Nee) break;

            //进行消元
            eliee = eliminatee->m_elieeList[i];
            elier = eliminater->m_elierList[jIndex];
            for (int k = (eliminatee->m_Ncol >> 5) - 8; k >= 0; k -= 8)
            {
                velier = _mm256_load_si256((__m256i*)(elier + k));
                veliee = _mm256_load_si256((__m256i*)(eliee + k));
                _mm256_store_si256((__m256i*)(eliee + k), _mm256_xor_si256(velier, veliee));
            }
        }
        pthread_barrier_wait(&finishBarrier);
    }
}

void pthreadAVX2Groebner(Eliminater*& m_eliminater, Eliminatee*& m_eliminatee, int m_threadNum)
{
    // 初始化
    jIndex = Nmat;
    eliminater = m_eliminater;
    eliminatee = m_eliminatee;
    threadNum = m_threadNum;
    pthread_barrier_init(&startBarrier, NULL, threadNum);
    pthread_barrier_init(&finishBarrier, NULL, threadNum);
    pthread_mutex_init(&mutexI, NULL);

    for (int th = 0; th < threadNum; th++)
    {
        thread_param_t[th].t_id = th;
        if (th != 0)pthread_create(&threads[th], NULL, subthreadAVX2Groebner, (void*)&thread_param_t[th]);
    }
    subthreadAVX2Groebner((void*)&thread_param_t[0]);

    //销毁
    pthread_barrier_destroy(&startBarrier);
    pthread_barrier_destroy(&finishBarrier);
    pthread_mutex_destroy(&mutexI);
}

#else
void* subthreadNeonGroebner(void* _params)
{
    groebnerPara* thread_param = (groebnerPara*)_params;
    int t_id = thread_param->t_id;
    uint32x4_t velier, veliee;
    uint* elier, * eliee;

    while (true)
    {
        if (t_id == 0)
        {
            eeIndex = 0;
            while (--jIndex >= 0)
            { // 遍历消元子
                elier = eliminater->m_elierList[jIndex];
                if (elier == nullptr)
                { // 不存在对应消元子，则找出第一个被消元行升格
                    for (int i = 0; i < Nee; i++)
                    { // 遍历被消元行
                        if (eliminatee->m_isElierList[i])
                            continue;
                        eliee = eliminatee->m_elieeList[i];
                        if (getBit(eliee, jIndex))
                        {
                            eliminater->m_elierList[jIndex] = eliee;
                            eliminatee->m_isElierList[i] = true;
                            ++jIndex;
                            break;
                        }
                    }
                    continue;
                }
                break;
            }
        }

        pthread_barrier_wait(&startBarrier);
        if (jIndex < 0)
        {
            if (t_id != 0)pthread_exit(NULL);
            return NULL;
        }

        while (eeIndex < Nee)
        {
            //动态划分
            pthread_mutex_lock(&mutexI);
            while (eeIndex < Nee && (eliminatee->m_isElierList[eeIndex] || !getBit(eliminatee->m_elieeList[eeIndex], jIndex)))++eeIndex;
            int i = eeIndex++;
            pthread_mutex_unlock(&mutexI);
            if (i >= Nee) break;

            //进行消元
            eliee = eliminatee->m_elieeList[i];
            elier = eliminater->m_elierList[jIndex];
            for (int k = (eliminatee->m_Ncol >> 5) - 4; k >= 0; k -= 8)
            {
                velier = vld1q_u32(elier + k);
                veliee = vld1q_u32(eliee + k);
                vst1q_u32(eliee + k, veorq_u32(velier, veliee));
            }
        }
        pthread_barrier_wait(&finishBarrier);
    }
}

void pthreadNeonGroebner(Eliminater*& m_eliminater, Eliminatee*& m_eliminatee, int m_threadNum)
{
    // 初始化
    jIndex = Nmat;
    eliminater = m_eliminater;
    eliminatee = m_eliminatee;
    threadNum = m_threadNum;
    pthread_barrier_init(&startBarrier, NULL, threadNum);
    pthread_barrier_init(&finishBarrier, NULL, threadNum);
    pthread_mutex_init(&mutexI, NULL);

    for (int th = 0; th < threadNum; th++)
    {
        thread_param_t[th].t_id = th;
        if (th != 0)pthread_create(&threads[th], NULL, subthreadNeonGroebner, (void*)&thread_param_t[th]);
    }
    subthreadNeonGroebner((void*)&thread_param_t[0]);

    //销毁
    pthread_barrier_destroy(&startBarrier);
    pthread_barrier_destroy(&finishBarrier);
    pthread_mutex_destroy(&mutexI);
}
#endif