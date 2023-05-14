#include "config.h"
#include "EliminateClass.hpp"

extern int Nmat;//矩阵阶数
extern int Nee;//被消元行行数
extern int Ner;//消元子行数

void openmpGroebner(Eliminater*& eliminater, Eliminatee*& eliminatee, int threadNum)
{
    uint* elier, * eliee;
#pragma omp parallel num_threads(threadNum) private(elier,eliee)
    for (int j = Nmat - 1; j >= 0; j--)
    { // 遍历消元子
        elier = eliminater->m_elierList[j];
#pragma omp single
        if (elier == nullptr)
        { // 如果不存在对应消元子则将被消元行升格
            for (int i = 0; i < Nee; ++i)
            {
                if (eliminatee->m_isElierList[i])
                    continue;
                uint* eliee = eliminatee->m_elieeList[i];
                if (getBit(eliee, j))
                {
                    eliminater->m_elierList[j] = eliee;
                    eliminatee->m_isElierList[i] = true;
                    elier = eliee;
                    break;
                }
            }
        }
#pragma omp barrier
#pragma omp for schedule(dynamic, 1)
        for (int i = 0; i < Nee; ++i)
        { // 遍历被消元行
            if (eliminatee->m_isElierList[i])
                continue;
            eliee = eliminatee->m_elieeList[i];
            elier = eliminater->m_elierList[j];
            if (getBit(eliee, j))
            { // 如果当前行需要消元
                int strLen = eliminatee->m_Ncol >> 5;
                for (int k = 0; k < strLen; ++k)
                    eliee[k] ^= elier[k];
            }
        }
    }
}

void openmpSIMDGroebner(Eliminater*& eliminater, Eliminatee*& eliminatee, int threadNum)
{
    uint* elier, * eliee;
#pragma omp parallel num_threads(threadNum) private(elier,eliee)
    for (int j = Nmat - 1; j >= 0; j--)
    { // 遍历消元子
        elier = eliminater->m_elierList[j];
#pragma omp single
        if (elier==nullptr)
        { // 如果不存在对应消元子则将被消元行升格
            for (int i = 0; i < Nee; ++i)
            {
                if (eliminatee->m_isElierList[i])
                    continue;
                uint* eliee = eliminatee->m_elieeList[i];
                if (getBit(eliee, j))
                {
                    eliminater->m_elierList[j] = eliee;
                    eliminatee->m_isElierList[i] = true;
                    elier = eliee;
                    break;
                }
            }
        }
#pragma omp barrier
#pragma omp for schedule(dynamic, 1)
        for (int i = 0; i < Nee; ++i)
        { // 遍历被消元行
            if (eliminatee->m_isElierList[i])
                continue;
            eliee = eliminatee->m_elieeList[i];
            elier = eliminater->m_elierList[j];
            if (getBit(eliee, j))
            { // 如果当前行需要消元
                int strLen = eliminatee->m_Ncol >> 5;
#pragma omp simd
                for (int k = 0; k < strLen; ++k)
                    eliee[k] ^= elier[k];
            }
        }
    }
}