#include "config.h"
#include "EliminateClass.hpp"

extern int Nmat;//�������
extern int Nee;//����Ԫ������
extern int Ner;//��Ԫ������

void openmpGroebner(Eliminater*& eliminater, Eliminatee*& eliminatee, int threadNum)
{
    uint* elier, * eliee;
#pragma omp parallel num_threads(threadNum) private(elier,eliee)
    for (int j = Nmat - 1; j >= 0; j--)
    { // ������Ԫ��
        elier = eliminater->m_elierList[j];
#pragma omp single
        if (elier == nullptr)
        { // ��������ڶ�Ӧ��Ԫ���򽫱���Ԫ������
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
        { // ��������Ԫ��
            if (eliminatee->m_isElierList[i])
                continue;
            eliee = eliminatee->m_elieeList[i];
            elier = eliminater->m_elierList[j];
            if (getBit(eliee, j))
            { // �����ǰ����Ҫ��Ԫ
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
    { // ������Ԫ��
        elier = eliminater->m_elierList[j];
#pragma omp single
        if (elier==nullptr)
        { // ��������ڶ�Ӧ��Ԫ���򽫱���Ԫ������
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
        { // ��������Ԫ��
            if (eliminatee->m_isElierList[i])
                continue;
            eliee = eliminatee->m_elieeList[i];
            elier = eliminater->m_elierList[j];
            if (getBit(eliee, j))
            { // �����ǰ����Ҫ��Ԫ
                int strLen = eliminatee->m_Ncol >> 5;
#pragma omp simd
                for (int k = 0; k < strLen; ++k)
                    eliee[k] ^= elier[k];
            }
        }
    }
}