#include "config.h"
#include "EliminateClass.hpp"

int Nmat;//矩阵阶数
int Nee;//被消元行行数
int Ner;//消元子行数

void dataInit(string DATA, Eliminater*& eliminater, Eliminatee*& eliminatee)
{
    stringstream ss(DATA);
    char c;
    ss >> Nmat >> c >> Nmat >> c >> Ner >> c >> Nee;
    if (DATA == "8_23045_18748_14325") Nmat = 23075;
    if (DATA == "9_37960_29304_14921") Nee = 14291;
    eliminater = new Eliminater(DATAPATH + DATA + Para1, Nmat, Ner);
    eliminatee = new Eliminatee(DATAPATH + DATA + Para2, Nmat, Nee);
}

//测试函数
void timingAll(string nameArr[], void(*funArr[])(Eliminater*&, Eliminatee*&, int), int funNum)
{
    Eliminater* globalEliminater;//全局消元子类
    Eliminatee* globalEliminatee;//全局被消元行类
    Eliminater* eliminater;//消元子类
    Eliminatee* eliminatee;//被消元行类
    for (int i = 0; i < funNum; ++i)cout << nameArr[i] << '\t';
    cout << endl;

    for (int dataIndex = 0; dataIndex < dataSize; ++dataIndex)
    {
        string DATA = dataArr[dataIndex];
        dataInit(DATA, globalEliminater, globalEliminatee);
        for (int funIndex = 0; funIndex < funNum; ++funIndex)
        {
            double totalTime = 0.0f;
            for (int step = 0; step < REPEAT_NUM; ++step)
            {
                eliminater = new Eliminater(*globalEliminater);
                eliminatee = new Eliminatee(*globalEliminatee);
                struct timeval start, end;
                gettimeofday(&start, NULL);

                //程序的主体部分
                funArr[funIndex](eliminater, eliminatee, 8);

                gettimeofday(&end, NULL);

#ifdef CHECK
                eliminatee->check(DATAPATH + DATA + Para3);
#endif

                totalTime += (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
                delete eliminater;
                delete eliminatee;
            }
            cout << totalTime / (double)REPEAT_NUM << '\t';
#ifdef CHECK
            cout << "check result: Right!\n";
#endif 
        }
        delete globalEliminater;
        delete globalEliminatee;
        cout << '\n';
    }
}

void timing(void(*fun)(Eliminater*&, Eliminatee*&, int),int dataIndex)
{
    Eliminater* eliminater;//消元子类
    Eliminatee* eliminatee;//被消元行类

    string DATA = dataArr[dataIndex - 1];
    dataInit(DATA, eliminater, eliminatee);
    cout << "数据加载完成\n";
    cout << "Nmat:" << Nmat << " Ner:" << Ner << " Nee:" << Nee << endl;

    double totalTime = 0.0f;
    struct timeval start, end;
    gettimeofday(&start, NULL);

    //程序的主体部分
    fun(eliminater, eliminatee, 8);

    gettimeofday(&end, NULL);

#ifdef CHECK
    eliminatee->check(DATAPATH + DATA + "\\消元结果.txt");
#endif

    totalTime += (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    delete eliminater;
    delete eliminatee;
    cout << totalTime << '\n';

#ifdef CHECK
    cout<<"check result: Right!\n";
#endif 
}

// void serialGroebner(Eliminater*& eliminater, Eliminatee*& eliminatee, int threadNum=1)
// {
//     uint* elier, * eliee;
//     for (int i = 0; i < Nee; ++i)
//     {
//         eliee = eliminatee->m_elieeList[i];
//         for (int j = Nmat - 1; j >= 0; --j)
//         {
//             if (getBit(eliee, j))
//             {
//                 elier = eliminater->m_elierList[j];
//                 if (elier != nullptr)
//                 {
//                     int strLen = eliminatee->m_Ncol >> 5;
//                     for (int k = 0; k < strLen; ++k) eliee[k] ^= elier[k];
//                 }
//                 else
//                 {
//                     eliminater->m_elierList[j] = eliee;
//                     eliminatee->m_isElierList[i] = true;
//                     break;
//                 }
//             }
//         }
//     }
// }

void serialGroebner(Eliminater*& eliminater, Eliminatee*& eliminatee, int threadNum = 1)
{
   uint* elier, * eliee;
   for (int j = Nmat - 1; j >= 0; j--)
   { // 遍历消元子
       elier = eliminater->m_elierList[j];
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
       for (int i = 0; i < Nee; ++i)
       { // 遍历被消元行
           if (eliminatee->m_isElierList[i])
               continue;
           eliee = eliminatee->m_elieeList[i];
           if (getBit(eliee, j))
           { // 如果当前行需要消元
               int strLen = eliminatee->m_Ncol >> 5;
               for (int k = 0; k < strLen; ++k)
                   eliee[k] ^= elier[k];
           }
       }
   }
}

#ifndef __ARM_NEON
void avx2Groebner(Eliminater*& eliminater, Eliminatee*& eliminatee, int threadNum = 1)
{
    uint* elier, * eliee;
    __m256i velier, veliee;
    for (int i = 0; i < Nee; ++i)
    {
        eliee = eliminatee->m_elieeList[i];
        for (int j = Nmat - 1; j >= 0; --j)
        {
            if (getBit(eliee, j))
            {
                elier = eliminater->m_elierList[j];
                if (elier != nullptr)
                {
                    for (int k = (eliminatee->m_Ncol>>5)-8; k >= 0 ; k-=8)
                    {
                        velier = _mm256_load_si256((__m256i*)(elier + k));
                        veliee = _mm256_load_si256((__m256i*)(eliee + k));
                        _mm256_store_si256((__m256i*)(eliee + k), _mm256_xor_si256(velier, veliee));
                    }
                }
                else
                {
                    eliminater->m_elierList[j] = eliee;
                    eliminatee->m_isElierList[i] = true;
                    break;
                }
            }
        }
    }
}

#else
void neonGroebner(Eliminater*& eliminater, Eliminatee*& eliminatee, int threadNum = 1)
{
    uint* elier, * eliee;
    uint32x4_t velier, veliee;
    for (int i = 0; i < Nee; ++i)
    {
        eliee = eliminatee->m_elieeList[i];
        for (int j = Nmat - 1; j >= 0; --j)
        {
            if (getBit(eliee, j))
            {
                elier = eliminater->m_elierList[j];
                if (elier != nullptr)
                {
                    int strLen = eliminatee->m_Ncol >> 7;
                   for (int k = 0; k < strLen; ++k)
                   {
                       velier = vld1q_u32(elier + (k << 2));
                       veliee = vld1q_u32(eliee + (k << 2));
                       vst1q_u32(eliee + (k << 2), veorq_u32(velier, veliee));
                   }
                }
                else
                {
                    eliminater->m_elierList[j] = eliee;
                    eliminatee->m_isElierList[i] = true;
                    break;
                }
            }
        }
    }
}
#endif