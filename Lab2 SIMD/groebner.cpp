#include "config.h"
#include "groebnerBase.hpp"

int Nmat;//矩阵阶数
int Ner;//消元子行数
int Nee;//被消元行行数
string path = "/home/data/Groebner/" + string(DATA);

void dataInit(Eliminater*& eliminater, Eliminatee*& eliminatee)
{
    stringstream ss(DATA);
    char c;
    ss >> Nmat >> c >> Nmat >> c >> Ner >> c >> Nee;
    if(strcmp(DATA,"8_23045_18748_14325")==0) Nmat=23075;
    if(strcmp(DATA,"9_37960_29304_14921")==0) Nee=14291;
    eliminater = new Eliminater(path + "/消元子.txt", Nmat, Ner);
    eliminatee = new Eliminatee(path + "/被消元行.txt", Nmat, Nee);
}

void serialGroebner(Eliminater*& eliminater, Eliminatee*& eliminatee)
{
    uint* elier, * eliee;
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
                    int strLen = eliminatee->m_Ncol >> 5;
                    for (int k = 0; k < strLen; ++k) eliee[k] ^= elier[k];
                }
                else
                {
                    eliminater->m_elierList[j] = eliee;
                    break;
                }
            }
        }
    }
}

#ifdef NEON_GROEBNER
void neonGroebner(Eliminater*& eliminater, Eliminatee*& eliminatee)
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
                        velier = vld1q_u32(elier + (k<<2));
                        veliee = vld1q_u32(eliee + (k<<2));
                        vst1q_u32(eliee + (k<<2), veorq_u32(velier, veliee));
                    }
                }
                else
                {
                    eliminater->m_elierList[j] = eliee;
                    break;
                }
            }
        }
    }
}
#endif

#ifdef SSE_GROEBNER
void sseGroebner(Eliminater*& eliminater, Eliminatee*& eliminatee)
{
    uint* elier, * eliee;
    __m128i velier, veliee;
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
                        velier = _mm_load_si128((__m128i*)elier + (k<<2));
                        veliee = _mm_load_si128((__m128i*)eliee + (k<<2));
                        _mm_store_si128((__m128i*)eliee + (k<<2), _mm_xor_si128(velier, veliee));
                    }
                }
                else
                {
                    eliminater->m_elierList[j] = eliee;
                    break;
                }
            }
        }
    }
}
#endif

#ifdef AVX2_GROEBNER
void avx2Groebner(Eliminater*& eliminater, Eliminatee*& eliminatee)
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
                    int strLen = eliminatee->m_Ncol >> 8;
                    for (int k = 0; k < strLen; ++k)
                    {
                        velier = _mm256_load_si256((__m256i*)elier + (k<<3));
                        veliee = _mm256_load_si256((__m256i*)eliee + (k<<3));
                        _mm256_store_si256((__m256i*)eliee + (k<<3), _mm256_xor_si256(velier, veliee));
                    }
                }
                else
                {
                    eliminater->m_elierList[j] = eliee;
                    break;
                }
            }
        }
    }
}
#endif

int main()
{
    Eliminater* eliminater;//消元子类
    Eliminatee* eliminatee;//被消元行类
    dataInit(eliminater, eliminatee);
    struct timeval start, end;
    gettimeofday(&start, NULL);
    cout<<DATA<<'\t';

    #ifdef SERIAL_GROEBNER
    cout<<"serialGroebner\t";
    serialGroebner(eliminater, eliminatee);

    #elif NEON_GROEBNER
    cout<<"neonGroebner\t";
    neonGroebner(eliminater, eliminatee);

    #elif SSE_GROEBNER
    cout<<"sseGroebner\t";
    sseGroebner(eliminater, eliminatee);

    #elif AVX2_GROEBNER
    cout<<"avx2Groebner\t";
    avx2Groebner(eliminater, eliminatee);
    #endif

    gettimeofday(&end, NULL);

    #ifdef CHECK
    eliminatee->check(path + "/消元结果.txt");
    #endif

    double elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    cout << elapsed_time <<'\t';
    cout<<'\n';
    return 0;
}
