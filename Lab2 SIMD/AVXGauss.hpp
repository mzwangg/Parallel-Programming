#include "config.h"

#ifdef AVX2_UNALIGNED
//不对齐的AVX2高斯消元
void AVX2Unaligned(float mat[N][N])
{
    __m256 vt,va,vaik,vakj,vaij;
    for(int k=0;k<N;++k)
    {
        //归一化
        vt = _mm256_set1_ps(mat[k][k]);
        int j= k + 1;
        for (; j + 8 <= N; j += 8)
        {
            va = _mm256_loadu_ps(mat[k] + j);
            va = _mm256_div_ps(va, vt);
            _mm256_storeu_ps(mat[k] + j, va);
        }
        for(;j<N;j++) mat[k][j]/=mat[k][k];
        mat[k][k]=1.0;

        //消元
        for(int i=k+1;i<N;++i)
        {
            vaik=_mm256_set1_ps(mat[i][k]);
            int j = k + 1;
            for(;j+8<=N;j+=8)
            {
                vakj = _mm256_loadu_ps(mat[k] + j);
                vaij = _mm256_loadu_ps(mat[i] + j);
                vaij =_mm256_fnmadd_ps(vakj,vaik,vaij);
                _mm256_storeu_ps(mat[i]+j,vaij);
            }
            for(;j<N;j++) mat[i][j]-=mat[k][j]*mat[i][k];
            mat[i][k]=0.0;
        }
    }
}
#endif

#ifdef AVX2_ALIGNED
//对齐的AVX2高斯消元
void AVX2Aligned(float mat[N][N])
{
    __m256 vt,va,vaik,vakj,vaij;
    for(int k=0;k<N;++k)
    {
        //归一化
        vt = _mm256_set1_ps(mat[k][k]);
        int j= k + 1;
        for (; j < N; j++)
        {
            //32字节对齐
            if (!((ull(mat[k] + j)) & 0x1f)) break;
            mat[k][j] /= mat[k][k];
        }
        for (; j + 8 <= N; j += 8)
        {
            va = _mm256_load_ps(mat[k] + j);
            va = _mm256_div_ps(va, vt);
            _mm256_store_ps(mat[k] + j, va);
        }
        for(;j<N;j++) mat[k][j]/=mat[k][k];
        mat[k][k]=1.0;

        //消元
        for(int i=k+1;i<N;++i)
        {
            vaik=_mm256_set1_ps(mat[i][k]);
            int j = k + 1;
            for (; j < N; ++j)
            {
                //32字节对齐
                if (!((ull(mat[i] + j)) & 0x1f)) break;
                mat[i][j] -= mat[k][j] * mat[i][k];
            }
            for(;j+8<=N;j+=8)
            {
                vakj = _mm256_load_ps(mat[k] + j);
                vaij = _mm256_load_ps(mat[i] + j);
                vaij =_mm256_fnmadd_ps(vakj,vaik,vaij);
                _mm256_store_ps(mat[i]+j,vaij);
            }
            for(;j<N;j++) mat[i][j]-=mat[k][j]*mat[i][k];
            mat[i][k]=0.0;
        }
    }
}
#endif
