#include "config.h"

#ifdef SSE_UNALIGNED
//不对齐的SSE高斯消元
void sseUnaligned(float mat[N][N])
{
    __m128 vt,va,vaik,vakj,vaij,vamul;
    for(int k=0;k<N;++k)
    {
        //归一化
        vt = _mm_set1_ps(mat[k][k]);
        int j= k + 1;
        for (; j + 4 <= N; j += 4)
        {
            va = _mm_loadu_ps(mat[k] + j);
            va = _mm_div_ps(va, vt);
            _mm_storeu_ps(mat[k] + j, va);
        }
        for(;j<N;j++) mat[k][j]/=mat[k][k];
        mat[k][k]=1.0;

        //消元
        for(int i=k+1;i<N;++i)
        {
            vaik=_mm_set1_ps(mat[i][k]);
            int j = k + 1;
            for(;j+4<=N;j+=4)
            {
                vakj = _mm_loadu_ps(mat[k] + j);
                vaij = _mm_loadu_ps(mat[i] + j);
                vamul = _mm_mul_ps(vakj, vaik);
                vaij =_mm_sub_ps(vaij,vamul);
                _mm_storeu_ps(mat[i]+j,vaij);
            }
            for(;j<N;j++) mat[i][j]-=mat[k][j]*mat[i][k];
            mat[i][k]=0.0;
        }
    }
}
#endif


#ifdef SSE_ALIGNED
//对齐的SSE高斯消元
void sseAligned(float mat[N][N])
{
    __m128 vt,va,vaik,vakj,vaij,vamul;
    for(int k=0;k<N;++k)
    {
        //归一化
        vt = _mm_set1_ps(mat[k][k]);
        int j= k + 1;
        for (; j < N; j++)
        {
            //16字节对齐
            if (!((ull(mat[k] + j)) & 0xf)) break;
            mat[k][j] /= mat[k][k];
        }
        for (; j + 4 <= N; j += 4)
        {
            va = _mm_load_ps(mat[k] + j);
            va = _mm_div_ps(va, vt);
            _mm_store_ps(mat[k] + j, va);
        }
        for(;j<N;j++) mat[k][j]/=mat[k][k];
        mat[k][k]=1.0;

        //消元
        for(int i=k+1;i<N;++i)
        {
            vaik=_mm_set1_ps(mat[i][k]);
            int j = k + 1;
            for (; j < N; ++j)
            {
                //16字节对齐
                if (!((ull(mat[i] + j)) & 0xf)) break;
                mat[i][j] -= mat[k][j] * mat[i][k];
            }
            for(;j+4<=N;j+=4)
            {
                vakj = _mm_load_ps(mat[k] + j);
                vaij = _mm_load_ps(mat[i] + j);
                vamul = _mm_mul_ps(vakj, vaik);
                vaij =_mm_sub_ps(vaij,vamul);
                _mm_store_ps(mat[i]+j,vaij);
            }
            for(;j<N;j++) mat[i][j]-=mat[k][j]*mat[i][k];
            mat[i][k]=0.0;
        }
    }
}
#endif
