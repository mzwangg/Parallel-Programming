#include "config.h"

#ifdef NEON_UNALIGNED_UNMLS
//不对齐的、不融合乘减的、NEON高斯消元
void neonUnalignedUnMLS(float mat[N][N])
{
    float32x4_t inv_vt,va,vaik,vakj,vaij,vamul;
    for(int k=0;k<N;++k)
    {
        //归一化
        inv_vt = vdupq_n_f32(1.0f / mat[k][k]);
        int j=k+1;
        for(;j+4<=N;j+=4)
        {
            va=vld1q_f32(mat[k]+j);
            va=vmulq_f32(va,inv_vt);
            vst1q_f32(mat[k]+j,va);
        }
        for(;j<N;j++) mat[k][j]/=mat[k][k];
        mat[k][k]=1.0;

        //消元
        for(int i=k+1;i<N;++i)
        {
            vaik=vdupq_n_f32(mat[i][k]);
            int j=k+1;
            for(;j+4<=N;j+=4)
            {
                vakj = vld1q_f32(mat[k] + j);
                vaij = vld1q_f32(mat[i] + j);
                vamul=vmulq_f32(vakj,vaik);
                vaij=vsubq_f32(vaij,vamul);
                vst1q_f32(mat[i]+j,vaij);
            }
            for(;j<N;j++) mat[i][j]-=mat[k][j]*mat[i][k];
            mat[i][k]=0.0;
        }
    }
}
#endif

#ifdef NEON_ALIGNED_UNMLS
//对齐的、不融合乘减的、NEON高斯消元
void neonAlignedUnMLS(float mat[N][N])
{
    float32x4_t inv_vt,va,vaik,vakj,vaij,vamul;
    for(int k=0;k<N;++k)
    {
        //归一化
        inv_vt = vdupq_n_f32(1.0f / mat[k][k]);
        int j=k+1;
        for (; j < N; j++)
        {
            //16字节对齐
            if (!((ull(mat[k] + j)) & 0xf)) break;
            mat[k][j] /= mat[k][k];
        }
        for(;j+4<=N;j+=4)
        {
            va=vld1q_f32(mat[k]+j);
            va=vmulq_f32(va,inv_vt);
            vst1q_f32(mat[k]+j,va);
        }
        for(;j<N;j++) mat[k][j]/=mat[k][k];
        mat[k][k]=1.0;

        //消元
        for(int i=k+1;i<N;++i)
        {
            vaik=vdupq_n_f32(mat[i][k]);
            int j = k + 1;
            for (; j < N; ++j)
            {
                //16字节对齐
                if (!((ull(mat[i] + j)) & 0xf)) break;
                mat[i][j] -= mat[k][j] * mat[i][k];
            }
            for(;j+4<=N;j+=4)
            {
                vakj = vld1q_f32(mat[k] + j);
                vaij = vld1q_f32(mat[i] + j);
                vamul=vmulq_f32(vakj,vaik);
                vaij=vsubq_f32(vaij,vamul);
                vst1q_f32(mat[i]+j,vaij);
            }
            for(;j<N;j++) mat[i][j]-=mat[k][j]*mat[i][k];
            mat[i][k]=0.0;
        }
    }
}
#endif


#ifdef NEON_UNALIGNED_MLS
//不对齐的、融合乘减的、NEON高斯消元
void neonUnalignedMLS(float mat[N][N])
{
    float32x4_t inv_vt,va,vaik,vakj,vaij;
    for(int k=0;k<N;++k)
    {
        //归一化
        inv_vt = vdupq_n_f32(1.0f / mat[k][k]);
        int j=k+1;
        for(;j+4<=N;j+=4)
        {
            va=vld1q_f32(mat[k]+j);
            va=vmulq_f32(va,inv_vt);
            vst1q_f32(mat[k]+j,va);
        }
        for(;j<N;j++) mat[k][j]/=mat[k][k];
        mat[k][k]=1.0;

        //消元
        for(int i=k+1;i<N;++i)
        {
            vaik=vdupq_n_f32(mat[i][k]);
            int j=k+1;
            for(;j+4<=N;j+=4)
            {
                vakj = vld1q_f32(mat[k] + j);
                vaij = vld1q_f32(mat[i] + j);
                vaij = vmlsq_f32(vaij, vakj, vaik);
                vst1q_f32(mat[i] + j, vaij);
            }
            for(;j<N;j++) mat[i][j]-=mat[k][j]*mat[i][k];
            mat[i][k]=0.0;
        }
    }
}
#endif


#ifdef NEON_ALIGNED_MLS
//对齐的、融合乘减的、NEON高斯消元
void neonAlignedMLS(float mat[N][N])
{
float32x4_t inv_vt,va,vaik,vakj,vaij;
    for(int k=0;k<N;++k)
    {
        //归一化
        inv_vt = vdupq_n_f32(1.0f / mat[k][k]);
        int j=k+1;
        for (; j < N; j++)
        {
            //16字节对齐
            if (!((ull(mat[k] + j)) & 0xf)) break;
            mat[k][j] /= mat[k][k];
        }
        for(;j+4<=N;j+=4)
        {
            va=vld1q_f32(mat[k]+j);
            va=vmulq_f32(va,inv_vt);
            vst1q_f32(mat[k]+j,va);
        }
        for(;j<N;j++) mat[k][j]/=mat[k][k];
        mat[k][k]=1.0;

        //消元
        for(int i=k+1;i<N;++i)
        {
            vaik=vdupq_n_f32(mat[i][k]);
            int j = k + 1;
            for (; j < N; ++j)
            {
                //16字节对齐
                if (!((ull(mat[i] + j)) & 0xf)) break;
                mat[i][j] -= mat[k][j] * mat[i][k];
            }
            for(;j+4<=N;j+=4)
            {
                vakj = vld1q_f32(mat[k] + j);
                vaij = vld1q_f32(mat[i] + j);
                vaij = vmlsq_f32(vaij, vakj, vaik);
                vst1q_f32(mat[i] + j, vaij);
            }
            for(;j<N;j++) mat[i][j]-=mat[k][j]*mat[i][k];
            mat[i][k]=0.0;
        }
    }
}
#endif
