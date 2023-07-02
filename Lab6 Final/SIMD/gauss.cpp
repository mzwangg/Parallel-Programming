#include "config.h"
#include "SerialGauss.hpp"
#include "NeonGauss.hpp"
#include "SSEGauss.hpp"
#include "AVXGauss.hpp"

//数组进行32位字节对齐
alignas(32) float matrix[N][N];
alignas(32) float ans[N][N];

#ifdef CHECK
//打印结果，用于调试
void printMatrix(const char* name, float mat[N][N])
{
    printf("%s：\n", name);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            printf("%.1f\t", mat[i][j]);
        printf("\n");
    }
}

//检查结果
void check(float mat[N][N], float ans[N][N])
{
    //N太大时打印并不好观察，故不再打印结果
    if (N <= MAX_PRINT_N)
    {
        printMatrix("ans", ans);
        printMatrix("matrix", mat);
    }

    //只要mat矩阵的秩与ans相同，则代表结果正确
    for(int i=0;i<N;++i)
        for(int j=0;j<N;++j)
            if(abs(mat[i][j] - ans[i][j]) > eps)
            {
                printf("\ncheck result: Wrong!\n");
                printf("matrix[%d][%d]=%f\n",i,j,matrix[i][j]);
                printf("ans[%d][%d]=%f\n",i,j,ans[i][j]);
                return;
            }
    printf("check result: Right!\n");
}
#endif

//为避免出现计算结果为Naf或无穷等极端情况，采用先生成上三角矩阵再进行行变换的方法生成测试样例
void dataGenerate(float mat[N][N], float ans[N][N] = NULL)
{
    memset(mat, 0, N*N*sizeof(float));
    for (int i = 0; i < N; ++i)
    {
        mat[i][i] = 1.0;
        for (int j = i + 1; j < N; ++j)
            mat[i][j] = rand()&0xF;//使得生成数据在0到15之间,防止下面进行行变换生成数据时发生溢出
    }

    //当前的结果即为最终的答案
    memcpy(ans,mat,N*N*sizeof(float));

    int inv_rate = N >> 4;
    bool reGenerate = false;
    while (true)
    {
        //如果重新生成数据，则将mat复制为ans
    REGEN:
        if (reGenerate) memcpy(mat, ans, N * N * sizeof(float));
        for (int k = 0; k < N; ++k)
            for (int i = k + 1; i < N; ++i)
            {
                if (rand() % inv_rate) continue;
                for (int j = 0; j < N; ++j)
                {
                    mat[i][j] += mat[k][j];
                    if (mat[i][j] > 1e24)
                    {
                        //如果mat[i][j]过大，则数据不在精确，此时调低行交换概率，重新开始
                        reGenerate = true;
                        inv_rate<<1;
                        goto REGEN;
                    }
                }
            }
        break;
    }

#ifdef CHECK
    printf("mat[N-1][N-1]:%f\n", mat[N - 1][N - 1]);
    if (N <= MAX_PRINT_N) printMatrix("matrixOrigin", mat);
#endif

}

int main()
{
    dataGenerate(matrix,ans);
    struct timeval start, end;
    gettimeofday(&start, NULL);

    //程序的主体部分
    #ifdef SERIAL
    if(N==64) cout<<"serial\t";
    serial(matrix);

    #elif NEON_UNALIGNED_UNMLS
    if(N==64) cout<<"neonUnalignedUnMLS\t";
    neonUnalignedUnMLS(matrix);

    #elif NEON_ALIGNED_UNMLS
    if(N==64) cout<<"neonAlignedUnMLS\t";
    neonAlignedUnMLS(matrix);

    #elif NEON_UNALIGNED_MLS
    if(N==64) cout<<"neonUnalignedMLS\t";
    neonUnalignedMLS(matrix);

    #elif NEON_ALIGNED_MLS
    if(N==64) cout<<"neonAlignedMLS\t";
    neonAlignedMLS(matrix);

    #elif SSE_UNALIGNED
    if(N==64) cout<<"sseUnaligned\t";
    sseUnaligned(matrix);

    #elif SSE_ALIGNED
    if(N==64) cout<<"sseAligned\t";
    sseAligned(matrix);

    #elif AVX2_UNALIGNED
    if(N==64) cout<<"AVX2Unaligned\t";
    AVX2Unaligned(matrix);

    #elif AVX2_ALIGNED
    if(N==64) cout<<"AVX2Aligned\t";
    AVX2Aligned(matrix);

    #endif

    gettimeofday(&end, NULL);

    #ifdef CHECK
    check(matrix, ans);
    #endif

    double elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    cout <<elapsed_time<<'\t';
    if(N==2048)cout<<'\n';
    return 0;
}
