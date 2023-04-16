#include "config.h"

#ifdef SERIAL
//串行列主元高斯消去
void serial(float mat[N][N])
{
    for (int k = 0; k < N; ++k)
    {
        //对行进行归一化
        float pivot = mat[k][k];
        for (int j = k + 1; j < N; ++j)mat[k][j] /= pivot;
        mat[k][k] = 1.0;

        //对右下角的矩阵进行消元
        for (int i = k+1; i < N; ++i)
        {
            float pivot = mat[i][k];
            for (int j = k+1; j < N; ++j) mat[i][j] -= mat[k][j] * pivot;
            mat[i][k] = 0.0;
        }
    }
}
#endif
