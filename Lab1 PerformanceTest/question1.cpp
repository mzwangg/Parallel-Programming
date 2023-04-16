#include<iostream>
#include<chrono>
#include<algorithm>
#include<fstream>
#include<iomanip>
#include<string.h>
using namespace std;


typedef float TYPE;//可更改矩阵数据类型
const int N = 4096;//可设定矩阵阶数,数组里的N会按顺序执行
const int BLOCKSIZE = 6;//矩阵分块大小
const int repeatNum = 5;//重复计算的次数
const bool SAVE = true;//可指定是否保存矩阵结果
ofstream ansOut("ans.xls", ios::out | ios::trunc);


template<typename T>void column(T* a, T* b, T*& ans)
{
    double time = 0x3f3f3f3f3f3f3f3f;
    if(SAVE) ansOut <<"Column\tN=" << N << endl;
    for (int step = 0; step < repeatNum; ++step)
    {
        chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

        //函数主体
        memset(ans, 0, sizeof(ans[0]) * N);
        for (int j = 0; j < N; ++j)
            for (int i = 0; i < N; ++i)
                ans[j] += a[i * N + j] * b[i];

        //计算用时
        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
        chrono::duration<double, std::milli> time_span = t2 - t1;

        //计算并输出时间
        time = min(time, time_span.count() / CLOCKS_PER_SEC);
    }
    if (SAVE)for (int i = 0; i < N; ++i) ansOut << setprecision(8) << ans[i] << '\t';
    cout <<"time:" << time << endl;
}


template<typename T>void row(T* a, T* b, T*& ans)
{
    double time = 0x3f3f3f3f3f3f3f3f;
    if (SAVE) ansOut << "Row\tN=" << N << endl;
    for (int step = 0; step < repeatNum; ++step)
    {
        chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

        //函数主体
        memset(ans, 0, sizeof(ans[0]) * N);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                ans[j] += a[i * N + j] * b[i];

        //计算用时
        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
        chrono::duration<double, std::milli> time_span = t2 - t1;

        //计算并输出时间
        time = min(time, time_span.count() / CLOCKS_PER_SEC);
    }
    if (SAVE)for (int i = 0; i < N; ++i) ansOut << setprecision(8) << ans[i] << '\t';
    cout << "time:" << time << endl;
}


//以分块矩阵方法计算
template<typename T>void block(T* a, T* b, T*& ans)
{
    double time = 0x3f3f3f3f3f3f3f3f;
    if (SAVE) ansOut << "Block\tN=" << N << endl;
    for (int step = 0; step < repeatNum; ++step)
    {
        chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

        //函数主体
        memset(ans, 0, sizeof(ans[0]) * N);
        for (int k = 0; k < N; k += BLOCKSIZE)
            for (int j = 0; j < N; ++j)
            {
                int iLimit = min(k + BLOCKSIZE, N);
                for (int i = k; i < iLimit; ++i)
                    ans[j] += a[i * N + j] * b[i];
            }

        //计算用时
        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
        chrono::duration<double, std::milli> time_span = t2 - t1;

        //计算并输出时间
        time = min(time, time_span.count() / CLOCKS_PER_SEC);
    }
    if (SAVE)for (int i = 0; i < N; ++i) ansOut << setprecision(8) << ans[i] << '\t';
    cout << "time:" << time << endl;
}


int main()
{
    //声明变量
    TYPE* a = new TYPE[N * N];
    TYPE* b = new TYPE[N];
    TYPE* ans = new TYPE[N];

    //赋值
    for (int i = 0; i < N; ++i)for (int j = 0; j < N; ++j)a[i * N + j] = 1;
    for (int i = 0; i < N; ++i)b[i] = 1;

    //进行计算
    //column<TYPE>(a, b, ans);
    //row<TYPE>(a, b, ans);
    block<TYPE>(a, b, ans);
}
