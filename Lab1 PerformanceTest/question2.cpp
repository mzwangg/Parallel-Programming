#include<iostream>
#include<algorithm>
#include<string.h>
#include<chrono>
#include"mpi.h"
using namespace std;

const int repeatNum = 5;
const int N = 1e7;
const int unrollN = 8;

//平凡算法
template<typename T>void sumSerial(T* arr)
{
    double time = 0x3f3f3f3f3f3f3f3f;
    for (int step = 0; step < repeatNum; ++step)
    {
        chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

        //主程序
        T sum = 0;
        for (int i = 0; i < N; ++i) sum += arr[i];

        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
        chrono::duration<double, std::milli> time_span = t2 - t1;
        time = min(time, time_span.count() / 1000);
        cout << "sum:" << sum << endl;
    }
    cout << "time:" << time << endl;
}


//循环展开算法
template<typename T>void sumSerialUnroll(T* arr)
{
    double time = 0x3f3f3f3f3f3f3f3f;
    for (int step = 0; step < repeatNum; ++step)
    {
        T sumArr[9], sum = 0;
        memset(sumArr, 0, 9 * sizeof(float));
        chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

        //主程序
        switch (unrollN)
        {
        case 2:
        {
            for (int i = 0; i < N % unrollN; ++i) sum += arr[i];
            for (int i = N % unrollN; i < N; i += unrollN)
            {
                sumArr[0] += arr[i];
                sumArr[1] += arr[i + 1];
            }
            sum = sumArr[0] + sumArr[1];
            break;
        }
        case 3:
        {
            for (int i = 0; i < N % unrollN; ++i) sum += arr[i];
            for (int i = N % unrollN; i < N; i += unrollN)
            {
                sumArr[0] += arr[i];
                sumArr[1] += arr[i + 1];
                sumArr[2] += arr[i + 2];
            }
            sum = sumArr[0] + sumArr[1] + sumArr[2];
            break;
        }
        case 4:
        {
            for (int i = 0; i < N % unrollN; ++i) sum += arr[i];
            for (int i = N % unrollN; i < N; i += unrollN)
            {
                sumArr[0] += arr[i];
                sumArr[1] += arr[i + 1];
                sumArr[2] += arr[i + 2];
                sumArr[3] += arr[i + 3];
            }
            sum = sumArr[0] + sumArr[1] + sumArr[2] + sumArr[3];
            break;
        }
        case 5:
        {
            for (int i = 0; i < N % unrollN; ++i) sum += arr[i];
            for (int i = N % unrollN; i < N; i += unrollN)
            {
                sumArr[0] += arr[i];
                sumArr[1] += arr[i + 1];
                sumArr[2] += arr[i + 2];
                sumArr[3] += arr[i + 3];
                sumArr[4] += arr[i + 4];
            }
            sum = sumArr[0] + sumArr[1] + sumArr[2] + sumArr[3] + sumArr[4];
            break;
        }
        case 6:
        {
            for (int i = 0; i < N % unrollN; ++i) sum += arr[i];
            for (int i = N % unrollN; i < N; i += unrollN)
            {
                sumArr[0] += arr[i];
                sumArr[1] += arr[i + 1];
                sumArr[2] += arr[i + 2];
                sumArr[3] += arr[i + 3];
                sumArr[4] += arr[i + 4];
                sumArr[5] += arr[i + 5];
            }
            sum = sumArr[0] + sumArr[1] + sumArr[2] + sumArr[3] + sumArr[4] + sumArr[5];
            break;
        }
        case 7:
        {
            for (int i = 0; i < N % unrollN; ++i) sum += arr[i];
            for (int i = N % unrollN; i < N; i += unrollN)
            {
                sumArr[0] += arr[i];
                sumArr[1] += arr[i + 1];
                sumArr[2] += arr[i + 2];
                sumArr[3] += arr[i + 3];
                sumArr[4] += arr[i + 4];
                sumArr[5] += arr[i + 5];
                sumArr[6] += arr[i + 6];
            }
            sum = sumArr[0] + sumArr[1] + sumArr[2] + sumArr[3] + sumArr[4] + sumArr[5] + sumArr[6];
            break;
        }
        case 8:
        {
            for (int i = 0; i < N % unrollN; ++i) sum += arr[i];
            for (int i = N % unrollN; i < N; i += unrollN)
            {
                sumArr[0] += arr[i];
                sumArr[1] += arr[i + 1];
                sumArr[2] += arr[i + 2];
                sumArr[3] += arr[i + 3];
                sumArr[4] += arr[i + 4];
                sumArr[5] += arr[i + 5];
                sumArr[6] += arr[i + 6];
                sumArr[7] += arr[i + 7];
            }
            sum = sumArr[0] + sumArr[1] + sumArr[2] + sumArr[3] + sumArr[4] + sumArr[5] + sumArr[6] + sumArr[7];
            break;
        }
        case 9:
        {
            for (int i = 0; i < N % unrollN; ++i) sum += arr[i];
            for (int i = N % unrollN; i < N; i += unrollN)
            {
                sumArr[0] += arr[i];
                sumArr[1] += arr[i + 1];
                sumArr[2] += arr[i + 2];
                sumArr[3] += arr[i + 3];
                sumArr[4] += arr[i + 4];
                sumArr[5] += arr[i + 5];
                sumArr[6] += arr[i + 6];
                sumArr[7] += arr[i + 7];
                sumArr[8] += arr[i + 8];
            }
            sum = sumArr[0] + sumArr[1] + sumArr[2] + sumArr[3] + sumArr[4] + sumArr[5] + sumArr[6] + sumArr[7] + sumArr[8];
            break;
        }
        }
        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
        chrono::duration<double, std::milli> time_span = t2 - t1;
        time = min(time, time_span.count() / 1000);
        cout << "sum:" << sum << endl;
    }
    cout << "time:" << time << endl;
}


//并行算法
template<typename T>void sumParallel(T* arr)
{
    int myid, numprocs;
    int lenth, end;
    float sum, totalSum;
    double local_start, local_finish, local_time, time;
    MPI_Init(NULL, NULL);

    //计时开始
    MPI_Barrier(MPI_COMM_WORLD);
    local_start = MPI_Wtime();

    //主程序
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    lenth = N / numprocs;
    end = (myid == numprocs - 1 ? N : (myid + 1) * lenth);
    sum = 0;
    totalSum = 0;
    for (int i = end - lenth; i < end; ++i) sum += arr[i];
    MPI_Reduce(&sum, &totalSum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

    //计时结束，计算时间
    local_finish = MPI_Wtime();
    local_time = local_finish - local_start;
    MPI_Reduce(&local_time, &time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (myid == 0)
    {
        cout << "sum:" << totalSum << endl;
        cout << "time:" << time << endl;
    }

    MPI_Finalize();
}

int main()
{
    float* arr = new float[N];
    for (int i = 0; i < N; ++i)arr[i] = 1;
    
    //选择进行的实验

    //sumSerial<float>(arr);
    //sumSerialUnroll<float>(arr);
    sumParallel<float>(arr);
}
