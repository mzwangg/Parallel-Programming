#include <iostream>
#include <CL/sycl.hpp>
#include <cstring>
#include "sys/time.h"

using namespace std;
using namespace sycl;

#define BLOCK 64

//为避免出现计算结果为naf或无穷等极端情况，采用先生成上三角矩阵再进行行变换的方法生成测试样例
void dataGenerate(float* mat, float* ans, int N)
{
    memset(mat, 0, N * N * sizeof(float));//虽然mat是二层指针，但是一层指针其实是相邻的
    for (int i = 0; i < N; ++i)
    {
        mat[i * N + i] = 1.0f;
        for (int j = i + 1; j < N; ++j)
            mat[i * N + j] = rand() & 0xF;//使得生成数据在0到15之间,防止下面进行行变换生成数据时发生溢出
    }

    //当前的结果即为最终的答案
    memcpy(ans, mat, N * N * sizeof(float));//虽然ans是二层指针，但是一层指针其实是相邻的

    int inv_rate = N >> 3;
    bool reGenerate = false;
    while (true)
    {
        //如果重新生成数据，则将mat复制为ans
    REGEn:
        if (reGenerate) memcpy(mat, ans, N * N * sizeof(float));
        for (int k = 0; k < N; ++k)
            for (int i = k; i < N; ++i)
            {
                if (rand() % inv_rate) continue;
                for (int j = 0; j < N; ++j)
                {
                    mat[i * N + j] += mat[k * N + j];
                    if (mat[i * N + j] > 16777216)//16777216为浮点数能精确表示的最大数
                    {
                        //如果mat[i][j]过大，则数据不在精确，此时调低行交换概率，重新开始
                        reGenerate = true;
                        inv_rate <<= 1;
                        goto REGEn;
                    }
                }
            }
        break;
    }
}

//打印结果，用于调试
void printMatrix(const char* name, float* mat, int m, int N)
{
    printf("%s:\n", name);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < N; j++)
            printf("%.1f\t", mat[i * N + j]);
        printf("\n");
    }
}

//以一个元素为工作项进行GPU并行高斯消元
double gpu_item_kernel(float* mat, int N, sycl::queue& q) {
    struct timeval start, end;
    gettimeofday(&start, NULL);

    for (int k = 0; k < N; ++k)
    {
        float pivot = mat[k * N + k];
        for (int j = k + 1; j < N; ++j)mat[k * N + j] /= pivot;
        mat[k * N + k] = 1;

        auto e = q.submit([&](sycl::handler& h) {
            h.parallel_for(range<2>(N - k - 1, N - k - 1), [=](id<2> index) {//将工作项划分单个元素
                    int i = index[0] + k + 1;  // 第一个维度的循环变量
                    int j = index[1] + k + 1;  // 第二个维度的循环变量
                    mat[i * N + j] -= mat[k * N + j] * mat[i * N + k];
                });
            });
        e.wait();
        for (int i = k + 1; i < N; ++i)mat[i * N + k] = 0;
    }

    gettimeofday(&end, NULL);
    return (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_usec - start.tv_usec) / 1000.0;
}

//以一个元素为工作项进行GPU并行高斯消元，并将全部进行GPU并行化
double gpu_item_plus_kernel(float* mat, int N, sycl::queue& q) {
    struct timeval start, end;
    gettimeofday(&start, NULL);

    for (int k = 0; k < N; ++k)
    {
        float pivot = mat[k * N + k];
        auto e1 = q.submit([&](sycl::handler& h) {
            h.parallel_for(range<1>(N - k - 1), [=](id<1> i) {//将工作项划分单个元素
                mat[k * N + i + k + 1] /= pivot;
                });
            });
        e1.wait();
        mat[k * N + k] = 1;

        auto e2 = q.submit([&](sycl::handler& h) {
            h.parallel_for(range<2>(N - k - 1, N - k - 1), [=](id<2> index) {//将工作项划分单个元素
                int i = index[0] + k + 1;  // 第一个维度的循环变量
                int j = index[1] + k + 1;  // 第二个维度的循环变量
                mat[i * N + j] -= mat[k * N + j] * mat[i * N + k];
                });
            });
        e2.wait();

        auto e3 = q.submit([&](sycl::handler& h) {
            h.parallel_for(range<1>(N - k - 1), [=](id<1> i) {//将工作项划分单个元素
                mat[(i+k+1)*N+k] = 0;
                });
            });
        e3.wait();
    }

    gettimeofday(&end, NULL);
    return (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_usec - start.tv_usec) / 1000.0;
}

//以一行为工作项进行GPU并行高斯消元
double gpu_line_kernel(float* mat, int N, sycl::queue& q) {
        struct timeval start, end;
        gettimeofday(&start, NULL);
    
        for (int k = 0; k < N; ++k)
        {
            float pivot = mat[k * N + k];
            for (int j = k + 1; j < N; ++j)mat[k * N + j] /= pivot;
            mat[k * N + k] = 1;
    
            auto e = q.submit([&](sycl::handler& h) {
                h.parallel_for(range<1>(N - k - 1), [=](id<1> i) {//将工作项划分为长度为N的一维向量
                    float* row1 = mat + k * N;//枢轴行
                    float* row2 = mat + (i + k + 1) * N;//当前行
                    float pivot = row2[k];
                    for (int j = k + 1; j < N; ++j) {
                        row2[j] -= row1[j] * pivot;
                    }
                    row2[k] = 0.0f;
                });
            });
            e.wait();
        }
    
        gettimeofday(&end, NULL);
        return (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_usec - start.tv_usec) / 1000.0;
}

//以一个块为工作项进行GPU并行高斯消元
double gpu_block_kernel(float* mat, int N, sycl::queue& q) {
        struct timeval start, end;
        gettimeofday(&start, NULL);
    
        for (int k = 0; k < N; ++k)
        {
            float pivot = mat[k * N + k];
            for (int j = k + 1; j < N; ++j)mat[k * N + j] /= pivot;
            mat[k * N + k] = 1;
    
            auto e = q.submit([&](sycl::handler& h) {
                h.parallel_for(range<2>(N - k - 1,BLOCK), [=](id<2> index) {//将工作项划分为1*N/BLOCK大小的块
                    float* row1 = mat + k * N;//枢轴行
                    float* row2 = mat + (index[0] + k + 1) * N;//当前行
                    float pivot = row2[k];
                    for (int j = k + 1 +index[1]; j < N; j+=BLOCK) {//循环划分数据
                        row2[j] -= row1[j] * pivot;
                    }
                });
            });
            e.wait();
            for (int i = k + 1; i < N; ++i)mat[i * N + k] = 0;
        }
    
        gettimeofday(&end, NULL);
        return (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_usec - start.tv_usec) / 1000.0;
}

//通过CPU串行算法进行高斯消元
double cpu_kernel(float* mat, int N) {
    struct timeval start, end;
    gettimeofday(&start, NULL);

    for (int k = 0; k < N; ++k)
    {
        //对行进行归一化
        float pivot = mat[k * N + k];
        for (int j = k + 1; j < N; ++j)mat[k * N + j] /= pivot;
        mat[k * N + k] = 1.0f;

        //对右下角的矩阵进行消元
        for (int i = k + 1; i < N; ++i)
        {
            float pivot = mat[i * N + k];
            for (int j = k + 1; j < N; ++j)
                mat[i * N + j] -= mat[k * N + j] * pivot;
            mat[i * N + k] = 0.0f;
        }
    }

    gettimeofday(&end, NULL);
    return (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_usec - start.tv_usec) / 1000.0;
}

//比较答案是否正确
void verify(float* mat, float* ans, int n) {
    int length = n * n;
    for (int i = 0; i < length; i++) {
        if (fabs(ans[i] - mat[i]) > 1e-3) {
            printf("\nsize:%d\tcheck result: Wrong!\n", n);
            printf("mat[%d][%d]=%f\n", i / n, i % n, mat[i]);
            printf("ans[%d][%d]=%f\n", i / n, i % n, ans[i]);
            exit(0);
        }
    }
}

void gemm(const int N, int iterations, sycl::queue& q) {

    auto gpu_mat = malloc_shared<float>(N * N, q);
    auto cpu_mat = malloc_host<float>(N * N, q);
    auto mat = malloc_host<float>(N * N, q);
    auto ans = malloc_host<float>(N * N, q);
    dataGenerate(mat, ans, N);

    double duration_gpu_item = 0.0f;
    double duration_gpu_item_plus = 0.0f;
    double duration_gpu_line = 0.0f;
    double duration_gpu_block = 0.0f;
    double duration_cpu = 0.0f;

    // GPU计算
    int warmup = 5;//先进行五次热身，不计入耗时

    //以单个元素为工作项进行GPU并行高斯消元
    for (int run = 0; run < iterations + warmup; run++) {
        memcpy(gpu_mat, mat, N * N * sizeof(float));
        float duration = gpu_item_kernel(gpu_mat, N, q);
        if (run >= warmup) duration_gpu_item += duration;
    }
    duration_gpu_item = duration_gpu_item / iterations;

    //以单个元素为工作项进行GPU并行高斯消元，且所有计算均使用GPU进行
    for (int run = 0; run < iterations + warmup; run++) {
        memcpy(gpu_mat, mat, N * N * sizeof(float));
        float duration = gpu_item_plus_kernel(gpu_mat, N, q);
        if (run >= warmup) duration_gpu_item_plus += duration;
    }
    duration_gpu_item_plus = duration_gpu_item_plus / iterations;

    //以一行为工作项进行GPU并行高斯消元
    for (int run = 0; run < iterations + warmup; run++) {
        memcpy(gpu_mat, mat, N * N * sizeof(float));
        float duration = gpu_line_kernel(gpu_mat, N, q);
        if (run >= warmup) duration_gpu_line += duration;
    }
    duration_gpu_line = duration_gpu_line / iterations;

    //以一个块为工作项进行GPU并行高斯消元
    for (int run = 0; run < iterations + warmup; run++) {
        memcpy(gpu_mat, mat, N * N * sizeof(float));
        float duration = gpu_block_kernel(gpu_mat, N, q);
        if (run >= warmup) duration_gpu_block += duration;
    }
    duration_gpu_block = duration_gpu_block / iterations;

    // CPU计算
    for (int run = 0; run < iterations + warmup; run++) {
        memcpy(cpu_mat, mat, N * N * sizeof(float));
        float duration = cpu_kernel(cpu_mat, N);
        if (run >= warmup) duration_cpu += duration;
    }
    duration_cpu = duration_cpu / iterations;
    printf("%d\t%f\t%f\t%f\t%f\t%f\n", N, duration_gpu_item, duration_gpu_item_plus,
        duration_gpu_line, duration_gpu_block, duration_cpu);

    // 比较两者的答案是否正确
    //verify(ans, gpu_mat, N);
    //verify(ans, cpu_mat, N);

    free(gpu_mat, q);
    free(cpu_mat, q);
    free(mat, q);
    free(ans, q);
}

int main() {

    //将属性列增加enable_profiling属性，以测量时间
    auto propList = cl::sycl::property_list{ cl::sycl::property::queue::enable_profiling() };
    queue my_gpu_queue(cl::sycl::gpu_selector_v, propList);//生成gpu的队列

    //打印设备信息
    cout << "Select device: " << my_gpu_queue.get_device().get_info<info::device::name>() << "\n";
    cout << "N\tgpu_item\tgpu_item_plus\tgpu_line\tgpu_block\tcpu\n";

    for (int N = 64; N <= 2048; N += 64)
        gemm(N, 6, my_gpu_queue);

    return 0;
}