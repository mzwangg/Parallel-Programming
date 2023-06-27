#include <chrono>
#include <iostream>
#include <CL/sycl.hpp>

#define random_float() (rand() / double(RAND_MAX))

using namespace std;
using namespace sycl;

#define tileY 2 //工作项行数
#define tileX 2 //工作项列数
#define RATE 4 //私有变量增大比例

//GPU版本
double gpu_kernel(float* A, float* B, float* C,
    int M, int N, int K,
    int BLOCK, sycl::queue& q) {

    // 计算得到全局范围和局部范围
    auto grid_rows = M / tileY;
    auto grid_cols = N / tileX;
    auto local_ndrange = range<2>(BLOCK, BLOCK);
    auto global_ndrange = range<2>(grid_rows, grid_cols);

    double duration = 0.0f;

    auto e = q.submit([&](sycl::handler& h) {
        h.parallel_for<class k_name_t>(
            sycl::nd_range<2>(global_ndrange, local_ndrange), [=](sycl::nd_item<2> index) {

    int row = tileY * index.get_global_id(0);
    int col = tileX * index.get_global_id(1);

    float sum[tileY][tileX] = { 0.0f };
    float subA[tileY * RATE] = { 0.0f };//将sunA扩大RATE倍，存储tileY*RATE大小的矩阵
    float subB[tileX * RATE] = { 0.0f };//将sunB扩大RATE倍，存储RATE*tileX大小的矩阵

    //核心计算
    for (int k = 0; k < N; k += RATE) {
        //将对应的tileY*RATE大小矩阵映射到subA中
        for (int i = 0; i < tileY; ++i)
            for (int j = 0; j < RATE; ++j)
                subA[i * RATE + j] = A[(row + i) * N + k + j];

        //将对应的RATE*tileX大小矩阵映射到subB中
        for (int i = 0; i < RATE; ++i)
            for (int j = 0; j < tileX; ++j)
                subB[i * tileX + j] = B[(k + i) * N + col + j];

        //计算sum=subA*subB
        for (int i = 0; i < tileY; i++)
            for (int j = 0; j < tileX; j++)
                for (int step = 0; step < RATE; ++step)
                    sum[i][j] += subA[i * RATE + step] * subB[step * tileX + j];
    }

    //将sum映射回原矩阵中
    for (int m = 0; m < tileY; m++) {
        for (int p = 0; p < tileX; p++) {
            C[(row + m) * N + col + p] = sum[m][p];
        }
    }

            });
        });
    e.wait();//等待GPU运算完毕

    //测量时间
    duration += (e.get_profiling_info<info::event_profiling::command_end>() -
        e.get_profiling_info<info::event_profiling::command_start>()) / 1000.0f / 1000.0f;

    return(duration);
}

//CPU版本
double cpu_kernel(float* cA, float* cB, float* cC, int M, int N, int K) {

    double duration = 0.0;
    std::chrono::high_resolution_clock::time_point s, e;

    s = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            float sum = 0.0f;
            for (int k = 0; k < K; k++) {
                sum += cA[i * K + k] * cB[k * N + j];
            }
            cC[i * N + j] = sum;
        }
    }
    e = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration<float, std::milli>(e - s).count();

    return(duration);
}

//判断答案是否正确
int verify(float* cpu_res, float* gpu_res, int length) {
    int err = 0;
    for (int i = 0; i < length; i++) {
        if (fabs(cpu_res[i] - gpu_res[i]) > 1e-3) {
            err++;
            printf("\n%lf, %lf", cpu_res[i], gpu_res[i]);
        }
    }
    return(err);
}

int gemm(const int M,
    const int N,
    const int K,
    const int block_size,
    const int iterations,
    sycl::queue& q) {

    cout << "Problem size: c(" << M << "," << N << ") ="
        << " a(" << M << "," << K << ") *"
        << " b(" << K << "," << N << ")\n";
    cout << "Select device: " 
        << q.get_device().get_info<info::device::name>() << "\n";

    auto A = malloc_shared<float>(M * K, q);
    auto B = malloc_shared<float>(K * N, q);
    auto C = malloc_shared<float>(M * N, q);
    auto C_host = malloc_host<float>(M * N, q);

    // 初始化
    for (int i = 0; i < M * K; i++) {
        A[i] = random_float();
    }

    for (int i = 0; i < K * N; i++) {
        B[i] = random_float();
    }

    for (int i = 0; i < M * N; i++) {
        C[i] = 0.0f;
        C_host[i] = 0.0f;
    }

    double flopsPerMatrixMul
        = 2.0 * static_cast<double>(M) * static_cast<double>(N) * static_cast<double>(K);

    double duration_gpu = 0.0f;
    double duration_cpu = 0.0f;

    //gpu版本时间测量
    int warmup = 10;//先进行十次热身，不计入耗时，使得测量结果更加精确
    for (int run = 0; run < iterations + warmup; run++) {
        float duration = gpu_kernel(A, B, C, M, N, K, block_size, q);
        if (run >= warmup) duration_gpu += duration;
    }
    duration_gpu = duration_gpu / iterations;

    //cpu版本时间测量
    warmup = 2;//先进行两次热身，不计入耗时，使得测量结果更加精确
    for (int run = 0; run < iterations / 2 + warmup; run++) {
        float duration = cpu_kernel(A, B, C_host, M, N, K);
        if (run >= warmup) duration_cpu += duration;
    }
    duration_cpu = duration_cpu / (iterations / 2);

    //检查答案是否正确
    int errCode = 0;
    errCode = verify(C_host, C, M * N);
    if (errCode > 0) printf("\nThere are %d errors\n", errCode);

    printf("\nGEMM size M = %d, N = %d, K = %d", M, N, K);
    printf("\nWork-Group size = %d * %d, tile_X = %d, tile_Y = %d", block_size, block_size, tileX, tileY);
    printf("\nPerformance Flops = %lf, \n"
        "GPU Computation Time = %lf (ms); \n"
        "CPU Computaiton Time = %lf (ms); \n",
        flopsPerMatrixMul, duration_gpu, duration_cpu);

    free(A, q);
    free(B, q);
    free(C, q);
    free(C_host, q);

    return(errCode);
}

int main() {

    //将属性列增加enable_profiling属性，以测量时间
    auto propList = cl::sycl::property_list{ cl::sycl::property::queue::enable_profiling() };
    queue my_gpu_queue(cl::sycl::gpu_selector{}, propList);//生成gpu的队列

    int errCode = gemm(512, 512, 512, /* GEMM size, M, N, K */
        4,             /* workgroup size */
        10,            /* repeat time */
        my_gpu_queue);

    return(errCode);
}