#include <iostream>
#include <CL/sycl.hpp>
#include <cstring>
#include "sys/time.h"

using namespace std;
using namespace sycl;

#define BLOCK 64

//Ϊ������ּ�����Ϊnaf������ȼ�����������������������Ǿ����ٽ����б任�ķ������ɲ�������
void dataGenerate(float* mat, float* ans, int N)
{
    memset(mat, 0, N * N * sizeof(float));//��Ȼmat�Ƕ���ָ�룬����һ��ָ����ʵ�����ڵ�
    for (int i = 0; i < N; ++i)
    {
        mat[i * N + i] = 1.0f;
        for (int j = i + 1; j < N; ++j)
            mat[i * N + j] = rand() & 0xF;//ʹ������������0��15֮��,��ֹ��������б任��������ʱ�������
    }

    //��ǰ�Ľ����Ϊ���յĴ�
    memcpy(ans, mat, N * N * sizeof(float));//��Ȼans�Ƕ���ָ�룬����һ��ָ����ʵ�����ڵ�

    int inv_rate = N >> 3;
    bool reGenerate = false;
    while (true)
    {
        //��������������ݣ���mat����Ϊans
    REGEn:
        if (reGenerate) memcpy(mat, ans, N * N * sizeof(float));
        for (int k = 0; k < N; ++k)
            for (int i = k; i < N; ++i)
            {
                if (rand() % inv_rate) continue;
                for (int j = 0; j < N; ++j)
                {
                    mat[i * N + j] += mat[k * N + j];
                    if (mat[i * N + j] > 16777216)//16777216Ϊ�������ܾ�ȷ��ʾ�������
                    {
                        //���mat[i][j]���������ݲ��ھ�ȷ����ʱ�����н������ʣ����¿�ʼ
                        reGenerate = true;
                        inv_rate <<= 1;
                        goto REGEn;
                    }
                }
            }
        break;
    }
}

//��ӡ��������ڵ���
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

//��һ��Ԫ��Ϊ���������GPU���и�˹��Ԫ
double gpu_item_kernel(float* mat, int N, sycl::queue& q) {
    struct timeval start, end;
    gettimeofday(&start, NULL);

    for (int k = 0; k < N; ++k)
    {
        float pivot = mat[k * N + k];
        for (int j = k + 1; j < N; ++j)mat[k * N + j] /= pivot;
        mat[k * N + k] = 1;

        auto e = q.submit([&](sycl::handler& h) {
            h.parallel_for(range<2>(N - k - 1, N - k - 1), [=](id<2> index) {//��������ֵ���Ԫ��
                    int i = index[0] + k + 1;  // ��һ��ά�ȵ�ѭ������
                    int j = index[1] + k + 1;  // �ڶ���ά�ȵ�ѭ������
                    mat[i * N + j] -= mat[k * N + j] * mat[i * N + k];
                });
            });
        e.wait();
        for (int i = k + 1; i < N; ++i)mat[i * N + k] = 0;
    }

    gettimeofday(&end, NULL);
    return (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_usec - start.tv_usec) / 1000.0;
}

//��һ��Ԫ��Ϊ���������GPU���и�˹��Ԫ������ȫ������GPU���л�
double gpu_item_plus_kernel(float* mat, int N, sycl::queue& q) {
    struct timeval start, end;
    gettimeofday(&start, NULL);

    for (int k = 0; k < N; ++k)
    {
        float pivot = mat[k * N + k];
        auto e1 = q.submit([&](sycl::handler& h) {
            h.parallel_for(range<1>(N - k - 1), [=](id<1> i) {//��������ֵ���Ԫ��
                mat[k * N + i + k + 1] /= pivot;
                });
            });
        e1.wait();
        mat[k * N + k] = 1;

        auto e2 = q.submit([&](sycl::handler& h) {
            h.parallel_for(range<2>(N - k - 1, N - k - 1), [=](id<2> index) {//��������ֵ���Ԫ��
                int i = index[0] + k + 1;  // ��һ��ά�ȵ�ѭ������
                int j = index[1] + k + 1;  // �ڶ���ά�ȵ�ѭ������
                mat[i * N + j] -= mat[k * N + j] * mat[i * N + k];
                });
            });
        e2.wait();

        auto e3 = q.submit([&](sycl::handler& h) {
            h.parallel_for(range<1>(N - k - 1), [=](id<1> i) {//��������ֵ���Ԫ��
                mat[(i+k+1)*N+k] = 0;
                });
            });
        e3.wait();
    }

    gettimeofday(&end, NULL);
    return (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_usec - start.tv_usec) / 1000.0;
}

//��һ��Ϊ���������GPU���и�˹��Ԫ
double gpu_line_kernel(float* mat, int N, sycl::queue& q) {
        struct timeval start, end;
        gettimeofday(&start, NULL);
    
        for (int k = 0; k < N; ++k)
        {
            float pivot = mat[k * N + k];
            for (int j = k + 1; j < N; ++j)mat[k * N + j] /= pivot;
            mat[k * N + k] = 1;
    
            auto e = q.submit([&](sycl::handler& h) {
                h.parallel_for(range<1>(N - k - 1), [=](id<1> i) {//���������Ϊ����ΪN��һά����
                    float* row1 = mat + k * N;//������
                    float* row2 = mat + (i + k + 1) * N;//��ǰ��
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

//��һ����Ϊ���������GPU���и�˹��Ԫ
double gpu_block_kernel(float* mat, int N, sycl::queue& q) {
        struct timeval start, end;
        gettimeofday(&start, NULL);
    
        for (int k = 0; k < N; ++k)
        {
            float pivot = mat[k * N + k];
            for (int j = k + 1; j < N; ++j)mat[k * N + j] /= pivot;
            mat[k * N + k] = 1;
    
            auto e = q.submit([&](sycl::handler& h) {
                h.parallel_for(range<2>(N - k - 1,BLOCK), [=](id<2> index) {//���������Ϊ1*N/BLOCK��С�Ŀ�
                    float* row1 = mat + k * N;//������
                    float* row2 = mat + (index[0] + k + 1) * N;//��ǰ��
                    float pivot = row2[k];
                    for (int j = k + 1 +index[1]; j < N; j+=BLOCK) {//ѭ����������
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

//ͨ��CPU�����㷨���и�˹��Ԫ
double cpu_kernel(float* mat, int N) {
    struct timeval start, end;
    gettimeofday(&start, NULL);

    for (int k = 0; k < N; ++k)
    {
        //���н��й�һ��
        float pivot = mat[k * N + k];
        for (int j = k + 1; j < N; ++j)mat[k * N + j] /= pivot;
        mat[k * N + k] = 1.0f;

        //�����½ǵľ��������Ԫ
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

//�Ƚϴ��Ƿ���ȷ
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

    // GPU����
    int warmup = 5;//�Ƚ�����������������ʱ

    //�Ե���Ԫ��Ϊ���������GPU���и�˹��Ԫ
    for (int run = 0; run < iterations + warmup; run++) {
        memcpy(gpu_mat, mat, N * N * sizeof(float));
        float duration = gpu_item_kernel(gpu_mat, N, q);
        if (run >= warmup) duration_gpu_item += duration;
    }
    duration_gpu_item = duration_gpu_item / iterations;

    //�Ե���Ԫ��Ϊ���������GPU���и�˹��Ԫ�������м����ʹ��GPU����
    for (int run = 0; run < iterations + warmup; run++) {
        memcpy(gpu_mat, mat, N * N * sizeof(float));
        float duration = gpu_item_plus_kernel(gpu_mat, N, q);
        if (run >= warmup) duration_gpu_item_plus += duration;
    }
    duration_gpu_item_plus = duration_gpu_item_plus / iterations;

    //��һ��Ϊ���������GPU���и�˹��Ԫ
    for (int run = 0; run < iterations + warmup; run++) {
        memcpy(gpu_mat, mat, N * N * sizeof(float));
        float duration = gpu_line_kernel(gpu_mat, N, q);
        if (run >= warmup) duration_gpu_line += duration;
    }
    duration_gpu_line = duration_gpu_line / iterations;

    //��һ����Ϊ���������GPU���и�˹��Ԫ
    for (int run = 0; run < iterations + warmup; run++) {
        memcpy(gpu_mat, mat, N * N * sizeof(float));
        float duration = gpu_block_kernel(gpu_mat, N, q);
        if (run >= warmup) duration_gpu_block += duration;
    }
    duration_gpu_block = duration_gpu_block / iterations;

    // CPU����
    for (int run = 0; run < iterations + warmup; run++) {
        memcpy(cpu_mat, mat, N * N * sizeof(float));
        float duration = cpu_kernel(cpu_mat, N);
        if (run >= warmup) duration_cpu += duration;
    }
    duration_cpu = duration_cpu / iterations;
    printf("%d\t%f\t%f\t%f\t%f\t%f\n", N, duration_gpu_item, duration_gpu_item_plus,
        duration_gpu_line, duration_gpu_block, duration_cpu);

    // �Ƚ����ߵĴ��Ƿ���ȷ
    //verify(ans, gpu_mat, N);
    //verify(ans, cpu_mat, N);

    free(gpu_mat, q);
    free(cpu_mat, q);
    free(mat, q);
    free(ans, q);
}

int main() {

    //������������enable_profiling���ԣ��Բ���ʱ��
    auto propList = cl::sycl::property_list{ cl::sycl::property::queue::enable_profiling() };
    queue my_gpu_queue(cl::sycl::gpu_selector_v, propList);//����gpu�Ķ���

    //��ӡ�豸��Ϣ
    cout << "Select device: " << my_gpu_queue.get_device().get_info<info::device::name>() << "\n";
    cout << "N\tgpu_item\tgpu_item_plus\tgpu_line\tgpu_block\tcpu\n";

    for (int N = 64; N <= 2048; N += 64)
        gemm(N, 6, my_gpu_queue);

    return 0;
}