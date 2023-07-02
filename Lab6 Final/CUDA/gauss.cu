#include "GaussBaseFunc.h"

__global__ void division_kernel(float* mat, int k, int n) 
{
	float pivot = mat[k * n + k];
	//��i��ʼ��Ϊ�߳������������̸߳��СΪ��1��1�������Բ���Ҫ����gridDim�й���Ϣ
	for (int i = k + 1 + threadIdx.x; i  < n; i += blockDim.x)//ͨ��ѭ�����ֵķ�ʽ��������
		mat[k * n + i] /= pivot;
}

__global__ void eliminate_kernel(float* mat, int k, int n) 
{
	//��Ž��̸��𽫹�һ������δ����ĶԽ���Ԫ����Ϊ 1
	if (blockIdx.x == 0 && threadIdx.x == 0) mat[k * n + k] = 1;

	for (int i = k + 1 + blockIdx.x; i < n; i += gridDim.x){
		float pivot = mat[i * n + k];
		for (int j = k + 1 + threadIdx.x; j < n; j += blockDim.x) {
			mat[i * n + j] -= pivot * mat[k * n + j];
		}
		__syncthreads();//����ͬ��
		if (threadIdx.x == 0) mat[i * n + k] = 0;//����ֵ�·���Ӧλ����0
	}	
}

double cuda(int n, float* mat, int blockSize)
{
	size_t size = n * n * sizeof(float);

	//��ʼ��ʱ
	float elapsedTime;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	// ����GPU�ڴ沢�����ݴ������ڴ渴�Ƶ�GPU�ڴ�
	float* d_m;
	cudaMalloc((void**)&d_m, n * n * sizeof(float));
	cudaMemcpy(d_m, mat, size, cudaMemcpyHostToDevice);

	for (int k = 0; k < n; k++) {
		division_kernel << <1, blockSize >> > (d_m, k, n);
		cudaDeviceSynchronize();
		eliminate_kernel << <1, blockSize >> > (d_m, k, n);
		cudaDeviceSynchronize();
	}

	// �����ݴ�GPU�ڴ渴�Ƶ������ڴ�
	cudaMemcpy(mat, d_m, size, cudaMemcpyDeviceToHost);

	//������ʱ
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);

	//�ͷ��ڴ�
	cudaFree(d_m);

	return elapsedTime / 1000.0;
}

double cuda_plus(int n, float* mat, int blockSize) {
	size_t size = n * n * sizeof(float);

	//��ʼ��ʱ
	float elapsedTime;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	//�õ��ദ�����������Ӷ�ȷ���������ά��
	int deviceId;
	int numberOfSMs;
	cudaGetDevice(&deviceId);
	cudaDeviceGetAttribute(&numberOfSMs, cudaDevAttrMultiProcessorCount, deviceId);

	// ����GPU�ڴ沢�����ݴ������ڴ渴�Ƶ�GPU�ڴ�
	float* d_m;
	cudaMalloc((void**)&d_m, n * n * sizeof(float));
	cudaMemcpy(d_m, mat, size, cudaMemcpyHostToDevice);

	//ִ�к˺���
	for (int k = 0; k < n; k++) {
		division_kernel << <1, blockSize >> > (d_m, k, n);
		cudaDeviceSynchronize();
		eliminate_kernel << <numberOfSMs, blockSize >> > (d_m, k, n);
		cudaDeviceSynchronize();
	}

	// �����ݴ�GPU�ڴ渴�Ƶ������ڴ�
	cudaMemcpy(mat, d_m, size, cudaMemcpyDeviceToHost);

	//������ʱ
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);

	//�ͷ��ڴ�
	cudaFree(d_m);

	return elapsedTime / 1000.0;
}

int main()
{
	freopen("out.xls", "w", stdout);
	//timing(cuda, 2048, 1024);
	//timing(serial, 2048, 1024);

	string nameArr[] = { "cuda_plus" };
	double(*funcArr[])(int, float*, int) = {cuda_plus };
	//timingAllMatSize(nameArr, funcArr, 2, 2048, 1024);
	timingAllBlockSize(nameArr, funcArr, 1, 2048, 1024);
}
