#include "GaussBaseFunc.h"

__global__ void division_kernel(float* mat, int k, int n) 
{
	float pivot = mat[k * n + k];
	//将i初始化为线程索引，这里线程格大小为（1，1），所以不需要计算gridDim有关信息
	for (int i = k + 1 + threadIdx.x; i  < n; i += blockDim.x)//通过循环划分的方式划分数据
		mat[k * n + i] /= pivot;
}

__global__ void eliminate_kernel(float* mat, int k, int n) 
{
	//零号进程负责将归一化部分未处理的对角线元素设为 1
	if (blockIdx.x == 0 && threadIdx.x == 0) mat[k * n + k] = 1;

	for (int i = k + 1 + blockIdx.x; i < n; i += gridDim.x){
		float pivot = mat[i * n + k];
		for (int j = k + 1 + threadIdx.x; j < n; j += blockDim.x) {
			mat[i * n + j] -= pivot * mat[k * n + j];
		}
		__syncthreads();//块内同步
		if (threadIdx.x == 0) mat[i * n + k] = 0;//将枢值下方对应位置清0
	}	
}

double cuda(int n, float* mat, int blockSize)
{
	size_t size = n * n * sizeof(float);

	//开始计时
	float elapsedTime;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	// 申请GPU内存并将数据从主机内存复制到GPU内存
	float* d_m;
	cudaMalloc((void**)&d_m, n * n * sizeof(float));
	cudaMemcpy(d_m, mat, size, cudaMemcpyHostToDevice);

	for (int k = 0; k < n; k++) {
		division_kernel << <1, blockSize >> > (d_m, k, n);
		cudaDeviceSynchronize();
		eliminate_kernel << <1, blockSize >> > (d_m, k, n);
		cudaDeviceSynchronize();
	}

	// 将数据从GPU内存复制到主机内存
	cudaMemcpy(mat, d_m, size, cudaMemcpyDeviceToHost);

	//结束计时
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);

	//释放内存
	cudaFree(d_m);

	return elapsedTime / 1000.0;
}

double cuda_plus(int n, float* mat, int blockSize) {
	size_t size = n * n * sizeof(float);

	//开始计时
	float elapsedTime;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	//得到多处理器数量，从而确定最佳网格维度
	int deviceId;
	int numberOfSMs;
	cudaGetDevice(&deviceId);
	cudaDeviceGetAttribute(&numberOfSMs, cudaDevAttrMultiProcessorCount, deviceId);

	// 申请GPU内存并将数据从主机内存复制到GPU内存
	float* d_m;
	cudaMalloc((void**)&d_m, n * n * sizeof(float));
	cudaMemcpy(d_m, mat, size, cudaMemcpyHostToDevice);

	//执行核函数
	for (int k = 0; k < n; k++) {
		division_kernel << <1, blockSize >> > (d_m, k, n);
		cudaDeviceSynchronize();
		eliminate_kernel << <numberOfSMs, blockSize >> > (d_m, k, n);
		cudaDeviceSynchronize();
	}

	// 将数据从GPU内存复制到主机内存
	cudaMemcpy(mat, d_m, size, cudaMemcpyDeviceToHost);

	//结束计时
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);

	//释放内存
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
