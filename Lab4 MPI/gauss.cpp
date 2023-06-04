#define _CRT_SECURE_NO_WARNINGS
#include "config.h"

float* mat[MAXN], * ans[MAXN], * pivotRow, * tempMat, * tempAns;
int n, comm_sz, my_rank;

//打印结果，用于调试
void printMatrix(const char* name, float* mat, int m, int n)
{
	printf("%s:\n", name);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
			printf("%.1f\t", mat[i * n + j]);
		printf("\n");
	}
}

//检查结果
void check(int step)
{

	//判断是否需要检查
#ifndef CHECK
	return;
#endif

	//n太大时打印并不好观察，故不再打印结果
	if (n <= MAX_PRINT_N)
	{
		printMatrix("ans", ans[0], n, n);
		printMatrix("mat", mat[0], n, n);
	}

	//只要mat矩阵的秩与ans相同，则代表结果正确
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			if (abs(mat[i][j] - ans[i][j]) > eps)
			{
				printf("\nsize:%d\tcheck result: Wrong!\n", n);
				printf("mat[%d][%d]=%f\n", i, j, mat[i][j]);
				printf("ans[%d][%d]=%f\n", i, j, ans[i][j]);
				MPI_Finalize();
				exit(0);
			}

	if (step == REPEAT_NUM)
		cout << "\nsize:" << n << "\tcheck result : Right!\n";
}

//为避免出现计算结果为naf或无穷等极端情况，采用先生成上三角矩阵再进行行变换的方法生成测试样例
void dataGenerate()
{
	memset(mat[0], 0, n * n * sizeof(float));//虽然mat是二层指针，但是一层指针其实是相邻的
	for (int i = 0; i < n; ++i)
	{
		mat[i][i] = 1.0f;
		for (int j = i + 1; j < n; ++j)
			mat[i][j] = rand() & 0xF;//使得生成数据在0到15之间,防止下面进行行变换生成数据时发生溢出
	}

	//当前的结果即为最终的答案
	memcpy(ans[0], mat[0], n * n * sizeof(float));//虽然ans是二层指针，但是一层指针其实是相邻的

	int inv_rate = n >> 3;
	bool reGenerate = false;
	while (true)
	{
		//如果重新生成数据，则将mat复制为ans
	REGEn:
		if (reGenerate) memcpy(mat[0], ans[0], n * n * sizeof(float));
		for (int k = 0; k < n; ++k)
			for (int i = k; i < n; ++i)
			{
				if (rand() % inv_rate) continue;
				for (int j = 0; j < n; ++j)
				{
					mat[i][j] += mat[k][j];
					if (mat[i][j] > floatLimit)
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

//测试函数
void timing(double(*fun)(int), int threadNum= THREADNUM)
{
	for (n = 64; n <= MAXN; n += 64)
	{
		double totalTime = 0.0f;
		//直接将整个二维数组声明为一整块连续的内存
		if (my_rank == 0) {
			tempMat = (float*)malloc(n * n * sizeof(float));
			tempAns = (float*)malloc(n * n * sizeof(float));
			for (int i = 0; i < n; ++i)
			{
				mat[i] = tempMat + i * n;
				ans[i] = tempAns + i * n;
			}
		}
		else {
			tempMat = (float*)malloc(n * (n / comm_sz + 1) * sizeof(float));
			for (int i = n / comm_sz; i >= 0; i--)
				mat[i] = tempMat + i * n;
		}
		pivotRow = (float*)malloc(n * sizeof(float));//该一维向量用于消元时传递枢轴所在行

		for (int step = 1; step <= REPEAT_NUM; ++step)
		{
			if (my_rank == 0)dataGenerate();
			MPI_Barrier(MPI_COMM_WORLD);//其他进程等待数据初始化

			//程序的主体部分
			double localTime = fun(threadNum);

			if (my_rank == 0)
			{
				check(step);
				totalTime += localTime;
			}
		}
		if (my_rank == 0) std::cout << totalTime / (double)REPEAT_NUM << "\t";

		//释放堆中的数据
		std::free(pivotRow);
		std::free(tempMat);
		if (my_rank == 0) std::free(tempAns);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	if (my_rank == 0) std::cout << "\n";
}

//测试函数
void timingAll(string name[], double(*funArr[])(int), int funNum, int threadNum = THREADNUM)
{
	//输出表头信息
	if (my_rank == 0) {
		std::cout << "Type/N\t";
		for (n = 64; n <= MAXN; n += 64)
			std::cout << n << '\t';
		std::cout << "\n";
	}

	for (int funIndex = 0; funIndex < funNum; ++funIndex)
	{
		if (my_rank == 0)std::cout << name[funIndex] << '\t';
		for (n = 64; n <= MAXN; n += 64)
		{
			double totalTime = 0.0f;
			//直接将整个二维数组声明为一整块连续的内存
			if (my_rank == 0) {
				tempMat = (float*)malloc(n * n * sizeof(float));
				tempAns = (float*)malloc(n * n * sizeof(float));
				for (int i = 0; i < n; ++i)
				{
					mat[i] = tempMat + i * n;
					ans[i] = tempAns + i * n;
				}
			}
			else {
				tempMat = (float*)malloc(n * (n / comm_sz + 1) * sizeof(float));
				for (int i = n / comm_sz; i >= 0; i--)
					mat[i] = tempMat + i * n;
			}
			pivotRow = (float*)malloc(n * sizeof(float));//该一维向量用于消元时传递枢轴所在行

			for (int step = 0; step < REPEAT_NUM; ++step)
			{
				if (my_rank == 0)dataGenerate();
				MPI_Barrier(MPI_COMM_WORLD);//其他进程等待数据初始化

				//程序的主体部分
				double localTime = funArr[funIndex](threadNum);

				if (my_rank == 0)
				{
					check(step);
					totalTime += localTime;
				}
			}

			//释放堆中的数据
			std::free(pivotRow);
			std::free(tempMat);
			if (my_rank == 0) {
				std::free(tempAns);
				std::cout << totalTime / (double)REPEAT_NUM << "\t";
			}
		}
		if (my_rank == 0) std::cout << "\n";
	}
}

//测试函数
void timingAllThreadNum(string name[], double(*funArr[])(int), int funNum, const int threadArr[], int threadArrSize)
{
	//输出表头信息
	if (my_rank == 0) {
		std::cout << "Type/Thread\t";
		for (int threadIndex = 0; threadIndex < threadArrSize; ++threadIndex)
			std::cout << threadNumArr[threadIndex] << '\t';
		std::cout << "\n";
	}

	//在堆区声明变量
	n = MAXN;
	if (my_rank == 0) {
		tempMat = (float*)malloc(n * n * sizeof(float));
		tempAns = (float*)malloc(n * n * sizeof(float));
		for (int i = 0; i < n; ++i)
		{
			mat[i] = tempMat + i * n;
			ans[i] = tempAns + i * n;
		}
	}
	else {
		tempMat = (float*)malloc(n * (n / comm_sz + 1) * sizeof(float));
		for (int i = n / comm_sz; i >= 0; i--)
			mat[i] = tempMat + i * n;
	}
	pivotRow = (float*)malloc(n * sizeof(float));//该一维向量用于消元时传递枢轴所在行

	for (int funIndex = 0; funIndex < funNum; ++funIndex)
	{
		if (my_rank == 0) std::cout << name[funIndex] << '\t';
		for (int threadIndex = 0; threadIndex < threadArrSize; ++threadIndex)
		{
			int thread = threadNumArr[threadIndex];
			double totalTime = 0.0f;
			
			for (int step = 0; step < REPEAT_NUM; ++step)
			{
				if (my_rank == 0)dataGenerate();
				MPI_Barrier(MPI_COMM_WORLD);//其他进程等待数据初始化

				//程序的主体部分
				double localTime = funArr[funIndex](thread);

				if (my_rank == 0)
				{
					check(step);
					totalTime += localTime;
				}
			}

			if (my_rank == 0) std::cout << totalTime / (double)REPEAT_NUM << "\t";
		}
		if (my_rank == 0) std::cout << '\n';
	}

	std::free(pivotRow);
	std::free(tempMat);
	if (my_rank == 0) std::free(tempAns);
}

//串行高斯消去
double serialGauss(int threadNum= THREADNUM)
{
	if (my_rank != 0)//仅使用0号进程进行计算
		return 0;

	double start_time, end_time;
	start_time = MPI_Wtime();

	for (int k = 0; k < n; ++k)
	{
			
		//对行进行归一化
		float pivot = mat[k][k];
		for (int j = k + 1; j < n; ++j)mat[k][j] /= pivot;
		mat[k][k] = 1.0f;

		//对右下角的矩阵进行消元
		for (int i = k + 1; i < n; ++i)
		{
			float pivot = mat[i][k];
			for (int j = k + 1; j < n; ++j)
				mat[i][j] -= mat[k][j] * pivot;
			mat[i][k] = 0.0f;
		}
	}

	end_time = MPI_Wtime();
	return end_time - start_time;
}

// MPI循环划分并行算法
double cycleGauss(int threadNum= THREADNUM) {
	double start_time, end_time;
	MPI_Barrier(MPI_COMM_WORLD);//所有线程同时开始，便于计时
	start_time = MPI_Wtime();

	//计算分配的任务数，该方法将最后不能除尽而剩下的任务也均匀地划分给了前几个进程
	int task_num = my_rank < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
	
	//首先将每个进程要处理的数据放到数组最前面
	if (my_rank == 0) {// 0号进程负责任务的初始分发工作
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//每个进程最高处理n / comm_sz + 1行数据
		for (int p = 1; p < comm_sz; p++) {
			for (int i = p; i < n; i += comm_sz) {
				memcpy(buff + i / comm_sz * n, mat[i], n * sizeof(float));
			}
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Send(buff, count * n, MPI_FLOAT, p, 0, MPI_COMM_WORLD);//将数据发送给对应的进程
		}
		for (int i = 1; i < task_num; i++) memcpy(mat[i], mat[i * comm_sz], n * sizeof(float));//将零号线程要处理的数据移到最前面
		std::free(buff);
	}else {// 非0号进程负责任务的接收工作
		MPI_Recv(mat[0], task_num * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	
	// 做归一化操作
	for (int k = 0; k < n; k++) {
		// 如果除法操作是本进程负责的任务，并将除法结果广播
		if (k % comm_sz == my_rank) {
			int localIndexK = k / comm_sz, pivot = mat[localIndexK][k];
			for (int j = k + 1; j < n; j++) {
				mat[localIndexK][j] /= pivot;
			}
			mat[localIndexK][k] = 1;
			memcpy(pivotRow + k + 1, mat[localIndexK] + k + 1, (n - k - 1) * sizeof(float));
		}
		MPI_Bcast(pivotRow + k + 1, n - k - 1, MPI_FLOAT, k % comm_sz, MPI_COMM_WORLD);
		
		// 进行消元操作
		int end = my_rank > (k % comm_sz) ? k / comm_sz : k / comm_sz + 1;
		for (int i = task_num-1; i >= end; i--) {
			float pivot = mat[i][k];
			for (int j = k + 1; j < n; j++) {
				mat[i][j] -= pivot * pivotRow[j];
			}
			mat[i][k] = 0;
		}
	}

	//将分散在各个进程的答案聚合起来
	if (my_rank == 0) {
		for (int i = task_num - 1; i > 0; i--)
		{
			memcpy(mat[i * comm_sz], mat[i], n * sizeof(float));
		}
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//每个进程最高处理n / comm_sz + 1行数据
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Recv((void*)buff, count * n, MPI_FLOAT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//接受对应的进程的数据
			for (int i = p; i < n; i += comm_sz) {
				memcpy(mat[i], buff + i / comm_sz * n, n * sizeof(float));
			}
		}
		std::free(buff);
	}
	else {
		MPI_Send(mat[0], task_num * n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);//等待全部线程完成
	end_time = MPI_Wtime();
	return end_time - start_time;
}

// MPI循环划分+流水线式分发并行算法
double piplineGauss(int threadNum= THREADNUM) {
	double start_time, end_time;
	MPI_Barrier(MPI_COMM_WORLD);//所有线程同时开始，便于计时
	start_time = MPI_Wtime();

	//计算分配的任务数，该方法将最后不能除尽而剩下的任务也均匀地划分给了前几个进程
	int task_num = my_rank < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;

	//首先将每个进程要处理的数据放到数组最前面
	if (my_rank == 0) {// 0号进程负责任务的初始分发工作
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//每个进程最高处理n / comm_sz + 1行数据
		for (int p = 1; p < comm_sz; p++) {
			for (int i = p; i < n; i += comm_sz) {
				memcpy(buff + i / comm_sz * n, mat[i], n * sizeof(float));
			}
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Send(buff, count * n, MPI_FLOAT, p, 0, MPI_COMM_WORLD);//将数据发送给对应的进程
		}
		for (int i = 1; i < task_num; i++) memcpy(mat[i], mat[i * comm_sz], n * sizeof(float));//将零号线程要处理的数据移到最前面
		std::free(buff);
	}
	else {// 非0号进程负责任务的接收工作
		MPI_Recv(mat[0], task_num * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// 做归一化操作
	int pre_proc = (my_rank + (comm_sz - 1)) % comm_sz;
	int next_proc = (my_rank + 1) % comm_sz;
	for (int k = 0; k < n; k++) {
		// 如果除法操作是本进程负责的任务，并将除法结果广播
		if (k % comm_sz == my_rank) {
			int localIndexK = k / comm_sz, pivot = mat[localIndexK][k];
			for (int j = k + 1; j < n; j++) {
				mat[localIndexK][j] /= pivot;
			}
			mat[localIndexK][k] = 1;
			memcpy(pivotRow + k + 1, mat[localIndexK] + k + 1, (n - k - 1) * sizeof(float));
			MPI_Send(pivotRow + k + 1, n - k - 1, MPI_FLOAT, next_proc, 1, MPI_COMM_WORLD);
		}
		else {
			MPI_Recv(pivotRow + k + 1, n - k - 1, MPI_FLOAT, pre_proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (next_proc != k % comm_sz) {
				MPI_Send(pivotRow + k + 1, n - k - 1, MPI_FLOAT, next_proc, 1, MPI_COMM_WORLD);
			}
		}

		// 进行消元操作
		int end = my_rank > (k % comm_sz) ? k / comm_sz : k / comm_sz + 1;
		for (int i = task_num - 1; i >= end; i--) {
			float pivot = mat[i][k];
			for (int j = k + 1; j < n; j++) {
				mat[i][j] -= pivot * pivotRow[j];
			}
			mat[i][k] = 0;
		}
	}

	//将分散在各个进程的答案聚合起来
	if (my_rank == 0) {
		for (int i = task_num - 1; i > 0; i--)
		{
			memcpy(mat[i * comm_sz], mat[i], n * sizeof(float));
		}
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//每个进程最高处理n / comm_sz + 1行数据
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Recv((void*)buff, count * n, MPI_FLOAT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//接受对应的进程的数据
			for (int i = p; i < n; i += comm_sz) {
				memcpy(mat[i], buff + i / comm_sz * n, n * sizeof(float));
			}
		}
		std::free(buff);
	}
	else {
		MPI_Send(mat[0], task_num * n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);//等待全部线程完成
	end_time = MPI_Wtime();
	return end_time - start_time;
}

// MPI循环划分+循环分发并行算法
double one2MulGauss(int threadNum= THREADNUM) {
	double start_time, end_time;
	MPI_Barrier(MPI_COMM_WORLD);//所有线程同时开始，便于计时
	start_time = MPI_Wtime();

	//计算分配的任务数，该方法将最后不能除尽而剩下的任务也均匀地划分给了前几个进程
	int task_num = my_rank < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;

	//首先将每个进程要处理的数据放到数组最前面
	if (my_rank == 0) {// 0号进程负责任务的初始分发工作
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//每个进程最高处理n / comm_sz + 1行数据
		for (int p = 1; p < comm_sz; p++) {
			for (int i = p; i < n; i += comm_sz) {
				memcpy(buff + i / comm_sz * n, mat[i], n * sizeof(float));
			}
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Send(buff, count * n, MPI_FLOAT, p, 0, MPI_COMM_WORLD);//将数据发送给对应的进程
		}
		for (int i = 1; i < task_num; i++) memcpy(mat[i], mat[i * comm_sz], n * sizeof(float));//将零号线程要处理的数据移到最前面
		std::free(buff);
	}
	else {// 非0号进程负责任务的接收工作
		MPI_Recv(mat[0], task_num * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// 做归一化操作
	for (int k = 0; k < n; k++) {
		// 如果除法操作是本进程负责的任务，则进行计算，并将除法结果广播
		if (k % comm_sz == my_rank) {
			int localIndexK = k / comm_sz, pivot = mat[localIndexK][k];
			for (int j = k + 1; j < n; j++) {
				mat[localIndexK][j] /= pivot;
			}
			mat[localIndexK][k] = 1;
			memcpy(pivotRow + k + 1, mat[localIndexK] + k + 1, (n - k - 1) * sizeof(float));
			for (int p = 0; p < comm_sz; p++) {
				if (p != my_rank) {
					MPI_Send(pivotRow + k + 1, n - k - 1, MPI_FLOAT, p, 1, MPI_COMM_WORLD);
				}
			}
		} else {// 其余进程接收除法行的结果
			MPI_Recv(pivotRow + k + 1, n - k - 1, MPI_FLOAT, k % comm_sz, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		// 进行消元操作
		int end = my_rank > (k % comm_sz) ? k / comm_sz : k / comm_sz + 1;
		for (int i = task_num - 1; i >= end; i--) {
			float pivot = mat[i][k];
			for (int j = k + 1; j < n; j++) {
				mat[i][j] -= pivot * pivotRow[j];
			}
			mat[i][k] = 0;
		}
	}

	//将分散在各个进程的答案聚合起来
	if (my_rank == 0) {
		for (int i = task_num - 1; i > 0; i--)
		{
			memcpy(mat[i * comm_sz], mat[i], n * sizeof(float));
		}
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//每个进程最高处理n / comm_sz + 1行数据
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Recv((void*)buff, count * n, MPI_FLOAT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//接受对应的进程的数据
			for (int i = p; i < n; i += comm_sz) {
				memcpy(mat[i], buff + i / comm_sz * n, n * sizeof(float));
			}
		}
		std::free(buff);
	}
	else {
		MPI_Send(mat[0], task_num * n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);//等待全部线程完成
	end_time = MPI_Wtime();
	return end_time - start_time;
}

// MPI块划分并行算法
double blockGauss(int threadNum= THREADNUM) {
	double start_time, end_time;
	MPI_Barrier(MPI_COMM_WORLD);//所有线程同时开始，便于计时
	start_time = MPI_Wtime();

	//计算每个进程分配的任务数、开始、结束等
	int* startArr = (int*)malloc(comm_sz * sizeof(int));
	startArr[0] = 0;
	for (int i = 1; i < comm_sz; i++) startArr[i] = startArr[i - 1] + (i-1 < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz);
	int start = startArr[my_rank];
	int task_num = my_rank < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
	int end = start + task_num;
	
	//首先将每个进程要处理的数据放到数组最前面
	if (my_rank == 0) {// 0号进程负责任务的初始分发工作
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Send(mat[startArr[p]], count * n, MPI_FLOAT, p, 0, MPI_COMM_WORLD);
		}
	}else {// 非0号进程负责任务的接收工作
		MPI_Recv(mat[0], task_num * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// 做归一化操作
	for (int k = 0; k < n; k++) {
		// 如果除法操作是本进程负责的任务则进行归一化，并将结果广播
		if (start <= k && k < end) {
			int localIndexK = k - start, pivot = mat[localIndexK][k];
			for (int j = k + 1; j < n; j++) {
				mat[localIndexK][j] /= pivot;
			}
			mat[localIndexK][k] = 1;
			memcpy(pivotRow + k + 1, mat[localIndexK] + k + 1, (n - k - 1) * sizeof(float));
		}
		int p, temp1 = n / comm_sz + 1, temp2 = n / comm_sz;
		p = (k < temp1* (n% comm_sz) ? k / temp1 : (k - (n % comm_sz)) / temp2);
		MPI_Bcast(pivotRow + k + 1, n - k - 1, MPI_FLOAT, p, MPI_COMM_WORLD);
		
		// 进行消元操作
		int startI = (k + 1 < start ? 0 : k + 1 - start);
		for (int i = startI; i < task_num; i++) {
			float pivot = mat[i][k];
			for (int j = k + 1; j < n; j++) {
				mat[i][j] -= pivot * pivotRow[j];
			}
			mat[i][k] = 0;
		}
	}

	//将分散在各个进程的答案聚合起来
	if (my_rank == 0) {
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Recv(mat[startArr[p]], count * n, MPI_FLOAT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	} else {
		MPI_Send(mat[0], task_num * n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}

	//释放堆中的数据
	std::free(startArr);

	MPI_Barrier(MPI_COMM_WORLD);//等待全部线程完成
	end_time = MPI_Wtime();
	return end_time - start_time;
}

// MPI循环划分+OMP并行算法
double mpiOmpGauss(int threadNum= THREADNUM) {
	double start_time, end_time;
	MPI_Barrier(MPI_COMM_WORLD);//所有线程同时开始，便于计时
	start_time = MPI_Wtime();

	//计算分配的任务数，该方法将最后不能除尽而剩下的任务也均匀地划分给了前几个进程
	int task_num = my_rank < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;

	//首先将每个进程要处理的数据放到数组最前面
	if (my_rank == 0) {// 0号进程负责任务的初始分发工作
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//每个进程最高处理n / comm_sz + 1行数据
		for (int p = 1; p < comm_sz; p++) {
			for (int i = p; i < n; i += comm_sz) {
				memcpy(buff + i / comm_sz * n, mat[i], n * sizeof(float));
			}
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Send(buff, count * n, MPI_FLOAT, p, 0, MPI_COMM_WORLD);//将数据发送给对应的进程
		}
		for (int i = 1; i < task_num; i++) memcpy(mat[i], mat[i * comm_sz], n * sizeof(float));//将零号线程要处理的数据移到最前面
		std::free(buff);
	}
	else {// 非0号进程负责任务的接收工作
		MPI_Recv(mat[0], task_num * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// 做归一化操作
	int i, j, k;
	float pivot;
#pragma omp parallel num_threads(threadNum) default(none) private(i, j, k, pivot, pivotRow, my_rank,task_num) shared(mat, n, comm_sz)
	for (k = 0; k < n; k++) {
		// 如果除法操作是本进程负责的任务，并将除法结果广播
		if (k % comm_sz == my_rank) {
			int localIndexK = k / comm_sz;
			pivot = mat[localIndexK][k];
#pragma omp for
			for (j = k + 1; j < n; j++)
			{
				mat[localIndexK][j] /= pivot;
			}
#pragma omp single
			{
				mat[localIndexK][k] = 1;
				memcpy(pivotRow + k + 1, mat[localIndexK] + k + 1, (n - k - 1) * sizeof(float));
			}	
		}
#pragma omp single
		MPI_Bcast(pivotRow + k + 1, n - k - 1, MPI_FLOAT, k % comm_sz, MPI_COMM_WORLD);

		int end = my_rank > (k % comm_sz) ? k / comm_sz : k / comm_sz + 1;
#pragma omp for schedule(dynamic, 1)
		for (i = task_num - 1; i >= end; i--) {
			pivot = mat[i][k];
			for (j = k + 1; j < n; j++) {
				mat[i][j] -= pivot * pivotRow[j];
			}
			mat[i][k] = 0;
		}
	}

	//将分散在各个进程的答案聚合起来
	if (my_rank == 0) {
		for (int i = task_num - 1; i > 0; i--)
		{
			memcpy(mat[i * comm_sz], mat[i], n * sizeof(float));
		}
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//每个进程最高处理n / comm_sz + 1行数据
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Recv((void*)buff, count * n, MPI_FLOAT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//接受对应的进程的数据
			for (int i = p; i < n; i += comm_sz) {
				memcpy(mat[i], buff + i / comm_sz * n, n * sizeof(float));
			}
		}
		std::free(buff);
	}
	else {
		MPI_Send(mat[0], task_num * n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);//等待全部线程完成
	end_time = MPI_Wtime();
	return end_time - start_time;
}

#ifndef __ARM_NEON
//MPI循环划分+SSE并行算法
double mpiSseGauss(int threadNum = THREADNUM)
{
	double start_time, end_time;
	__m128 vt, va, vaik, vakj, vaij, vamul;
	MPI_Barrier(MPI_COMM_WORLD);//所有线程同时开始，便于计时
	start_time = MPI_Wtime();

	//计算分配的任务数，该方法将最后不能除尽而剩下的任务也均匀地划分给了前几个进程
	int task_num = my_rank < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;

	//首先将每个进程要处理的数据放到数组最前面
	if (my_rank == 0) {// 0号进程负责任务的初始分发工作
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//每个进程最高处理n / comm_sz + 1行数据
		for (int p = 1; p < comm_sz; p++) {
			for (int i = p; i < n; i += comm_sz) {
				memcpy(buff + i / comm_sz * n, mat[i], n * sizeof(float));
			}
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Send(buff, count * n, MPI_FLOAT, p, 0, MPI_COMM_WORLD);//将数据发送给对应的进程
		}
		for (int i = 1; i < task_num; i++) memcpy(mat[i], mat[i * comm_sz], n * sizeof(float));//将零号线程要处理的数据移到最前面
		std::free(buff);
	}
	else {// 非0号进程负责任务的接收工作
		MPI_Recv(mat[0], task_num * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// 做归一化操作
	for (int k = 0; k < n; k++) {
		// 如果除法操作是本进程负责的任务，并将除法结果广播
		if (k % comm_sz == my_rank) {
			int localIndexK = k / comm_sz;
			vt = _mm_set1_ps(mat[localIndexK][k]);
			int j = k + 1;
			for (; j < n; j++)
			{
				//16字节对齐
				if (!((ull(mat[localIndexK] + j)) & 0xf)) break;
				mat[localIndexK][j] /= mat[localIndexK][k];
			}
			for (; j + 4 <= n; j += 4)
			{
				va = _mm_load_ps(mat[localIndexK] + j);
				va = _mm_div_ps(va, vt);
				_mm_store_ps(mat[localIndexK] + j, va);
			}
			for (; j < n; j++) mat[localIndexK][j] /= mat[localIndexK][k];
			mat[localIndexK][k] = 1;
			memcpy(pivotRow + k + 1, mat[localIndexK] + k + 1, (n - k - 1) * sizeof(float));
		}
		MPI_Bcast(pivotRow + k + 1, n - k - 1, MPI_FLOAT, k % comm_sz, MPI_COMM_WORLD);

		// 进行消元操作
		int end = my_rank > (k % comm_sz) ? k / comm_sz : k / comm_sz + 1;
		for (int i = task_num - 1; i >= end; i--) {
			vaik = _mm_set1_ps(mat[i][k]);
			int j = k + 1;
			for (; j < n; ++j)
			{
				//16字节对齐
				if (!((ull(mat[i] + j)) & 0xf)) break;
				mat[i][j] -= pivotRow[j] * mat[i][k];
			}
			for (; j + 4 <= n; j += 4)
			{
				vakj = _mm_load_ps(pivotRow + j);
				vaij = _mm_load_ps(mat[i] + j);
				vamul = _mm_mul_ps(vakj, vaik);
				vaij = _mm_sub_ps(vaij, vamul);
				_mm_store_ps(mat[i] + j, vaij);
			}
			for (; j < n; j++) mat[i][j] -= pivotRow[j] * mat[i][k];
			mat[i][k] = 0;
		}
	}

	//将分散在各个进程的答案聚合起来
	if (my_rank == 0) {
		for (int i = task_num - 1; i > 0; i--)
			memcpy(mat[i * comm_sz], mat[i], n * sizeof(float));
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//每个进程最高处理n / comm_sz + 1行数据
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Recv((void*)buff, count * n, MPI_FLOAT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//接受对应的进程的数据
			for (int i = p; i < n; i += comm_sz) {
				memcpy(mat[i], buff + i / comm_sz * n, n * sizeof(float));
			}
		}
		std::free(buff);
	}
	else {
		MPI_Send(mat[0], task_num * n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);//等待全部线程完成
	end_time = MPI_Wtime();
	return end_time - start_time;
}

// MPI循环划分+OMP+SSE并行算法
double mpiSseOmpGauss(int threadNum= THREADNUM) {
	double start_time, end_time;
	MPI_Barrier(MPI_COMM_WORLD);//所有线程同时开始，便于计时
	start_time = MPI_Wtime();

	//计算分配的任务数，该方法将最后不能除尽而剩下的任务也均匀地划分给了前几个进程
	int task_num = my_rank < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;

	//首先将每个进程要处理的数据放到数组最前面
	if (my_rank == 0) {// 0号进程负责任务的初始分发工作
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//每个进程最高处理n / comm_sz + 1行数据
		for (int p = 1; p < comm_sz; p++) {
			for (int i = p; i < n; i += comm_sz) {
				memcpy(buff + i / comm_sz * n, mat[i], n * sizeof(float));
			}
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Send(buff, count * n, MPI_FLOAT, p, 0, MPI_COMM_WORLD);//将数据发送给对应的进程
		}
		for (int i = 1; i < task_num; i++) memcpy(mat[i], mat[i * comm_sz], n * sizeof(float));//将零号线程要处理的数据移到最前面
		std::free(buff);
	}
	else {// 非0号进程负责任务的接收工作
		MPI_Recv(mat[0], task_num * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// 做归一化操作
	int i, j, k;
	float pivot;
	__m128 vt, va, vaik, vakj, vaij, vamul;
#pragma omp parallel num_threads(threadNum) default(none) private(i, j, k, pivot, vt, va, vaik, vakj, vaij, vamul, pivotRow, my_rank, task_num) shared(mat, n, comm_sz)
	for (k = 0; k < n; k++) {
		// 如果除法操作是本进程负责的任务，并将除法结果广播
		if (k % comm_sz == my_rank) {
			int localIndexK = k / comm_sz;
			pivot = mat[localIndexK][k];
			vt = _mm_set1_ps(mat[localIndexK][k]);
#pragma omp for
			for (j = k + 1; j <= n - 4; j += 4)
			{
				va = _mm_loadu_ps(mat[localIndexK] + j);
				va = _mm_div_ps(va, vt);
				_mm_storeu_ps(mat[localIndexK] + j, va);
			}
#pragma omp single
			{
				for (j = n - ((n - k - 1) % 4); j < n; j++)
				{
					mat[localIndexK][j] /= pivot;
				}
				mat[localIndexK][k] = 1;
				memcpy(pivotRow + k + 1, mat[localIndexK] + k + 1, (n - k - 1) * sizeof(float));
			}
		}
#pragma omp single
		MPI_Bcast(pivotRow + k + 1, n - k - 1, MPI_FLOAT, k % comm_sz, MPI_COMM_WORLD);

		//消元
		int end = my_rank > (k % comm_sz) ? k / comm_sz : k / comm_sz + 1;
#pragma omp for schedule(dynamic, 1)
		for (i = task_num - 1; i >= end; i--) {
			vaik = _mm_set1_ps(mat[i][k]);
			j = k + 1;
			for (; j + 4 <= n; j += 4)
			{
				vakj = _mm_load_ps(pivotRow + j);
				vaij = _mm_load_ps(mat[i] + j);
				vamul = _mm_mul_ps(vakj, vaik);
				vaij = _mm_sub_ps(vaij, vamul);
				_mm_store_ps(mat[i] + j, vaij);
			}
			for (; j < n; j++) mat[i][j] -= pivotRow[j] * mat[i][k];
			mat[i][k] = 0;
		}
	}

	//将分散在各个进程的答案聚合起来
	if (my_rank == 0) {
		for (int i = task_num - 1; i > 0; i--)
		{
			memcpy(mat[i * comm_sz], mat[i], n * sizeof(float));
		}
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//每个进程最高处理n / comm_sz + 1行数据
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Recv((void*)buff, count * n, MPI_FLOAT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//接受对应的进程的数据
			for (int i = p; i < n; i += comm_sz) {
				memcpy(mat[i], buff + i / comm_sz * n, n * sizeof(float));
			}
		}
		std::free(buff);
	}
	else {
		MPI_Send(mat[0], task_num * n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);//等待全部线程完成
	end_time = MPI_Wtime();
	return end_time - start_time;
}

//MPI循环划分+AVX2并行算法
double mpiAvx2Gauss(int threadNum = THREADNUM)
{
	double start_time, end_time;
	__m256 vt, va, vaik, vakj, vaij;
	MPI_Barrier(MPI_COMM_WORLD);//所有线程同时开始，便于计时
	start_time = MPI_Wtime();

	//计算分配的任务数，该方法将最后不能除尽而剩下的任务也均匀地划分给了前几个进程
	int task_num = my_rank < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;

	//首先将每个进程要处理的数据放到数组最前面
	if (my_rank == 0) {// 0号进程负责任务的初始分发工作
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//每个进程最高处理n / comm_sz + 1行数据
		for (int p = 1; p < comm_sz; p++) {
			for (int i = p; i < n; i += comm_sz) {
				memcpy(buff + i / comm_sz * n, mat[i], n * sizeof(float));
			}
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Send(buff, count * n, MPI_FLOAT, p, 0, MPI_COMM_WORLD);//将数据发送给对应的进程
		}
		for (int i = 1; i < task_num; i++) memcpy(mat[i], mat[i * comm_sz], n * sizeof(float));//将零号线程要处理的数据移到最前面
		std::free(buff);
	}
	else {// 非0号进程负责任务的接收工作
		MPI_Recv(mat[0], task_num * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// 做归一化操作
	for (int k = 0; k < n; k++) {
		// 如果除法操作是本进程负责的任务，并将除法结果广播
		if (k % comm_sz == my_rank) {
			int localIndexK = k / comm_sz;
			vt = _mm256_set1_ps(mat[localIndexK][k]);
			int j = k + 1;
			for (; j < n; j++)
			{
				//32字节对齐
				if (!((ull(mat[localIndexK] + j)) & 0x1f)) break;
				mat[localIndexK][j] /= mat[localIndexK][k];
			}
			for (; j + 8 <= n; j += 8)
			{
				va = _mm256_load_ps(mat[localIndexK] + j);
				va = _mm256_div_ps(va, vt);
				_mm256_store_ps(mat[localIndexK] + j, va);
			}
			for (; j < n; j++) mat[localIndexK][j] /= mat[localIndexK][k];
			mat[localIndexK][k] = 1;
			memcpy(pivotRow + k + 1, mat[localIndexK] + k + 1, (n - k - 1) * sizeof(float));
		}
		MPI_Bcast(pivotRow + k + 1, n - k - 1, MPI_FLOAT, k % comm_sz, MPI_COMM_WORLD);

		// 进行消元操作
		int end = my_rank > (k % comm_sz) ? k / comm_sz : k / comm_sz + 1;
		for (int i = task_num - 1; i >= end; i--) {
			vaik = _mm256_set1_ps(mat[i][k]);
			int j = k + 1;
			for (; j < n; ++j)
			{
				//32字节对齐
				if (!((ull(mat[i] + j)) & 0x1f)) break;
				mat[i][j] -= pivotRow[j] * mat[i][k];
			}
			for (; j + 8 <= n; j += 8)
			{
				vakj = _mm256_load_ps(pivotRow + j);
				vaij = _mm256_load_ps(mat[i] + j);
				vaij = _mm256_fnmadd_ps(vakj, vaik, vaij);
				_mm256_store_ps(mat[i] + j, vaij);
			}
			for (; j < n; j++) mat[i][j] -= pivotRow[j] * mat[i][k];
			mat[i][k] = 0;
		}
	}

	//将分散在各个进程的答案聚合起来
	if (my_rank == 0) {
		for (int i = task_num - 1; i > 0; i--)
			memcpy(mat[i * comm_sz], mat[i], n * sizeof(float));
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//每个进程最高处理n / comm_sz + 1行数据
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Recv((void*)buff, count * n, MPI_FLOAT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//接受对应的进程的数据
			for (int i = p; i < n; i += comm_sz) {
				memcpy(mat[i], buff + i / comm_sz * n, n * sizeof(float));
			}
		}
		std::free(buff);
	}
	else {
		MPI_Send(mat[0], task_num * n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);//等待全部线程完成
	end_time = MPI_Wtime();
	return end_time - start_time;
}

// MPI循环划分+OMP+AVX2并行算法
double mpiAvx2OmpGauss(int threadNum = THREADNUM) {
	double start_time, end_time;
	MPI_Barrier(MPI_COMM_WORLD);//所有线程同时开始，便于计时
	start_time = MPI_Wtime();

	//计算分配的任务数，该方法将最后不能除尽而剩下的任务也均匀地划分给了前几个进程
	int task_num = my_rank < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;

	//首先将每个进程要处理的数据放到数组最前面
	if (my_rank == 0) {// 0号进程负责任务的初始分发工作
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//每个进程最高处理n / comm_sz + 1行数据
		for (int p = 1; p < comm_sz; p++) {
			for (int i = p; i < n; i += comm_sz) {
				memcpy(buff + i / comm_sz * n, mat[i], n * sizeof(float));
			}
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Send(buff, count * n, MPI_FLOAT, p, 0, MPI_COMM_WORLD);//将数据发送给对应的进程
		}
		for (int i = 1; i < task_num; i++) memcpy(mat[i], mat[i * comm_sz], n * sizeof(float));//将零号线程要处理的数据移到最前面
		std::free(buff);
	}
	else {// 非0号进程负责任务的接收工作
		MPI_Recv(mat[0], task_num * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// 做归一化操作
	int i, j, k;
	float pivot;
	__m256 vt, va, vaik, vakj, vaij;
#pragma omp parallel num_threads(threadNum) default(none) private(i, j, k, pivot, vt, va, vaik, vakj, vaij, pivotRow, my_rank, task_num) shared(mat, n, comm_sz)
	for (k = 0; k < n; k++) {
		// 如果除法操作是本进程负责的任务，并将除法结果广播
		if (k % comm_sz == my_rank) {
			int localIndexK = k / comm_sz;
			pivot = mat[localIndexK][k];
			vt = _mm256_set1_ps(mat[localIndexK][k]);
#pragma omp for
			for (j = k + 1; j <= n - 8; j += 8)
			{
				va = _mm256_loadu_ps(mat[localIndexK] + j);
				va = _mm256_div_ps(va, vt);
				_mm256_storeu_ps(mat[localIndexK] + j, va);
			}
#pragma omp single
			{
				for (j = n - ((n - k - 1) % 8); j < n; j++)
				{
					mat[localIndexK][j] /= pivot;
				}
				mat[localIndexK][k] = 1;
				memcpy(pivotRow + k + 1, mat[localIndexK] + k + 1, (n - k - 1) * sizeof(float));
			}
		}
#pragma omp single
		MPI_Bcast(pivotRow + k + 1, n - k - 1, MPI_FLOAT, k % comm_sz, MPI_COMM_WORLD);

		//消元
		int end = my_rank > (k % comm_sz) ? k / comm_sz : k / comm_sz + 1;
#pragma omp for schedule(dynamic, 1)
		for (i = task_num - 1; i >= end; i--) {
			vaik = _mm256_set1_ps(mat[i][k]);
			int j = k + 1;
			for (; j + 8 <= n; j += 8)
			{
				vakj = _mm256_loadu_ps(pivotRow + j);
				vaij = _mm256_loadu_ps(mat[i] + j);
				vaij = _mm256_fnmadd_ps(vakj, vaik, vaij);
				_mm256_storeu_ps(mat[i] + j, vaij);
			}
			for (; j < n; j++) mat[i][j] -= pivotRow[j] * mat[i][k];
			mat[i][k] = 0;
		}
	}

	//将分散在各个进程的答案聚合起来
	if (my_rank == 0) {
		for (int i = task_num - 1; i > 0; i--)
		{
			memcpy(mat[i * comm_sz], mat[i], n * sizeof(float));
		}
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//每个进程最高处理n / comm_sz + 1行数据
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Recv((void*)buff, count * n, MPI_FLOAT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//接受对应的进程的数据
			for (int i = p; i < n; i += comm_sz) {
				memcpy(mat[i], buff + i / comm_sz * n, n * sizeof(float));
			}
		}
		std::free(buff);
	}
	else {
		MPI_Send(mat[0], task_num * n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);//等待全部线程完成
	end_time = MPI_Wtime();
	return end_time - start_time;
}
#else
//MPI循环划分+Neon并行算法
double mpiNeonGauss(int threadNum = THREADNUM)
{
	double start_time, end_time;
	float32x4_t inv_vt, va, vaik, vakj, vaij;
	MPI_Barrier(MPI_COMM_WORLD);//所有线程同时开始，便于计时
	start_time = MPI_Wtime();

	//计算分配的任务数，该方法将最后不能除尽而剩下的任务也均匀地划分给了前几个进程
	int task_num = my_rank < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;

	//首先将每个进程要处理的数据放到数组最前面
	if (my_rank == 0) {// 0号进程负责任务的初始分发工作
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//每个进程最高处理n / comm_sz + 1行数据
		for (int p = 1; p < comm_sz; p++) {
			for (int i = p; i < n; i += comm_sz) {
				memcpy(buff + i / comm_sz * n, mat[i], n * sizeof(float));
			}
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Send(buff, count * n, MPI_FLOAT, p, 0, MPI_COMM_WORLD);//将数据发送给对应的进程
		}
		for (int i = 1; i < task_num; i++) memcpy(mat[i], mat[i * comm_sz], n * sizeof(float));//将零号线程要处理的数据移到最前面
		std::free(buff);
	}
	else {// 非0号进程负责任务的接收工作
		MPI_Recv(mat[0], task_num * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// 做归一化操作
	for (int k = 0; k < n; k++) {
		// 如果除法操作是本进程负责的任务，并将除法结果广播
		if (k % comm_sz == my_rank) {
			int localIndexK = k / comm_sz;
			inv_vt = vdupq_n_f32(1.0f / mat[localIndexK][k]);
			int j = k + 1;
			for (; j < n; j++)
			{
				//16字节对齐
				if (!((ull(mat[localIndexK] + j)) & 0xf)) break;
				mat[localIndexK][j] /= mat[localIndexK][k];
			}
			for (; j + 4 <= n; j += 4)
			{
				va = vld1q_f32(mat[localIndexK] + j);
				va = vmulq_f32(va, inv_vt);
				vst1q_f32(mat[localIndexK] + j, va);
			}
			for (; j < n; j++) mat[localIndexK][j] /= mat[localIndexK][k];
			mat[localIndexK][k] = 1;
			memcpy(pivotRow + k + 1, mat[localIndexK] + k + 1, (n - k - 1) * sizeof(float));
		}
		MPI_Bcast(pivotRow + k + 1, n - k - 1, MPI_FLOAT, k % comm_sz, MPI_COMM_WORLD);

		// 进行消元操作
		int end = my_rank > (k % comm_sz) ? k / comm_sz : k / comm_sz + 1;
		for (int i = task_num - 1; i >= end; i--) {
			vaik = vdupq_n_f32(mat[i][k]);
			int j = k + 1;
			for (; j < n; ++j)
			{
				//16字节对齐
				if (!((ull(mat[i] + j)) & 0xf)) break;
				mat[i][j] -= pivotRow[j] * mat[i][k];
			}
			for (; j + 4 <= n; j += 4)
			{
				vakj = vld1q_f32(pivotRow + j);
				vaij = vld1q_f32(mat[i] + j);
				vaij = vmlsq_f32(vaij, vakj, vaik);
				vst1q_f32(mat[i] + j, vaij);
			}
			for (; j < n; j++) mat[i][j] -= pivotRow[j] * mat[i][k];
			mat[i][k] = 0;
		}
	}

	//将分散在各个进程的答案聚合起来
	if (my_rank == 0) {
		for (int i = task_num - 1; i > 0; i--)
			memcpy(mat[i * comm_sz], mat[i], n * sizeof(float));
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//每个进程最高处理n / comm_sz + 1行数据
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Recv((void*)buff, count * n, MPI_FLOAT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//接受对应的进程的数据
			for (int i = p; i < n; i += comm_sz) {
				memcpy(mat[i], buff + i / comm_sz * n, n * sizeof(float));
			}
		}
		std::free(buff);
	}
	else {
		MPI_Send(mat[0], task_num * n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);//等待全部线程完成
	end_time = MPI_Wtime();
	return end_time - start_time;
}

// MPI循环划分+OMP+Neon并行算法
double mpiNeonOmpGauss(int threadNum= THREADNUM) {
	double start_time, end_time;
	MPI_Barrier(MPI_COMM_WORLD);//所有线程同时开始，便于计时
	start_time = MPI_Wtime();

	//计算分配的任务数，该方法将最后不能除尽而剩下的任务也均匀地划分给了前几个进程
	int task_num = my_rank < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;

	//首先将每个进程要处理的数据放到数组最前面
	if (my_rank == 0) {// 0号进程负责任务的初始分发工作
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//每个进程最高处理n / comm_sz + 1行数据
		for (int p = 1; p < comm_sz; p++) {
			for (int i = p; i < n; i += comm_sz) {
				memcpy(buff + i / comm_sz * n, mat[i], n * sizeof(float));
			}
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Send(buff, count * n, MPI_FLOAT, p, 0, MPI_COMM_WORLD);//将数据发送给对应的进程
		}
		for (int i = 1; i < task_num; i++) memcpy(mat[i], mat[i * comm_sz], n * sizeof(float));//将零号线程要处理的数据移到最前面
		std::free(buff);
	}
	else {// 非0号进程负责任务的接收工作
		MPI_Recv(mat[0], task_num * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// 做归一化操作
	int i, j, k;
	float pivot;
	float32x4_t inv_vt, va, vaik, vakj, vaij;
#pragma omp parallel num_threads(threadNum) default(none) private(i, j, k, pivot, inv_vt, va, vaik, vakj, vaij,pivotRow, my_rank, task_num) shared(mat, n, comm_sz)
	for (k = 0; k < n; k++) {
		// 如果除法操作是本进程负责的任务，并将除法结果广播
		if (k % comm_sz == my_rank) {
			int localIndexK = k / comm_sz;
			pivot = mat[localIndexK][k];
			inv_vt = vdupq_n_f32(1.0f / mat[localIndexK][k]);
#pragma omp for
			for (j = k + 1; j <= n - 8; j += 8)
			{
				va = vld1q_f32(mat[localIndexK] + j);
				va = vmulq_f32(va, inv_vt);
				vst1q_f32(mat[localIndexK] + j, va);
			}
#pragma omp single
			{
				for (j = n - ((n - k - 1) % 8); j < n; j++)
				{
					mat[localIndexK][j] /= pivot;
				}
				mat[localIndexK][k] = 1;
				memcpy(pivotRow + k + 1, mat[localIndexK] + k + 1, (n - k - 1) * sizeof(float));
			}
		}
#pragma omp single
		MPI_Bcast(pivotRow + k + 1, n - k - 1, MPI_FLOAT, k % comm_sz, MPI_COMM_WORLD);

		int end = my_rank > (k % comm_sz) ? k / comm_sz : k / comm_sz + 1;
#pragma omp for schedule(dynamic, 1)
		for (i = task_num - 1; i >= end; i--) {
			vaik = vdupq_n_f32(mat[i][k]);
			j = k + 1;
			for (; j + 4 <= n; j += 4)
			{
				vakj = vld1q_f32(pivotRow + j);
				vaij = vld1q_f32(mat[i] + j);
				vaij = vmlsq_f32(vaij, vakj, vaik);
				vst1q_f32(mat[i] + j, vaij);
			}
			for (; j < n; j++) mat[i][j] -= pivotRow[j] * mat[i][k];
			mat[i][k] = 0.0;
		}
	}

	//将分散在各个进程的答案聚合起来
	if (my_rank == 0) {
		for (int i = task_num - 1; i > 0; i--)
		{
			memcpy(mat[i * comm_sz], mat[i], n * sizeof(float));
		}
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//每个进程最高处理n / comm_sz + 1行数据
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Recv((void*)buff, count * n, MPI_FLOAT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//接受对应的进程的数据
			for (int i = p; i < n; i += comm_sz) {
				memcpy(mat[i], buff + i / comm_sz * n, n * sizeof(float));
			}
		}
		std::free(buff);
	}
	else {
		MPI_Send(mat[0], task_num * n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);//等待全部线程完成
	end_time = MPI_Wtime();
	return end_time - start_time;
}
#endif

int main()
{
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	if (my_rank == 0)cout << "comm_sz=" << comm_sz << endl;

	/*double(*funArr1[])(int) = { serialGauss ,cycleGauss ,piplineGauss ,one2MulGauss ,blockGauss,mpiNeonGauss };
	string nameArr1[] = { "serialGauss", "cycleGauss" ,"piplineGauss","one2MulGauss","blockGauss","mpiNeonGauss"};
	timingAll(nameArr1, funArr1, 6);*/

	//double(*funArr2[])(int) = { serialGauss ,cycleGauss ,mpiSseGauss ,mpiAvx2Gauss,mpiAvx2OmpGauss };
	//string nameArr2[] = { "serialGauss", "cycleGauss" ,"mpiSseGauss","mpiAvx2Gauss","mpiAvx2OmpGauss" };
	//timingAll(nameArr2, funArr2, 5);

	double(*funArr3[])(int) = { mpiAvx2OmpGauss, };
	string nameArr3[] = { "mpiAvx2OmpGauss"};
	timingAllThreadNum(nameArr3, funArr3, 1, threadNumArr, 19);

	MPI_Finalize();
}