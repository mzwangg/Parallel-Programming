#define _CRT_SECURE_NO_WARNINGS
#include "config.h"

float* mat[MAXN], * ans[MAXN], * pivotRow, * tempMat, * tempAns;
int n, comm_sz, my_rank;

//��ӡ��������ڵ���
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

//�����
void check(int step)
{

	//�ж��Ƿ���Ҫ���
#ifndef CHECK
	return;
#endif

	//n̫��ʱ��ӡ�����ù۲죬�ʲ��ٴ�ӡ���
	if (n <= MAX_PRINT_N)
	{
		printMatrix("ans", ans[0], n, n);
		printMatrix("mat", mat[0], n, n);
	}

	//ֻҪmat���������ans��ͬ�����������ȷ
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

//Ϊ������ּ�����Ϊnaf������ȼ�����������������������Ǿ����ٽ����б任�ķ������ɲ�������
void dataGenerate()
{
	memset(mat[0], 0, n * n * sizeof(float));//��Ȼmat�Ƕ���ָ�룬����һ��ָ����ʵ�����ڵ�
	for (int i = 0; i < n; ++i)
	{
		mat[i][i] = 1.0f;
		for (int j = i + 1; j < n; ++j)
			mat[i][j] = rand() & 0xF;//ʹ������������0��15֮��,��ֹ��������б任��������ʱ�������
	}

	//��ǰ�Ľ����Ϊ���յĴ�
	memcpy(ans[0], mat[0], n * n * sizeof(float));//��Ȼans�Ƕ���ָ�룬����һ��ָ����ʵ�����ڵ�

	int inv_rate = n >> 3;
	bool reGenerate = false;
	while (true)
	{
		//��������������ݣ���mat����Ϊans
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

//���Ժ���
void timing(double(*fun)(int), int threadNum= THREADNUM)
{
	for (n = 64; n <= MAXN; n += 64)
	{
		double totalTime = 0.0f;
		//ֱ�ӽ�������ά��������Ϊһ�����������ڴ�
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
		pivotRow = (float*)malloc(n * sizeof(float));//��һά����������Ԫʱ��������������

		for (int step = 1; step <= REPEAT_NUM; ++step)
		{
			if (my_rank == 0)dataGenerate();
			MPI_Barrier(MPI_COMM_WORLD);//�������̵ȴ����ݳ�ʼ��

			//��������岿��
			double localTime = fun(threadNum);

			if (my_rank == 0)
			{
				check(step);
				totalTime += localTime;
			}
		}
		if (my_rank == 0) std::cout << totalTime / (double)REPEAT_NUM << "\t";

		//�ͷŶ��е�����
		std::free(pivotRow);
		std::free(tempMat);
		if (my_rank == 0) std::free(tempAns);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	if (my_rank == 0) std::cout << "\n";
}

//���Ժ���
void timingAll(string name[], double(*funArr[])(int), int funNum, int threadNum = THREADNUM)
{
	//�����ͷ��Ϣ
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
			//ֱ�ӽ�������ά��������Ϊһ�����������ڴ�
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
			pivotRow = (float*)malloc(n * sizeof(float));//��һά����������Ԫʱ��������������

			for (int step = 0; step < REPEAT_NUM; ++step)
			{
				if (my_rank == 0)dataGenerate();
				MPI_Barrier(MPI_COMM_WORLD);//�������̵ȴ����ݳ�ʼ��

				//��������岿��
				double localTime = funArr[funIndex](threadNum);

				if (my_rank == 0)
				{
					check(step);
					totalTime += localTime;
				}
			}

			//�ͷŶ��е�����
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

//���Ժ���
void timingAllThreadNum(string name[], double(*funArr[])(int), int funNum, const int threadArr[], int threadArrSize)
{
	//�����ͷ��Ϣ
	if (my_rank == 0) {
		std::cout << "Type/Thread\t";
		for (int threadIndex = 0; threadIndex < threadArrSize; ++threadIndex)
			std::cout << threadNumArr[threadIndex] << '\t';
		std::cout << "\n";
	}

	//�ڶ�����������
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
	pivotRow = (float*)malloc(n * sizeof(float));//��һά����������Ԫʱ��������������

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
				MPI_Barrier(MPI_COMM_WORLD);//�������̵ȴ����ݳ�ʼ��

				//��������岿��
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

//���и�˹��ȥ
double serialGauss(int threadNum= THREADNUM)
{
	if (my_rank != 0)//��ʹ��0�Ž��̽��м���
		return 0;

	double start_time, end_time;
	start_time = MPI_Wtime();

	for (int k = 0; k < n; ++k)
	{
			
		//���н��й�һ��
		float pivot = mat[k][k];
		for (int j = k + 1; j < n; ++j)mat[k][j] /= pivot;
		mat[k][k] = 1.0f;

		//�����½ǵľ��������Ԫ
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

// MPIѭ�����ֲ����㷨
double cycleGauss(int threadNum= THREADNUM) {
	double start_time, end_time;
	MPI_Barrier(MPI_COMM_WORLD);//�����߳�ͬʱ��ʼ�����ڼ�ʱ
	start_time = MPI_Wtime();

	//�����������������÷���������ܳ�����ʣ�µ�����Ҳ���ȵػ��ָ���ǰ��������
	int task_num = my_rank < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
	
	//���Ƚ�ÿ������Ҫ��������ݷŵ�������ǰ��
	if (my_rank == 0) {// 0�Ž��̸�������ĳ�ʼ�ַ�����
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//ÿ��������ߴ���n / comm_sz + 1������
		for (int p = 1; p < comm_sz; p++) {
			for (int i = p; i < n; i += comm_sz) {
				memcpy(buff + i / comm_sz * n, mat[i], n * sizeof(float));
			}
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Send(buff, count * n, MPI_FLOAT, p, 0, MPI_COMM_WORLD);//�����ݷ��͸���Ӧ�Ľ���
		}
		for (int i = 1; i < task_num; i++) memcpy(mat[i], mat[i * comm_sz], n * sizeof(float));//������߳�Ҫ����������Ƶ���ǰ��
		std::free(buff);
	}else {// ��0�Ž��̸�������Ľ��չ���
		MPI_Recv(mat[0], task_num * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	
	// ����һ������
	for (int k = 0; k < n; k++) {
		// ������������Ǳ����̸�������񣬲�����������㲥
		if (k % comm_sz == my_rank) {
			int localIndexK = k / comm_sz, pivot = mat[localIndexK][k];
			for (int j = k + 1; j < n; j++) {
				mat[localIndexK][j] /= pivot;
			}
			mat[localIndexK][k] = 1;
			memcpy(pivotRow + k + 1, mat[localIndexK] + k + 1, (n - k - 1) * sizeof(float));
		}
		MPI_Bcast(pivotRow + k + 1, n - k - 1, MPI_FLOAT, k % comm_sz, MPI_COMM_WORLD);
		
		// ������Ԫ����
		int end = my_rank > (k % comm_sz) ? k / comm_sz : k / comm_sz + 1;
		for (int i = task_num-1; i >= end; i--) {
			float pivot = mat[i][k];
			for (int j = k + 1; j < n; j++) {
				mat[i][j] -= pivot * pivotRow[j];
			}
			mat[i][k] = 0;
		}
	}

	//����ɢ�ڸ������̵Ĵ𰸾ۺ�����
	if (my_rank == 0) {
		for (int i = task_num - 1; i > 0; i--)
		{
			memcpy(mat[i * comm_sz], mat[i], n * sizeof(float));
		}
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//ÿ��������ߴ���n / comm_sz + 1������
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Recv((void*)buff, count * n, MPI_FLOAT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//���ܶ�Ӧ�Ľ��̵�����
			for (int i = p; i < n; i += comm_sz) {
				memcpy(mat[i], buff + i / comm_sz * n, n * sizeof(float));
			}
		}
		std::free(buff);
	}
	else {
		MPI_Send(mat[0], task_num * n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);//�ȴ�ȫ���߳����
	end_time = MPI_Wtime();
	return end_time - start_time;
}

// MPIѭ������+��ˮ��ʽ�ַ������㷨
double piplineGauss(int threadNum= THREADNUM) {
	double start_time, end_time;
	MPI_Barrier(MPI_COMM_WORLD);//�����߳�ͬʱ��ʼ�����ڼ�ʱ
	start_time = MPI_Wtime();

	//�����������������÷���������ܳ�����ʣ�µ�����Ҳ���ȵػ��ָ���ǰ��������
	int task_num = my_rank < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;

	//���Ƚ�ÿ������Ҫ��������ݷŵ�������ǰ��
	if (my_rank == 0) {// 0�Ž��̸�������ĳ�ʼ�ַ�����
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//ÿ��������ߴ���n / comm_sz + 1������
		for (int p = 1; p < comm_sz; p++) {
			for (int i = p; i < n; i += comm_sz) {
				memcpy(buff + i / comm_sz * n, mat[i], n * sizeof(float));
			}
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Send(buff, count * n, MPI_FLOAT, p, 0, MPI_COMM_WORLD);//�����ݷ��͸���Ӧ�Ľ���
		}
		for (int i = 1; i < task_num; i++) memcpy(mat[i], mat[i * comm_sz], n * sizeof(float));//������߳�Ҫ����������Ƶ���ǰ��
		std::free(buff);
	}
	else {// ��0�Ž��̸�������Ľ��չ���
		MPI_Recv(mat[0], task_num * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// ����һ������
	int pre_proc = (my_rank + (comm_sz - 1)) % comm_sz;
	int next_proc = (my_rank + 1) % comm_sz;
	for (int k = 0; k < n; k++) {
		// ������������Ǳ����̸�������񣬲�����������㲥
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

		// ������Ԫ����
		int end = my_rank > (k % comm_sz) ? k / comm_sz : k / comm_sz + 1;
		for (int i = task_num - 1; i >= end; i--) {
			float pivot = mat[i][k];
			for (int j = k + 1; j < n; j++) {
				mat[i][j] -= pivot * pivotRow[j];
			}
			mat[i][k] = 0;
		}
	}

	//����ɢ�ڸ������̵Ĵ𰸾ۺ�����
	if (my_rank == 0) {
		for (int i = task_num - 1; i > 0; i--)
		{
			memcpy(mat[i * comm_sz], mat[i], n * sizeof(float));
		}
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//ÿ��������ߴ���n / comm_sz + 1������
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Recv((void*)buff, count * n, MPI_FLOAT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//���ܶ�Ӧ�Ľ��̵�����
			for (int i = p; i < n; i += comm_sz) {
				memcpy(mat[i], buff + i / comm_sz * n, n * sizeof(float));
			}
		}
		std::free(buff);
	}
	else {
		MPI_Send(mat[0], task_num * n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);//�ȴ�ȫ���߳����
	end_time = MPI_Wtime();
	return end_time - start_time;
}

// MPIѭ������+ѭ���ַ������㷨
double one2MulGauss(int threadNum= THREADNUM) {
	double start_time, end_time;
	MPI_Barrier(MPI_COMM_WORLD);//�����߳�ͬʱ��ʼ�����ڼ�ʱ
	start_time = MPI_Wtime();

	//�����������������÷���������ܳ�����ʣ�µ�����Ҳ���ȵػ��ָ���ǰ��������
	int task_num = my_rank < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;

	//���Ƚ�ÿ������Ҫ��������ݷŵ�������ǰ��
	if (my_rank == 0) {// 0�Ž��̸�������ĳ�ʼ�ַ�����
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//ÿ��������ߴ���n / comm_sz + 1������
		for (int p = 1; p < comm_sz; p++) {
			for (int i = p; i < n; i += comm_sz) {
				memcpy(buff + i / comm_sz * n, mat[i], n * sizeof(float));
			}
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Send(buff, count * n, MPI_FLOAT, p, 0, MPI_COMM_WORLD);//�����ݷ��͸���Ӧ�Ľ���
		}
		for (int i = 1; i < task_num; i++) memcpy(mat[i], mat[i * comm_sz], n * sizeof(float));//������߳�Ҫ����������Ƶ���ǰ��
		std::free(buff);
	}
	else {// ��0�Ž��̸�������Ľ��չ���
		MPI_Recv(mat[0], task_num * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// ����һ������
	for (int k = 0; k < n; k++) {
		// ������������Ǳ����̸������������м��㣬������������㲥
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
		} else {// ������̽��ճ����еĽ��
			MPI_Recv(pivotRow + k + 1, n - k - 1, MPI_FLOAT, k % comm_sz, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		// ������Ԫ����
		int end = my_rank > (k % comm_sz) ? k / comm_sz : k / comm_sz + 1;
		for (int i = task_num - 1; i >= end; i--) {
			float pivot = mat[i][k];
			for (int j = k + 1; j < n; j++) {
				mat[i][j] -= pivot * pivotRow[j];
			}
			mat[i][k] = 0;
		}
	}

	//����ɢ�ڸ������̵Ĵ𰸾ۺ�����
	if (my_rank == 0) {
		for (int i = task_num - 1; i > 0; i--)
		{
			memcpy(mat[i * comm_sz], mat[i], n * sizeof(float));
		}
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//ÿ��������ߴ���n / comm_sz + 1������
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Recv((void*)buff, count * n, MPI_FLOAT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//���ܶ�Ӧ�Ľ��̵�����
			for (int i = p; i < n; i += comm_sz) {
				memcpy(mat[i], buff + i / comm_sz * n, n * sizeof(float));
			}
		}
		std::free(buff);
	}
	else {
		MPI_Send(mat[0], task_num * n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);//�ȴ�ȫ���߳����
	end_time = MPI_Wtime();
	return end_time - start_time;
}

// MPI�黮�ֲ����㷨
double blockGauss(int threadNum= THREADNUM) {
	double start_time, end_time;
	MPI_Barrier(MPI_COMM_WORLD);//�����߳�ͬʱ��ʼ�����ڼ�ʱ
	start_time = MPI_Wtime();

	//����ÿ�����̷��������������ʼ��������
	int* startArr = (int*)malloc(comm_sz * sizeof(int));
	startArr[0] = 0;
	for (int i = 1; i < comm_sz; i++) startArr[i] = startArr[i - 1] + (i-1 < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz);
	int start = startArr[my_rank];
	int task_num = my_rank < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
	int end = start + task_num;
	
	//���Ƚ�ÿ������Ҫ��������ݷŵ�������ǰ��
	if (my_rank == 0) {// 0�Ž��̸�������ĳ�ʼ�ַ�����
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Send(mat[startArr[p]], count * n, MPI_FLOAT, p, 0, MPI_COMM_WORLD);
		}
	}else {// ��0�Ž��̸�������Ľ��չ���
		MPI_Recv(mat[0], task_num * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// ����һ������
	for (int k = 0; k < n; k++) {
		// ������������Ǳ����̸������������й�һ������������㲥
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
		
		// ������Ԫ����
		int startI = (k + 1 < start ? 0 : k + 1 - start);
		for (int i = startI; i < task_num; i++) {
			float pivot = mat[i][k];
			for (int j = k + 1; j < n; j++) {
				mat[i][j] -= pivot * pivotRow[j];
			}
			mat[i][k] = 0;
		}
	}

	//����ɢ�ڸ������̵Ĵ𰸾ۺ�����
	if (my_rank == 0) {
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Recv(mat[startArr[p]], count * n, MPI_FLOAT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	} else {
		MPI_Send(mat[0], task_num * n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}

	//�ͷŶ��е�����
	std::free(startArr);

	MPI_Barrier(MPI_COMM_WORLD);//�ȴ�ȫ���߳����
	end_time = MPI_Wtime();
	return end_time - start_time;
}

// MPIѭ������+OMP�����㷨
double mpiOmpGauss(int threadNum= THREADNUM) {
	double start_time, end_time;
	MPI_Barrier(MPI_COMM_WORLD);//�����߳�ͬʱ��ʼ�����ڼ�ʱ
	start_time = MPI_Wtime();

	//�����������������÷���������ܳ�����ʣ�µ�����Ҳ���ȵػ��ָ���ǰ��������
	int task_num = my_rank < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;

	//���Ƚ�ÿ������Ҫ��������ݷŵ�������ǰ��
	if (my_rank == 0) {// 0�Ž��̸�������ĳ�ʼ�ַ�����
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//ÿ��������ߴ���n / comm_sz + 1������
		for (int p = 1; p < comm_sz; p++) {
			for (int i = p; i < n; i += comm_sz) {
				memcpy(buff + i / comm_sz * n, mat[i], n * sizeof(float));
			}
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Send(buff, count * n, MPI_FLOAT, p, 0, MPI_COMM_WORLD);//�����ݷ��͸���Ӧ�Ľ���
		}
		for (int i = 1; i < task_num; i++) memcpy(mat[i], mat[i * comm_sz], n * sizeof(float));//������߳�Ҫ����������Ƶ���ǰ��
		std::free(buff);
	}
	else {// ��0�Ž��̸�������Ľ��չ���
		MPI_Recv(mat[0], task_num * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// ����һ������
	int i, j, k;
	float pivot;
#pragma omp parallel num_threads(threadNum) default(none) private(i, j, k, pivot, pivotRow, my_rank,task_num) shared(mat, n, comm_sz)
	for (k = 0; k < n; k++) {
		// ������������Ǳ����̸�������񣬲�����������㲥
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

	//����ɢ�ڸ������̵Ĵ𰸾ۺ�����
	if (my_rank == 0) {
		for (int i = task_num - 1; i > 0; i--)
		{
			memcpy(mat[i * comm_sz], mat[i], n * sizeof(float));
		}
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//ÿ��������ߴ���n / comm_sz + 1������
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Recv((void*)buff, count * n, MPI_FLOAT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//���ܶ�Ӧ�Ľ��̵�����
			for (int i = p; i < n; i += comm_sz) {
				memcpy(mat[i], buff + i / comm_sz * n, n * sizeof(float));
			}
		}
		std::free(buff);
	}
	else {
		MPI_Send(mat[0], task_num * n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);//�ȴ�ȫ���߳����
	end_time = MPI_Wtime();
	return end_time - start_time;
}

#ifndef __ARM_NEON
//MPIѭ������+SSE�����㷨
double mpiSseGauss(int threadNum = THREADNUM)
{
	double start_time, end_time;
	__m128 vt, va, vaik, vakj, vaij, vamul;
	MPI_Barrier(MPI_COMM_WORLD);//�����߳�ͬʱ��ʼ�����ڼ�ʱ
	start_time = MPI_Wtime();

	//�����������������÷���������ܳ�����ʣ�µ�����Ҳ���ȵػ��ָ���ǰ��������
	int task_num = my_rank < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;

	//���Ƚ�ÿ������Ҫ��������ݷŵ�������ǰ��
	if (my_rank == 0) {// 0�Ž��̸�������ĳ�ʼ�ַ�����
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//ÿ��������ߴ���n / comm_sz + 1������
		for (int p = 1; p < comm_sz; p++) {
			for (int i = p; i < n; i += comm_sz) {
				memcpy(buff + i / comm_sz * n, mat[i], n * sizeof(float));
			}
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Send(buff, count * n, MPI_FLOAT, p, 0, MPI_COMM_WORLD);//�����ݷ��͸���Ӧ�Ľ���
		}
		for (int i = 1; i < task_num; i++) memcpy(mat[i], mat[i * comm_sz], n * sizeof(float));//������߳�Ҫ����������Ƶ���ǰ��
		std::free(buff);
	}
	else {// ��0�Ž��̸�������Ľ��չ���
		MPI_Recv(mat[0], task_num * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// ����һ������
	for (int k = 0; k < n; k++) {
		// ������������Ǳ����̸�������񣬲�����������㲥
		if (k % comm_sz == my_rank) {
			int localIndexK = k / comm_sz;
			vt = _mm_set1_ps(mat[localIndexK][k]);
			int j = k + 1;
			for (; j < n; j++)
			{
				//16�ֽڶ���
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

		// ������Ԫ����
		int end = my_rank > (k % comm_sz) ? k / comm_sz : k / comm_sz + 1;
		for (int i = task_num - 1; i >= end; i--) {
			vaik = _mm_set1_ps(mat[i][k]);
			int j = k + 1;
			for (; j < n; ++j)
			{
				//16�ֽڶ���
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

	//����ɢ�ڸ������̵Ĵ𰸾ۺ�����
	if (my_rank == 0) {
		for (int i = task_num - 1; i > 0; i--)
			memcpy(mat[i * comm_sz], mat[i], n * sizeof(float));
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//ÿ��������ߴ���n / comm_sz + 1������
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Recv((void*)buff, count * n, MPI_FLOAT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//���ܶ�Ӧ�Ľ��̵�����
			for (int i = p; i < n; i += comm_sz) {
				memcpy(mat[i], buff + i / comm_sz * n, n * sizeof(float));
			}
		}
		std::free(buff);
	}
	else {
		MPI_Send(mat[0], task_num * n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);//�ȴ�ȫ���߳����
	end_time = MPI_Wtime();
	return end_time - start_time;
}

// MPIѭ������+OMP+SSE�����㷨
double mpiSseOmpGauss(int threadNum= THREADNUM) {
	double start_time, end_time;
	MPI_Barrier(MPI_COMM_WORLD);//�����߳�ͬʱ��ʼ�����ڼ�ʱ
	start_time = MPI_Wtime();

	//�����������������÷���������ܳ�����ʣ�µ�����Ҳ���ȵػ��ָ���ǰ��������
	int task_num = my_rank < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;

	//���Ƚ�ÿ������Ҫ��������ݷŵ�������ǰ��
	if (my_rank == 0) {// 0�Ž��̸�������ĳ�ʼ�ַ�����
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//ÿ��������ߴ���n / comm_sz + 1������
		for (int p = 1; p < comm_sz; p++) {
			for (int i = p; i < n; i += comm_sz) {
				memcpy(buff + i / comm_sz * n, mat[i], n * sizeof(float));
			}
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Send(buff, count * n, MPI_FLOAT, p, 0, MPI_COMM_WORLD);//�����ݷ��͸���Ӧ�Ľ���
		}
		for (int i = 1; i < task_num; i++) memcpy(mat[i], mat[i * comm_sz], n * sizeof(float));//������߳�Ҫ����������Ƶ���ǰ��
		std::free(buff);
	}
	else {// ��0�Ž��̸�������Ľ��չ���
		MPI_Recv(mat[0], task_num * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// ����һ������
	int i, j, k;
	float pivot;
	__m128 vt, va, vaik, vakj, vaij, vamul;
#pragma omp parallel num_threads(threadNum) default(none) private(i, j, k, pivot, vt, va, vaik, vakj, vaij, vamul, pivotRow, my_rank, task_num) shared(mat, n, comm_sz)
	for (k = 0; k < n; k++) {
		// ������������Ǳ����̸�������񣬲�����������㲥
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

		//��Ԫ
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

	//����ɢ�ڸ������̵Ĵ𰸾ۺ�����
	if (my_rank == 0) {
		for (int i = task_num - 1; i > 0; i--)
		{
			memcpy(mat[i * comm_sz], mat[i], n * sizeof(float));
		}
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//ÿ��������ߴ���n / comm_sz + 1������
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Recv((void*)buff, count * n, MPI_FLOAT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//���ܶ�Ӧ�Ľ��̵�����
			for (int i = p; i < n; i += comm_sz) {
				memcpy(mat[i], buff + i / comm_sz * n, n * sizeof(float));
			}
		}
		std::free(buff);
	}
	else {
		MPI_Send(mat[0], task_num * n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);//�ȴ�ȫ���߳����
	end_time = MPI_Wtime();
	return end_time - start_time;
}

//MPIѭ������+AVX2�����㷨
double mpiAvx2Gauss(int threadNum = THREADNUM)
{
	double start_time, end_time;
	__m256 vt, va, vaik, vakj, vaij;
	MPI_Barrier(MPI_COMM_WORLD);//�����߳�ͬʱ��ʼ�����ڼ�ʱ
	start_time = MPI_Wtime();

	//�����������������÷���������ܳ�����ʣ�µ�����Ҳ���ȵػ��ָ���ǰ��������
	int task_num = my_rank < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;

	//���Ƚ�ÿ������Ҫ��������ݷŵ�������ǰ��
	if (my_rank == 0) {// 0�Ž��̸�������ĳ�ʼ�ַ�����
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//ÿ��������ߴ���n / comm_sz + 1������
		for (int p = 1; p < comm_sz; p++) {
			for (int i = p; i < n; i += comm_sz) {
				memcpy(buff + i / comm_sz * n, mat[i], n * sizeof(float));
			}
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Send(buff, count * n, MPI_FLOAT, p, 0, MPI_COMM_WORLD);//�����ݷ��͸���Ӧ�Ľ���
		}
		for (int i = 1; i < task_num; i++) memcpy(mat[i], mat[i * comm_sz], n * sizeof(float));//������߳�Ҫ����������Ƶ���ǰ��
		std::free(buff);
	}
	else {// ��0�Ž��̸�������Ľ��չ���
		MPI_Recv(mat[0], task_num * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// ����һ������
	for (int k = 0; k < n; k++) {
		// ������������Ǳ����̸�������񣬲�����������㲥
		if (k % comm_sz == my_rank) {
			int localIndexK = k / comm_sz;
			vt = _mm256_set1_ps(mat[localIndexK][k]);
			int j = k + 1;
			for (; j < n; j++)
			{
				//32�ֽڶ���
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

		// ������Ԫ����
		int end = my_rank > (k % comm_sz) ? k / comm_sz : k / comm_sz + 1;
		for (int i = task_num - 1; i >= end; i--) {
			vaik = _mm256_set1_ps(mat[i][k]);
			int j = k + 1;
			for (; j < n; ++j)
			{
				//32�ֽڶ���
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

	//����ɢ�ڸ������̵Ĵ𰸾ۺ�����
	if (my_rank == 0) {
		for (int i = task_num - 1; i > 0; i--)
			memcpy(mat[i * comm_sz], mat[i], n * sizeof(float));
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//ÿ��������ߴ���n / comm_sz + 1������
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Recv((void*)buff, count * n, MPI_FLOAT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//���ܶ�Ӧ�Ľ��̵�����
			for (int i = p; i < n; i += comm_sz) {
				memcpy(mat[i], buff + i / comm_sz * n, n * sizeof(float));
			}
		}
		std::free(buff);
	}
	else {
		MPI_Send(mat[0], task_num * n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);//�ȴ�ȫ���߳����
	end_time = MPI_Wtime();
	return end_time - start_time;
}

// MPIѭ������+OMP+AVX2�����㷨
double mpiAvx2OmpGauss(int threadNum = THREADNUM) {
	double start_time, end_time;
	MPI_Barrier(MPI_COMM_WORLD);//�����߳�ͬʱ��ʼ�����ڼ�ʱ
	start_time = MPI_Wtime();

	//�����������������÷���������ܳ�����ʣ�µ�����Ҳ���ȵػ��ָ���ǰ��������
	int task_num = my_rank < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;

	//���Ƚ�ÿ������Ҫ��������ݷŵ�������ǰ��
	if (my_rank == 0) {// 0�Ž��̸�������ĳ�ʼ�ַ�����
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//ÿ��������ߴ���n / comm_sz + 1������
		for (int p = 1; p < comm_sz; p++) {
			for (int i = p; i < n; i += comm_sz) {
				memcpy(buff + i / comm_sz * n, mat[i], n * sizeof(float));
			}
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Send(buff, count * n, MPI_FLOAT, p, 0, MPI_COMM_WORLD);//�����ݷ��͸���Ӧ�Ľ���
		}
		for (int i = 1; i < task_num; i++) memcpy(mat[i], mat[i * comm_sz], n * sizeof(float));//������߳�Ҫ����������Ƶ���ǰ��
		std::free(buff);
	}
	else {// ��0�Ž��̸�������Ľ��չ���
		MPI_Recv(mat[0], task_num * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// ����һ������
	int i, j, k;
	float pivot;
	__m256 vt, va, vaik, vakj, vaij;
#pragma omp parallel num_threads(threadNum) default(none) private(i, j, k, pivot, vt, va, vaik, vakj, vaij, pivotRow, my_rank, task_num) shared(mat, n, comm_sz)
	for (k = 0; k < n; k++) {
		// ������������Ǳ����̸�������񣬲�����������㲥
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

		//��Ԫ
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

	//����ɢ�ڸ������̵Ĵ𰸾ۺ�����
	if (my_rank == 0) {
		for (int i = task_num - 1; i > 0; i--)
		{
			memcpy(mat[i * comm_sz], mat[i], n * sizeof(float));
		}
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//ÿ��������ߴ���n / comm_sz + 1������
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Recv((void*)buff, count * n, MPI_FLOAT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//���ܶ�Ӧ�Ľ��̵�����
			for (int i = p; i < n; i += comm_sz) {
				memcpy(mat[i], buff + i / comm_sz * n, n * sizeof(float));
			}
		}
		std::free(buff);
	}
	else {
		MPI_Send(mat[0], task_num * n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);//�ȴ�ȫ���߳����
	end_time = MPI_Wtime();
	return end_time - start_time;
}
#else
//MPIѭ������+Neon�����㷨
double mpiNeonGauss(int threadNum = THREADNUM)
{
	double start_time, end_time;
	float32x4_t inv_vt, va, vaik, vakj, vaij;
	MPI_Barrier(MPI_COMM_WORLD);//�����߳�ͬʱ��ʼ�����ڼ�ʱ
	start_time = MPI_Wtime();

	//�����������������÷���������ܳ�����ʣ�µ�����Ҳ���ȵػ��ָ���ǰ��������
	int task_num = my_rank < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;

	//���Ƚ�ÿ������Ҫ��������ݷŵ�������ǰ��
	if (my_rank == 0) {// 0�Ž��̸�������ĳ�ʼ�ַ�����
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//ÿ��������ߴ���n / comm_sz + 1������
		for (int p = 1; p < comm_sz; p++) {
			for (int i = p; i < n; i += comm_sz) {
				memcpy(buff + i / comm_sz * n, mat[i], n * sizeof(float));
			}
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Send(buff, count * n, MPI_FLOAT, p, 0, MPI_COMM_WORLD);//�����ݷ��͸���Ӧ�Ľ���
		}
		for (int i = 1; i < task_num; i++) memcpy(mat[i], mat[i * comm_sz], n * sizeof(float));//������߳�Ҫ����������Ƶ���ǰ��
		std::free(buff);
	}
	else {// ��0�Ž��̸�������Ľ��չ���
		MPI_Recv(mat[0], task_num * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// ����һ������
	for (int k = 0; k < n; k++) {
		// ������������Ǳ����̸�������񣬲�����������㲥
		if (k % comm_sz == my_rank) {
			int localIndexK = k / comm_sz;
			inv_vt = vdupq_n_f32(1.0f / mat[localIndexK][k]);
			int j = k + 1;
			for (; j < n; j++)
			{
				//16�ֽڶ���
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

		// ������Ԫ����
		int end = my_rank > (k % comm_sz) ? k / comm_sz : k / comm_sz + 1;
		for (int i = task_num - 1; i >= end; i--) {
			vaik = vdupq_n_f32(mat[i][k]);
			int j = k + 1;
			for (; j < n; ++j)
			{
				//16�ֽڶ���
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

	//����ɢ�ڸ������̵Ĵ𰸾ۺ�����
	if (my_rank == 0) {
		for (int i = task_num - 1; i > 0; i--)
			memcpy(mat[i * comm_sz], mat[i], n * sizeof(float));
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//ÿ��������ߴ���n / comm_sz + 1������
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Recv((void*)buff, count * n, MPI_FLOAT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//���ܶ�Ӧ�Ľ��̵�����
			for (int i = p; i < n; i += comm_sz) {
				memcpy(mat[i], buff + i / comm_sz * n, n * sizeof(float));
			}
		}
		std::free(buff);
	}
	else {
		MPI_Send(mat[0], task_num * n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);//�ȴ�ȫ���߳����
	end_time = MPI_Wtime();
	return end_time - start_time;
}

// MPIѭ������+OMP+Neon�����㷨
double mpiNeonOmpGauss(int threadNum= THREADNUM) {
	double start_time, end_time;
	MPI_Barrier(MPI_COMM_WORLD);//�����߳�ͬʱ��ʼ�����ڼ�ʱ
	start_time = MPI_Wtime();

	//�����������������÷���������ܳ�����ʣ�µ�����Ҳ���ȵػ��ָ���ǰ��������
	int task_num = my_rank < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;

	//���Ƚ�ÿ������Ҫ��������ݷŵ�������ǰ��
	if (my_rank == 0) {// 0�Ž��̸�������ĳ�ʼ�ַ�����
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//ÿ��������ߴ���n / comm_sz + 1������
		for (int p = 1; p < comm_sz; p++) {
			for (int i = p; i < n; i += comm_sz) {
				memcpy(buff + i / comm_sz * n, mat[i], n * sizeof(float));
			}
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Send(buff, count * n, MPI_FLOAT, p, 0, MPI_COMM_WORLD);//�����ݷ��͸���Ӧ�Ľ���
		}
		for (int i = 1; i < task_num; i++) memcpy(mat[i], mat[i * comm_sz], n * sizeof(float));//������߳�Ҫ����������Ƶ���ǰ��
		std::free(buff);
	}
	else {// ��0�Ž��̸�������Ľ��չ���
		MPI_Recv(mat[0], task_num * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// ����һ������
	int i, j, k;
	float pivot;
	float32x4_t inv_vt, va, vaik, vakj, vaij;
#pragma omp parallel num_threads(threadNum) default(none) private(i, j, k, pivot, inv_vt, va, vaik, vakj, vaij,pivotRow, my_rank, task_num) shared(mat, n, comm_sz)
	for (k = 0; k < n; k++) {
		// ������������Ǳ����̸�������񣬲�����������㲥
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

	//����ɢ�ڸ������̵Ĵ𰸾ۺ�����
	if (my_rank == 0) {
		for (int i = task_num - 1; i > 0; i--)
		{
			memcpy(mat[i * comm_sz], mat[i], n * sizeof(float));
		}
		float* buff = (float*)malloc((n / comm_sz + 1) * n * sizeof(float));//ÿ��������ߴ���n / comm_sz + 1������
		for (int p = 1; p < comm_sz; p++) {
			int count = p < (n% comm_sz) ? n / comm_sz + 1 : n / comm_sz;
			MPI_Recv((void*)buff, count * n, MPI_FLOAT, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//���ܶ�Ӧ�Ľ��̵�����
			for (int i = p; i < n; i += comm_sz) {
				memcpy(mat[i], buff + i / comm_sz * n, n * sizeof(float));
			}
		}
		std::free(buff);
	}
	else {
		MPI_Send(mat[0], task_num * n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);//�ȴ�ȫ���߳����
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