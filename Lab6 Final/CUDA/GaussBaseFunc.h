#pragma once
#include "config.h"
//��ӡ��������ڵ���
void printMatrix(int n, const char* name, float* mat)
{
	printf("%s:\n", name);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			printf("%.1f\t", mat[i * n + j]);
		printf("\n");
	}
}

//�����
void check(int n, float* mat, float* ans)
{
#ifndef CHECK
	return;
#endif

	//n̫��ʱ��ӡ�����ù۲죬�ʲ��ٴ�ӡ���
	if (n <= MAX_PRINT_N)
	{
		printMatrix(n, "ans", ans);
		printMatrix(n, "matrix", mat);
	}

	//ֻҪmat���������ans��ͬ�����������ȷ
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			if (abs(mat[i * n + j] - ans[i * n + j]) > eps)
			{
				printf("\nsize:%d\tcheck result: Wrong!\n", n);
				printf("matrix[%d][%d]=%f\n", i, j, mat[i * n + j]);
				printf("ans[%d][%d]=%f\n", i, j, ans[i * n + j]);
				exit(0);
			}

	//��ʱ˵�������ȷ
	cout << "\nsize:" << n << "\tcheck result : Right!\n";
}

//Ϊ������ּ�����Ϊnaf������ȼ�����������������������Ǿ����ٽ����б任�ķ������ɲ�������
void dataGenerate(int n, float* mat, float* ans)
{
	memset(mat, 0, n * n * sizeof(float));
	for (int i = 0; i < n; ++i)
	{
		mat[i * n + i] = 1.0f;
		for (int j = i + 1; j < n; ++j)
			mat[i * n + j] = rand() & 0xF;//ʹ������������0��15֮��,��ֹ��������б任��������ʱ�������
	}

	//��ǰ�Ľ����Ϊ���յĴ�
	memcpy(ans, mat, n * n * sizeof(float));

	int inv_rate = n >> 3;
	bool reGenerate = false;
	while (true)
	{
		//��������������ݣ���mat����Ϊans
	REGEn:
		if (reGenerate) memcpy(mat, ans, n * n * sizeof(float));
		for (int k = 0; k < n; ++k)
			for (int i = k; i < n; ++i)
			{
				if (rand() % inv_rate) continue;
				for (int j = 0; j < n; ++j)
				{
					mat[i * n + j] += mat[k * n + j];
					if (mat[i * n + j] > floatLimit)
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

//���и�˹��ȥ
double serial(int n, float* mat, int blockSize = 0)
{
	auto start = chrono::high_resolution_clock::now();

	for (int k = 0; k < n; ++k)
	{
		//���н��й�һ��
		float pivot = mat[k * n + k];
		for (int j = k + 1; j < n; ++j)mat[k * n + j] /= pivot;
		mat[k * n + k] = 1.0f;

		//�����½ǵľ��������Ԫ
		for (int i = k + 1; i < n; ++i)
		{
			float pivot = mat[i * n + k];
			for (int j = k + 1; j < n; ++j)
				mat[i * n + j] -= mat[k * n + j] * pivot;
			mat[i * n + k] = 0.0f;
		}
	}

	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::microseconds>(stop - start) / 1000000.0;
	return duration.count();
}

//��ʱ����
void timing(double(*func)(int, float*, int), int n = 2048, int blockSize = 1024)
{
	double totalTime = 0.0f;
	float* matrix = (float*)malloc(n * n * sizeof(float));
	float* ans = (float*)malloc(n * n * sizeof(float));

	//�Ƚ���WARMUP�������������ʱ��ʹ�ò���ʱ����Ӿ�ȷ
	for (int run = 0; run < REPEAT_NUM + WARMUP; run++) {
		dataGenerate(n, matrix, ans);
		float duration = func(n, matrix, blockSize);
		if (run >= WARMUP) totalTime += duration;
	}
	cout << totalTime / REPEAT_NUM << '\t';

	//�����
	check(n, matrix, ans);

	//�ͷ��ڴ�
	free(matrix);
	free(ans);
}

//��ʱ����
void timingAllMatSize(string nameArr[], double(*funcArr[])(int, float*, int), int funcNum, int maxN = 2048, int blockSize = 1024)
{
	//�����ͷ
	cout << "N\t";
	for (int i = 0; i < funcNum; ++i)
		cout << nameArr[i] << '\t';
	cout << endl;

	//���в���
	for (int n = 64; n <= maxN; n += 64) {
		cout << n << '\t';
		for (int funcIndex = 0; funcIndex < funcNum; ++funcIndex) {
			timing(funcArr[funcIndex], n, blockSize);
		}
		cout << endl;
	}
}

//��ʱ����
void timingAllBlockSize(string nameArr[], double(*funcArr[])(int, float*, int), int funcNum, int n = 2048, int maxBlockSize = 1024)
{
	//�����ͷ
	cout << "blockSize\t";
	for (int i = 0; i < funcNum; ++i)
		cout << nameArr[i] << '\t';
	cout << endl;

	//���в���
	for (int blockSize = 256; blockSize <= maxBlockSize; blockSize += 256) {
		cout << blockSize << '\t';
		for (int funcIndex = 0; funcIndex < funcNum; ++funcIndex) {
			timing(funcArr[funcIndex], n, blockSize);
		}
		cout << endl;
	}
}