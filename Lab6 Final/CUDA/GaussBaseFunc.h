#pragma once
#include "config.h"
//打印结果，用于调试
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

//检查结果
void check(int n, float* mat, float* ans)
{
#ifndef CHECK
	return;
#endif

	//n太大时打印并不好观察，故不再打印结果
	if (n <= MAX_PRINT_N)
	{
		printMatrix(n, "ans", ans);
		printMatrix(n, "matrix", mat);
	}

	//只要mat矩阵的秩与ans相同，则代表结果正确
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			if (abs(mat[i * n + j] - ans[i * n + j]) > eps)
			{
				printf("\nsize:%d\tcheck result: Wrong!\n", n);
				printf("matrix[%d][%d]=%f\n", i, j, mat[i * n + j]);
				printf("ans[%d][%d]=%f\n", i, j, ans[i * n + j]);
				exit(0);
			}

	//此时说明结果正确
	cout << "\nsize:" << n << "\tcheck result : Right!\n";
}

//为避免出现计算结果为naf或无穷等极端情况，采用先生成上三角矩阵再进行行变换的方法生成测试样例
void dataGenerate(int n, float* mat, float* ans)
{
	memset(mat, 0, n * n * sizeof(float));
	for (int i = 0; i < n; ++i)
	{
		mat[i * n + i] = 1.0f;
		for (int j = i + 1; j < n; ++j)
			mat[i * n + j] = rand() & 0xF;//使得生成数据在0到15之间,防止下面进行行变换生成数据时发生溢出
	}

	//当前的结果即为最终的答案
	memcpy(ans, mat, n * n * sizeof(float));

	int inv_rate = n >> 3;
	bool reGenerate = false;
	while (true)
	{
		//如果重新生成数据，则将mat复制为ans
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

//串行高斯消去
double serial(int n, float* mat, int blockSize = 0)
{
	auto start = chrono::high_resolution_clock::now();

	for (int k = 0; k < n; ++k)
	{
		//对行进行归一化
		float pivot = mat[k * n + k];
		for (int j = k + 1; j < n; ++j)mat[k * n + j] /= pivot;
		mat[k * n + k] = 1.0f;

		//对右下角的矩阵进行消元
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

//计时函数
void timing(double(*func)(int, float*, int), int n = 2048, int blockSize = 1024)
{
	double totalTime = 0.0f;
	float* matrix = (float*)malloc(n * n * sizeof(float));
	float* ans = (float*)malloc(n * n * sizeof(float));

	//先进行WARMUP次热身，不计入耗时，使得测量时间更加精确
	for (int run = 0; run < REPEAT_NUM + WARMUP; run++) {
		dataGenerate(n, matrix, ans);
		float duration = func(n, matrix, blockSize);
		if (run >= WARMUP) totalTime += duration;
	}
	cout << totalTime / REPEAT_NUM << '\t';

	//检验答案
	check(n, matrix, ans);

	//释放内存
	free(matrix);
	free(ans);
}

//计时函数
void timingAllMatSize(string nameArr[], double(*funcArr[])(int, float*, int), int funcNum, int maxN = 2048, int blockSize = 1024)
{
	//输出表头
	cout << "N\t";
	for (int i = 0; i < funcNum; ++i)
		cout << nameArr[i] << '\t';
	cout << endl;

	//进行测量
	for (int n = 64; n <= maxN; n += 64) {
		cout << n << '\t';
		for (int funcIndex = 0; funcIndex < funcNum; ++funcIndex) {
			timing(funcArr[funcIndex], n, blockSize);
		}
		cout << endl;
	}
}

//计时函数
void timingAllBlockSize(string nameArr[], double(*funcArr[])(int, float*, int), int funcNum, int n = 2048, int maxBlockSize = 1024)
{
	//输出表头
	cout << "blockSize\t";
	for (int i = 0; i < funcNum; ++i)
		cout << nameArr[i] << '\t';
	cout << endl;

	//进行测量
	for (int blockSize = 256; blockSize <= maxBlockSize; blockSize += 256) {
		cout << blockSize << '\t';
		for (int funcIndex = 0; funcIndex < funcNum; ++funcIndex) {
			timing(funcArr[funcIndex], n, blockSize);
		}
		cout << endl;
	}
}