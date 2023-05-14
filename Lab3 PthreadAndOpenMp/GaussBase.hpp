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
#ifdef CHECK
	if (n <= MAX_PRINT_N)printMatrix(n, "original:", mat);
#endif
}

//串行高斯消去
void serial(int n, float* mat, int threadNum = 1)
{
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
}

//测试函数
void timingAll(string name[], void(*funArr[])(int, float*, int), int funNum, int threadNum=8)
{
	for(int funIndex=0;funIndex<funNum;++funIndex)
	{
		cout << name[funIndex]<< '\t';
		for (int n = 64; n <= MAXN; n += 64)
		{
			double totalTime = 0.0f;
			for (int step = 0; step < REPEAT_NUM; ++step)
			{
				float* matrix = (float*)malloc(n * n * sizeof(float));
				float* ans = (float*)malloc(n * n * sizeof(float));

				dataGenerate(n, matrix, ans);
				struct timeval start, end;
				gettimeofday(&start, NULL);

				//程序的主体部分
				funArr[funIndex](n, matrix, threadNum);

				gettimeofday(&end, NULL);

	#ifdef CHECK
				check(n, matrix, ans);
	#endif

				totalTime += (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
				free(matrix);
				free(ans);
			}
			cout << totalTime / (double)REPEAT_NUM << "\t";
	#ifdef CHECK
			cout << "\nsize:" << n << "\tcheck result : Right!\n";
	#endif
		}
		cout << "\n";
	}
}

//测试函数
void timingAllThreadNum(string name[], void(*funArr[])(int, float*, int), int funNum)
{
	cout << "不同线程数下的性能表\n\t";
	for (int threadIndex = 0; threadIndex < 19; ++threadIndex)
	{
		int threadNum = threadNumArr[threadIndex];
		cout << "threadNum=" + to_string(threadNum) << '\t';
	}
	cout << "\n";

	for (int funIndex = 0; funIndex < funNum; ++funIndex)
	{
		cout << name[funIndex] << '\t';
		for (int threadIndex = 0; threadIndex < 19; ++threadIndex)
		{
			int threadNum = threadNumArr[threadIndex];
			double totalTime = 0.0f;
			int n = MAXN;
			for (int step = 0; step < REPEAT_NUM; ++step)
			{
				float* matrix = (float*)malloc(n * n * sizeof(float));
				float* ans = (float*)malloc(n * n * sizeof(float));

				dataGenerate(n, matrix, ans);
				struct timeval start, end;
				gettimeofday(&start, NULL);

				//程序的主体部分
				funArr[funIndex](n, matrix, threadNum);

				gettimeofday(&end, NULL);

#ifdef CHECK
				check(n, matrix, ans);
#endif

				totalTime += (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
				free(matrix);
				free(ans);
			}
			cout << totalTime / (double)REPEAT_NUM << "\t";
#ifdef CHECK
			cout << "\nsize:" << n << "\tcheck result : Right!\n";
#endif
		}
		cout << '\n';
	}
}

//测试函数
void timing(string name, void(*fun)(int, float*, int), int n=2048, int threadNum=8)
{
	std::cout << name + "_" + to_string(threadNum) << '\t';
	double totalTime = 0.0f;
	float* matrix = (float*)malloc(n * n * sizeof(float));
	float* ans = (float*)malloc(n * n * sizeof(float));

	dataGenerate(n, matrix, ans);
	struct timeval start, end;
	gettimeofday(&start, NULL);

	//程序的主体部分
	fun(n, matrix, threadNum);

	gettimeofday(&end, NULL);

#ifdef CHECK
	check(n, matrix, ans);
#endif

	totalTime += (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
	free(matrix);
	free(ans);
#ifdef CHECK
	cout << "\nsize:"<<n<<"\tcheck result : Right!\n";
#endif
	cout << totalTime << endl;
}

#ifndef __ARM_NEON
//对齐的AVX2高斯消元
void AVX2Aligned(int n, float* mat, int threadNum = 1)
{
	__m256 vt, va, vaik, vakj, vaij;
	for (int k = 0; k < n; ++k)
	{
		//归一化
		vt = _mm256_set1_ps(mat[k * n + k]);
		int j = k + 1;
		for (; j < n; j++)
		{
			//32字节对齐
			if (!((ull(&mat[k * n + j])) & 0x1f)) break;
			mat[k * n + j] /= mat[k * n + k];
		}
		for (; j + 8 <= n; j += 8)
		{
			va = _mm256_load_ps(&mat[k * n + j]);
			va = _mm256_div_ps(va, vt);
			_mm256_store_ps(&mat[k * n + j], va);
		}
		for (; j < n; j++) mat[k * n + j] /= mat[k * n + k];
		mat[k * n + k] = 1.0f;

		//消元
		for (int i = k + 1; i < n; ++i)
		{
			vaik = _mm256_set1_ps(mat[i * n + k]);
			int j = k + 1;
			for (; j < n; ++j)
			{
				//32字节对齐
				if (!((ull(&mat[i * n + j])) & 0x1f)) break;
				mat[i * n + j] -= mat[k * n + j] * mat[i * n + k];
			}
			for (; j + 8 <= n; j += 8)
			{
				vakj = _mm256_load_ps(&mat[k * n + j]);
				vaij = _mm256_load_ps(&mat[i * n + j]);
				vaij = _mm256_fnmadd_ps(vakj, vaik, vaij);
				_mm256_store_ps(&mat[i * n + j], vaij);
			}
			for (; j < n; j++) mat[i * n + j] -= mat[k * n + j] * mat[i * n + k];
			mat[i * n + k] = 0.0f;
		}
	}
}

#else
//对齐的、融合乘减的、NEON高斯消元
void neonAlignedMLS(int n, float* mat, int threadNum = 1)
{
	float32x4_t vt, va, vaik, vakj, vaij;
	for (int k = 0; k < n; ++k)
	{
		//归一化
		vt = vdupq_n_f32(mat[k * n + k]);
		int j = k + 1;
		for (; j < n; j++)
		{
			//16字节对齐
			if (!((ull(&mat[k * n + j])) & 0xf)) break;
			mat[k * n + j] /= mat[k * n + k];
		}
		for (; j + 4 <= n; j += 4)
		{
			va = vld1q_f32(&mat[k * n + j]);
			va = vmulq_f32(va, vt);
			vst1q_f32(&mat[k * n + j], va);
		}
		for (; j < n; j++) mat[k * n + j] /= mat[k * n + k];
		mat[k * n + k] = 1.0f;

		//消元
		for (int i = k + 1; i < n; ++i)
		{
			vaik = vdupq_n_f32(mat[i * n + k]);
			int j = k + 1;
			for (; j < n; ++j)
			{
				//16字节对齐
				if (!((ull(&mat[i * n + j])) & 0xf)) break;
				mat[i * n + j] -= mat[k * n + j] * mat[i * n + k];
			}
			for (; j + 4 <= n; j += 4)
			{
				vakj = vld1q_f32(&mat[k * n + j]);
				vaij = vld1q_f32(&mat[i * n + j]);
				vaij = vmlsq_f32(vakj, vaik, vaij);
				vst1q_f32(&mat[i * n + j], vaij);
			}
			for (; j < n; j++) mat[i * n + j] -= mat[k * n + j] * mat[i * n + k];
			mat[i * n + k] = 0.0f;
		}
	}
}
#endif