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
#ifdef CHECK
	if (n <= MAX_PRINT_N)printMatrix(n, "original:", mat);
#endif
}

//���и�˹��ȥ
void serial(int n, float* mat, int threadNum = 1)
{
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
}

//���Ժ���
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

				//��������岿��
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

//���Ժ���
void timingAllThreadNum(string name[], void(*funArr[])(int, float*, int), int funNum)
{
	cout << "��ͬ�߳����µ����ܱ�\n\t";
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

				//��������岿��
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

//���Ժ���
void timing(string name, void(*fun)(int, float*, int), int n=2048, int threadNum=8)
{
	std::cout << name + "_" + to_string(threadNum) << '\t';
	double totalTime = 0.0f;
	float* matrix = (float*)malloc(n * n * sizeof(float));
	float* ans = (float*)malloc(n * n * sizeof(float));

	dataGenerate(n, matrix, ans);
	struct timeval start, end;
	gettimeofday(&start, NULL);

	//��������岿��
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
//�����AVX2��˹��Ԫ
void AVX2Aligned(int n, float* mat, int threadNum = 1)
{
	__m256 vt, va, vaik, vakj, vaij;
	for (int k = 0; k < n; ++k)
	{
		//��һ��
		vt = _mm256_set1_ps(mat[k * n + k]);
		int j = k + 1;
		for (; j < n; j++)
		{
			//32�ֽڶ���
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

		//��Ԫ
		for (int i = k + 1; i < n; ++i)
		{
			vaik = _mm256_set1_ps(mat[i * n + k]);
			int j = k + 1;
			for (; j < n; ++j)
			{
				//32�ֽڶ���
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
//����ġ��ںϳ˼��ġ�NEON��˹��Ԫ
void neonAlignedMLS(int n, float* mat, int threadNum = 1)
{
	float32x4_t vt, va, vaik, vakj, vaij;
	for (int k = 0; k < n; ++k)
	{
		//��һ��
		vt = vdupq_n_f32(mat[k * n + k]);
		int j = k + 1;
		for (; j < n; j++)
		{
			//16�ֽڶ���
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

		//��Ԫ
		for (int i = k + 1; i < n; ++i)
		{
			vaik = vdupq_n_f32(mat[i * n + k]);
			int j = k + 1;
			for (; j < n; ++j)
			{
				//16�ֽڶ���
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