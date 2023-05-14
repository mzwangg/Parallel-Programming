#include "config.h"
#include <omp.h>

// 循环划分+静态划分
void openmpCyclicStatic(int n, float* mat, int threadNum)
{
	int i, j, k;
	float pivot;
#pragma omp parallel num_threads(threadNum) default(none) private(i, j, k, pivot) shared(mat, n)
	for (k = 0; k < n; k++)
	{
		pivot = mat[k * n + k];
#pragma omp for
		for (j = k + 1; j < n; j++)
		{
			mat[k * n + j] /= pivot;
		}
#pragma omp single
		mat[k * n + k] = 1.0f;
#pragma omp for schedule(static, 1)
		for (i = k + 1; i < n; i++)
		{
			pivot = mat[i * n + k];
#pragma omp simd
			for (j = k + 1; j < n; j++)
			{
				mat[i * n + j] -= pivot * mat[k * n + j];
			}
			mat[i * n + k] = 0;
		}
	}
}

// 块循环划分+静态划分
void openmpBlockCyclicStatic(int n, float* mat, int threadNum)
{
	int i, j, k;
	float pivot;
#pragma omp parallel num_threads(threadNum) default(none) private(i, j, k, pivot) shared(mat, n)
	for (k = 0; k < n; k++)
	{
		pivot = mat[k * n + k];
#pragma omp for
		for (j = k + 1; j < n; j++)
		{
			mat[k * n + j] /= pivot;
		}
#pragma omp single
		mat[k * n + k] = 1.0f;
#pragma omp for schedule(static, 4)
		for (i = k + 1; i < n; i++)
		{
			pivot = mat[i * n + k];
#pragma omp simd
			for (j = k + 1; j < n; j++)
			{
				mat[i * n + j] -= pivot * mat[k * n + j];
			}
			mat[i * n + k] = 0;
		}
	}
}

// 块划分+静态划分
void openmpBlockStatic(int n, float* mat, int threadNum)
{
	int i, j, k;
	float pivot;
	int chunk;
#pragma omp parallel num_threads(threadNum) default(none) private(i, j, k, pivot) shared(mat, n, chunk, threadNum)
	for (k = 0; k < n; k++)
	{
		pivot = mat[k * n + k];
#pragma omp for
		for (j = k + 1; j < n; j++)
		{
			mat[k * n + j] /= pivot;
		}
#pragma omp single
		{
			mat[k * n + k] = 1.0f;
			chunk = (int)ceil((n - k - 1) / double(threadNum));
		}
#pragma omp for schedule(static, chunk)
		for (i = k + 1; i < n; i++)
		{
			pivot = mat[i * n + k];
#pragma omp simd
			for (j = k + 1; j < n; j++)
			{
				mat[i * n + j] -= pivot * mat[k * n + j];
			}
			mat[i * n + k] = 0;
		}
	}
}

// 动态划分
void openmpDynamic(int n, float* mat, int threadNum)
{
	int i, j, k;
	float pivot;
#pragma omp parallel num_threads(threadNum) default(none) private(i, j, k, pivot) shared(mat, n)
	for (k = 0; k < n; k++)
	{
		pivot = mat[k * n + k];
#pragma omp for
		for (j = k + 1; j < n; j++)
		{
			mat[k * n + j] /= pivot;
		}
#pragma omp single
		mat[k * n + k] = 1.0f;
#pragma omp for schedule(dynamic, 1)
		for (i = k + 1; i < n; i++)
		{
			pivot = mat[i * n + k];
#pragma omp simd
			for (j = k + 1; j < n; j++)
			{
				mat[i * n + j] -= pivot * mat[k * n + j];
			}
			mat[i * n + k] = 0;
		}
	}
}

// 块循环划分+动态划分
void openmpBlockCyclicDynamic(int n, float* mat, int threadNum)
{
	int i, j, k;
	float pivot;
#pragma omp parallel num_threads(threadNum) default(none) private(i, j, k, pivot) shared(mat, n)
	for (k = 0; k < n; k++)
	{

		pivot = mat[k * n + k];
#pragma omp for
		for (j = k + 1; j < n; j++)
		{
			mat[k * n + j] /= pivot;
		}
#pragma omp single
		mat[k * n + k] = 1.0f;
#pragma omp for schedule(dynamic, 4)
		for (i = k + 1; i < n; i++)
		{
			pivot = mat[i * n + k];
#pragma omp simd
			for (j = k + 1; j < n; j++)
			{
				mat[i * n + j] -= pivot * mat[k * n + j];
			}
			mat[i * n + k] = 0;
		}
	}
}

// 列划分
void openmpColumn(int n, float* mat, int threadNum)
{
	int i, j, k;
#pragma omp parallel num_threads(threadNum), default(none), private(i, j, k), shared(mat, n)
	for (k = 0; k < n; k++)
	{
#pragma omp for schedule(dynamic)
		for (j = k + 1; j < n; j++)
		{
			mat[k * n + j] /= mat[k * n + k];
			for (i = k + 1; i < n; i++)
			{
				mat[i * n + j] -= mat[i * n + k] * mat[k * n + j];
			}
		}
#pragma omp single
		{
			mat[k * n + k] = 1.0f;
			for (i = k + 1; i < n; i++)
			{
				mat[i * n + k] = 0.0f;
			}
		}
	}
}

// guided划分
void openmpGuided(int n, float* mat, int threadNum)
{
	int i, j, k;
	float pivot;
#pragma omp parallel num_threads(threadNum) default(none) private(i, j, k, pivot) shared(mat, n)
	for (k = 0; k < n; k++)
	{
		pivot = mat[k * n + k];
#pragma omp for
		for (j = k + 1; j < n; j++)
		{
			mat[k * n + j] /= pivot;
		}
#pragma omp single
		mat[k * n + k] = 1.0f;
#pragma omp for schedule(guided)
		for (i = k + 1; i < n; i++)
		{
			pivot = mat[i * n + k];
#pragma omp simd
			for (j = k + 1; j < n; j++)
			{
				mat[i * n + j] -= pivot * mat[k * n + j];
			}
			mat[i * n + k] = 0;
		}
	}
}

// 动态线程
void openmpdynamicThread(int n, float* mat, int threadNum)
{
	int i, j, k;
	float pivot;
	for (k = 0; k < n; k++)
	{
#pragma omp parallel num_threads(threadNum) default(none) private(i, j, pivot) shared(mat, k, n)
		{
		pivot = mat[k * n + k];
#pragma omp for
		for (j = k + 1; j < n; j++)
		{
			mat[k * n + j] /= pivot;
		}
#pragma omp single
		mat[k * n + k] = 1.0f;
#pragma omp for schedule(dynamic, 1)
			for (i = k + 1; i < n; i++)
			{
				pivot = mat[i * n + k];
#pragma omp simd
				for (j = k + 1; j < n; j++)
				{
					mat[i * n + j] -= pivot * mat[k * n + j];
				}
				mat[i * n + k] = 0.0f;
			}
		}
	}
}

#ifndef __ARM_NEON
// openmp + manual simd
void openmpManualSIMD(int n, float* mat, int threadNum)
{
	int i, j, k;
	float pivot;
	__m256 vt, va, vaik, vakj, vaij;
#pragma omp parallel num_threads(threadNum) default(none) private(i, j, k, pivot, vt, va, vaik, vakj, vaij) shared(mat, n)
	for (k = 0; k < n; k++)
	{
		vt = _mm256_set1_ps(mat[k * n + k]);
		pivot = mat[k * n + k];

#pragma omp for schedule(dynamic)
		for (int j = k + 1; j <= n - 8; j += 8)
		{
			va = _mm256_loadu_ps(&mat[k * n + j]);
			va = _mm256_div_ps(va, vt);
			_mm256_storeu_ps(&mat[k * n + j], va);
		}

#pragma omp single
		{
			for (int j = n - ((n - k - 1) % 8); j < n; j++)
			{
				mat[k * n + j] /= pivot;
			}
			mat[k * n + k] = 1;
		}

#pragma omp for schedule(dynamic)
		for (int i = k + 1; i < n; i++)
		{
			vaik = _mm256_set1_ps(mat[i * n + k]);
			int j = k + 1;
			for (; j < n; ++j)
			{
				//32锟街节讹拷锟斤拷
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
#endif

// openmp offloading
void calculate_openmp_offloading(int n, float* mat, int threadNum)
{
	int is_cpu = false;
#pragma omp target map(tofrom: mat[0:n*n]) map(from: is_cpu) map(to: n)
	{
		int i, j, k;
		is_cpu = omp_is_initial_device();
		for (k = 0; k < n; k++)
		 {
#pragma omp single
			{
				for (j = k + 1; j < n; j++) {
					mat[k * n + j] /= mat[k * n + k];
				}
				mat[k * n + k] = 1.0f;
			}
#pragma omp for
			for (i = k + 1; i < n; i++) 
			{
				for (j = k + 1; j < n; j++) 
				{
					mat[i * n + j] -= mat[i * n + k] * mat[k * n + j];
				}
				mat[i * n + k] = 0.0f;
			}
		}
	}
#ifdef CHECK
	cout << (is_cpu ? "CPU" : "GPU") << endl;
#endif
}