#include "config.h"
#include <pthread.h>
#include <semaphore.h>

pthread_t threads[MAXN];
gaussPara thread_param_t[MAXN];
pthread_barrier_t barrier_Divsion;
pthread_barrier_t barrier_Elimination;
sem_t sem_leader;
sem_t sem_Divsion[MAXN];
sem_t sem_Elimination[MAXN];
pthread_mutex_t mutexI;
int myIndex = 0;

#ifndef __ARM_NEON
// �黮���̺߳���
void* threadFunc_block(void* param)
{
	__m256 vt, va, vaik, vakj, vaij;
	gaussPara* thread_param = (gaussPara*)param;
	int n = thread_param->n;
	int t_id = thread_param->t_id;
	int threadNum = thread_param->threadNum;
	float* mat = thread_param->mat;

	for (int k = 0; k < n; k++)
	{
		// ���������߳̽��г�������
		int step = threadNum * 8;
		vt = _mm256_set1_ps(mat[k * n + k]);
		int j = k + 1 + t_id * 8;
		for (; j + 8 <= n; j += step)
		{
			va = _mm256_loadu_ps(&mat[k * n + j]);
			va = _mm256_div_ps(va, vt);
			_mm256_storeu_ps(&mat[k * n + j], va);
		}
		if (t_id == threadNum - 1)//�߳�threadNum-1�������������µ����ݣ���Ϊ���߳����п����ڴ�ʱ����
		{
			for (j = n - ((n - k - 1) % 8); j < n; j++)
				mat[k * n + j] /= mat[k * n + k];
		}

		pthread_barrier_wait(&barrier_Divsion);

		//����ֵ��ֵΪ1����������if (t_id == threadNum-1)����У����?���������̻߳�û��ɳ���?
		//֮�󲢲����õ�������ݣ������ڴ˿̸�ֵ�?1�����������?
		if (t_id == 0) mat[k * n + k] = 1.0f;

		// �黮�����񣬻ᵼ�¸��ز���
		int L = ceil((n - k) * 1.0f / (threadNum - 1));
		for (int i = k + 1 + t_id * L; i < n && i < k + 1 + (t_id + 1) * L; i++)
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

		// �����߳�׼��������һ��
		pthread_barrier_wait(&barrier_Elimination);
	}
	if (t_id != 0)pthread_exit(NULL);
	return NULL;
}

// �黮��
void pthreadBlock(int n, float* mat, int threadNum)
{
	// ��ʼ��
	pthread_barrier_init(&barrier_Divsion, NULL, threadNum);
	pthread_barrier_init(&barrier_Elimination, NULL, threadNum);

	// �����߳�
	for (int i = threadNum - 1; i >= 0; i--)//���򴴽�����֤���һ���߳����п�������ִ�������
	{
		thread_param_t[i].t_id = i;
		thread_param_t[i].mat = mat;
		thread_param_t[i].n = n;
		thread_param_t[i].threadNum = threadNum;
		if (i != 0)pthread_create(&threads[i], NULL, threadFunc_block, (void*)(&thread_param_t[i]));
	}
	threadFunc_block((void*)(&thread_param_t[0]));//���߳���Ϊ�߳�0���м���

	// �ȴ������߳�ִ�����?
	for (int i = 1; i < threadNum; i++)
	{
		pthread_join(threads[i], NULL);
	}

	// ����
	pthread_barrier_destroy(&barrier_Divsion);
	pthread_barrier_destroy(&barrier_Elimination);
}

// ѭ�������̺߳���
void* threadFunc_cyclic(void* param)
{
	__m256 vt, va, vaik, vakj, vaij;
	gaussPara* thread_param = (gaussPara*)param;
	int n = thread_param->n;
	int t_id = thread_param->t_id;
	int threadNum = thread_param->threadNum;
	float* mat = thread_param->mat;

	for (int k = 0; k < n; k++)
	{
		// ���������߳̽��г�������
		int step = threadNum * 8;
		vt = _mm256_set1_ps(mat[k * n + k]);
		int j = k + 1 + t_id * 8;
		for (; j + 8 <= n; j += step)
		{
			va = _mm256_loadu_ps(&mat[k * n + j]);
			va = _mm256_div_ps(va, vt);
			_mm256_storeu_ps(&mat[k * n + j], va);
		}
		if (t_id == threadNum - 1)//�߳�threadNum-1�������������µ����ݣ���Ϊ���߳����п����ڴ�ʱ����
		{
			for (j = n - ((n - k - 1) % 8); j < n; j++)
				mat[k * n + j] /= mat[k * n + k];
		}

		pthread_barrier_wait(&barrier_Divsion);

		//����ֵ��ֵΪ1����������if (t_id == threadNum-1)����У����?���������̻߳�û��ɳ���?
		//֮�󲢲����õ�������ݣ������ڴ˿̸�ֵ�?1�����������?
		if (t_id == 0) mat[k * n + k] = 1.0f;

		// ѭ����������ʹ���ش��¾���
		for (int i = k + 1 + t_id; i < n; i += threadNum)
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

		// �����߳�׼��������һ��
		pthread_barrier_wait(&barrier_Elimination);
	}
	if (t_id != 0)pthread_exit(NULL);
	return NULL;
}

// ѭ������
void pthreadCyclic(int n, float* mat, int threadNum)
{
	// ��ʼ��
	pthread_barrier_init(&barrier_Divsion, NULL, threadNum);
	pthread_barrier_init(&barrier_Elimination, NULL, threadNum);

	// �����߳�
	for (int i = threadNum - 1; i >= 0; i--)//���򴴽�����֤���һ���߳����п�������ִ�������
	{
		thread_param_t[i].t_id = i;
		thread_param_t[i].mat = mat;
		thread_param_t[i].n = n;
		thread_param_t[i].threadNum = threadNum;
		if (i != 0)pthread_create(&threads[i], NULL, threadFunc_cyclic, (void*)(&thread_param_t[i]));
	}
	threadFunc_cyclic((void*)(&thread_param_t[0]));//���߳���Ϊ�߳�0���м���

	// �ȴ������߳�ִ�����?
	for (int i = 1; i < threadNum; i++)
	{
		pthread_join(threads[i], NULL);
	}

	// ����
	pthread_barrier_destroy(&barrier_Divsion);
	pthread_barrier_destroy(&barrier_Elimination);
}

// SSE+ѭ�������̺߳���
void* threadFunc_SSE(void* param)
{
	__m128 vt, va, vaik, vakj, vaij, vamul;
	gaussPara* thread_param = (gaussPara*)param;
	int n = thread_param->n;
	int t_id = thread_param->t_id;
	int threadNum = thread_param->threadNum;
	float* mat = thread_param->mat;

	for (int k = 0; k < n; k++)
	{
		// ���������߳̽��г�������
		int step = threadNum * 4;
		vt = _mm_set1_ps(mat[k * n + k]);
		int j = k + 1 + t_id * 4;
		for (; j + 4 <= n; j += step)
		{
			va = _mm_loadu_ps(&mat[k * n + j]);
			va = _mm_div_ps(va, vt);
			_mm_storeu_ps(&mat[k * n + j], va);
		}
		if (t_id == threadNum - 1)//�߳�threadNum-1�������������µ����ݣ���Ϊ���߳����п����ڴ�ʱ����
		{
			for (j = n - ((n - k - 1) % 4); j < n; j++)
				mat[k * n + j] /= mat[k * n + k];
		}

		pthread_barrier_wait(&barrier_Divsion);

		//����ֵ��ֵΪ1����������if (t_id == threadNum-1)����У����?���������̻߳�û��ɳ���?
		//֮�󲢲����õ�������ݣ������ڴ˿̸�ֵ�?1�����������?
		if (t_id == 0) mat[k * n + k] = 1.0f;

		// ѭ����������ʹ���ش��¾���
		for (int i = k + 1 + t_id; i < n; i += threadNum)
		{
			vaik = _mm_set1_ps(mat[i * n + k]);
			int j = k + 1;
			for (; j < n; ++j)
			{
				//16�ֽڶ���
				if (!((ull(&mat[i * n + j])) & 0xf)) break;
				mat[i * n + j] -= mat[k * n + j] * mat[i * n + k];
			}
			for (; j + 4 <= n; j += 4)
			{
				vakj = _mm_load_ps(&mat[k * n + j]);
				vaij = _mm_load_ps(&mat[i * n + j]);
				vamul = _mm_mul_ps(vakj, vaik);
				vaij = _mm_sub_ps(vaij, vamul);
				_mm_store_ps(&mat[i * n + j], vaij);
			}
			for (; j < n; j++) mat[i * n + j] -= mat[k * n + j] * mat[i * n + k];
			mat[i * n + k] = 0.0f;
		}

		// �����߳�׼��������һ��
		pthread_barrier_wait(&barrier_Elimination);
	}
	if (t_id != 0)pthread_exit(NULL);
	return NULL;
}

// SSE+ѭ������
void pthreadSSE(int n, float* mat, int threadNum)
{
	// ��ʼ��
	pthread_barrier_init(&barrier_Divsion, NULL, threadNum);
	pthread_barrier_init(&barrier_Elimination, NULL, threadNum);

	// �����߳�
	for (int i = threadNum - 1; i >= 0; i--)//���򴴽�����֤���һ���߳����п�������ִ�������
	{
		thread_param_t[i].t_id = i;
		thread_param_t[i].mat = mat;
		thread_param_t[i].n = n;
		thread_param_t[i].threadNum = threadNum;
		if (i != 0)pthread_create(&threads[i], NULL, threadFunc_SSE, (void*)(&thread_param_t[i]));
	}
	threadFunc_SSE((void*)(&thread_param_t[0]));//���߳���Ϊ�߳�0���м���

	// �ȴ������߳�ִ�����?
	for (int i = 1; i < threadNum; i++)
	{
		pthread_join(threads[i], NULL);
	}

	// ����
	pthread_barrier_destroy(&barrier_Divsion);
	pthread_barrier_destroy(&barrier_Elimination);
}

// ѭ������+�ź����̺߳���
void* threadFunc_cyclicSem(void* param)
{
	__m256 vt, va, vaik, vakj, vaij;
	gaussPara* thread_param = (gaussPara*)param;
	int n = thread_param->n;
	int t_id = thread_param->t_id;
	int threadNum = thread_param->threadNum;
	float* mat = thread_param->mat;

	for (int k = 0; k < n; k++)
	{
		// ���������߳̽��г�������
		int step = threadNum * 8;
		vt = _mm256_set1_ps(mat[k * n + k]);
		int j = k + 1 + t_id * 8;
		for (; j + 8 <= n; j += step)
		{
			va = _mm256_loadu_ps(&mat[k * n + j]);
			va = _mm256_div_ps(va, vt);
			_mm256_storeu_ps(&mat[k * n + j], va);
		}
		//�߳�threadNum-1�������������µ����ݣ���Ϊ���߳����п����ڴ�ʱ����
		if (t_id == threadNum - 1){
			for (j = n - ((n - k - 1) % 8); j < n; j++)
				mat[k * n + j] /= mat[k * n + k];
			for (int i = 0; i < threadNum-1; ++i)
				sem_wait(&sem_leader);
			for (int i = 0; i < threadNum - 1; ++i)
				sem_post(&sem_Divsion[i]);
		}else{
			sem_post(&sem_leader);
			sem_wait(&sem_Divsion[t_id]);
		}

		//����ֵ��ֵΪ1����������if (t_id == threadNum-1)����У����?���������̻߳�û��ɳ���?
		//֮�󲢲����õ�������ݣ������ڴ˿̸�ֵ�?1�����������?
		if (t_id == 0) mat[k * n + k] = 1.0f;

		// ѭ����������ʹ���ش��¾���
		for (int i = k + 1 + t_id; i < n; i += threadNum)
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

		// �����߳�׼��������һ��
		if (t_id == threadNum - 1) {
			for (int i = 0; i < threadNum - 1; ++i)
				sem_wait(&sem_leader);
			for (int i = 0; i < threadNum - 1; ++i)
				sem_post(&sem_Elimination[i]);
		}
		else {
			sem_post(&sem_leader);
			sem_wait(&sem_Elimination[t_id]);
		}
	}
	if (t_id != 0)pthread_exit(NULL);
	return NULL;
}

// ѭ������+�ź���
void pthreadCyclicSem(int n, float* mat, int threadNum)
{
	// ��ʼ��
	sem_init(&sem_leader, 0, 0);
	for (int i = 0; i < threadNum; ++i)
	{
		sem_init(&sem_Divsion[i], 0, 0);
		sem_init(&sem_Elimination[i], 0, 0);
	}

	// �����߳�
	for (int i = threadNum - 1; i >= 0; i--)//���򴴽�����֤���һ���߳����п�������ִ�������
	{
		thread_param_t[i].t_id = i;
		thread_param_t[i].mat = mat;
		thread_param_t[i].n = n;
		thread_param_t[i].threadNum = threadNum;
		if (i != 0)pthread_create(&threads[i], NULL, threadFunc_cyclicSem, (void*)(&thread_param_t[i]));
	}
	threadFunc_cyclicSem((void*)(&thread_param_t[0]));//���߳���Ϊ�߳�0���м���

	// �ȴ������߳�ִ�����?
	for (int i = 1; i < threadNum; i++)
	{
		pthread_join(threads[i], NULL);
	}

	// ����
	sem_destroy(&sem_leader);
	for (int i = 0; i < threadNum; ++i)
	{
		sem_destroy(&sem_Divsion[i]);
		sem_destroy(&sem_Elimination[i]);
	}
}

// ��ѭ�������̺߳���
void* threadFunc_blockCyclic(void* param)
{
	__m256 vt, va, vaik, vakj, vaij;
	gaussPara* thread_param = (gaussPara*)param;
	int n = thread_param->n;
	int t_id = thread_param->t_id;
	int threadNum = thread_param->threadNum;
	float* mat = thread_param->mat;
	int chunk = 4;

	for (int k = 0; k < n; k++)
	{
		// ���������߳̽��г�������
		int step = threadNum * 8;
		vt = _mm256_set1_ps(mat[k * n + k]);
		int j = k + 1 + t_id * 8;
		for (; j + 8 <= n; j += step)
		{
			va = _mm256_loadu_ps(&mat[k * n + j]);
			va = _mm256_div_ps(va, vt);
			_mm256_storeu_ps(&mat[k * n + j], va);
		}
		if (t_id == threadNum - 1)//�߳�threadNum-1�������������µ����ݣ���Ϊ���߳����п����ڴ�ʱ����
		{
			for (j = n - ((n - k - 1) % 8); j < n; j++)
				mat[k * n + j] /= mat[k * n + k];
		}

		pthread_barrier_wait(&barrier_Divsion);

		//����ֵ��ֵΪ1����������if (t_id == threadNum-1)����У����?���������̻߳�û��ɳ���?
		//֮�󲢲����õ�������ݣ������ڴ˿̸�ֵ�?1�����������?
		if (t_id == 0) mat[k * n + k] = 1.0f;

		// ��ѭ����������ʹ���ش��¾���
		for (int outi = k + 1 + t_id * chunk; outi < n; outi += threadNum * chunk)
		{
			for(int i=outi;i<n&&i<outi+chunk;++i)
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

		// �����߳�׼��������һ��
		pthread_barrier_wait(&barrier_Elimination);
	}
	if (t_id != 0)pthread_exit(NULL);
	return NULL;
}

// ��ѭ����̬����
void pthreadBlockCyclic(int n, float* mat, int threadNum)
{
	// ��ʼ��
	pthread_barrier_init(&barrier_Divsion, NULL, threadNum);
	pthread_barrier_init(&barrier_Elimination, NULL, threadNum);

	// �����߳�
	for (int i = threadNum - 1; i >= 0; i--)//���򴴽�����֤���һ���߳����п�������ִ�������
	{
		thread_param_t[i].t_id = i;
		thread_param_t[i].mat = mat;
		thread_param_t[i].n = n;
		thread_param_t[i].threadNum = threadNum;
		if (i != 0)pthread_create(&threads[i], NULL, threadFunc_blockCyclic, (void*)(&thread_param_t[i]));
	}
	threadFunc_blockCyclic((void*)(&thread_param_t[0]));//���߳���Ϊ�߳�0���м���

	// �ȴ������߳�ִ�����?
	for (int i = 1; i < threadNum; i++)
	{
		pthread_join(threads[i], NULL);
	}

	// ����
	pthread_barrier_destroy(&barrier_Divsion);
	pthread_barrier_destroy(&barrier_Elimination);
}

// ��̬�����̺߳���
void* threadFunc_dynamic(void* param)
{
	__m256 vt, va, vaik, vakj, vaij;
	gaussPara* thread_param = (gaussPara*)param;
	int n = thread_param->n;
	int t_id = thread_param->t_id;
	int threadNum = thread_param->threadNum;
	float* mat = thread_param->mat;

	for (int k = 0; k < n; k++)
	{
		// ���������߳̽��г�������
		int step = threadNum * 8;
		vt = _mm256_set1_ps(mat[k * n + k]);
		int j = k + 1 + t_id * 8;
		for (; j + 8 <= n; j += step)
		{
			va = _mm256_loadu_ps(&mat[k * n + j]);
			va = _mm256_div_ps(va, vt);
			_mm256_storeu_ps(&mat[k * n + j], va);
		}
		if (t_id == threadNum - 1)//�߳�threadNum-1�������������µ����ݣ���Ϊ���߳����п����ڴ�ʱ����
		{
			for (j = n - ((n - k - 1) % 8); j < n; j++)
				mat[k * n + j] /= mat[k * n + k];
			myIndex = k + 1;//���ڶ�̬��������
		}

		pthread_barrier_wait(&barrier_Divsion);

		//����ֵ��ֵΪ1����������if (t_id == threadNum-1)����У����?���������̻߳�û��ɳ���?
		//֮�󲢲����õ�������ݣ������ڴ˿̸�ֵ�?1�����������?
		if (t_id == 0) mat[k * n + k] = 1.0f;

		// ��̬��������ʹ���ش��¾���
		while(myIndex<n)
		{
			pthread_mutex_lock(&mutexI);
			int i = myIndex++;
			pthread_mutex_unlock(&mutexI);
			if (i >= n) break;

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

		// �����߳�׼��������һ��
		pthread_barrier_wait(&barrier_Elimination);
	}
	if (t_id != 0)pthread_exit(NULL);
	return NULL;
}

// ��̬����
void pthreadDynamic(int n, float* mat, int threadNum)
{
	// ��ʼ��
	pthread_barrier_init(&barrier_Divsion, NULL, threadNum);
	pthread_barrier_init(&barrier_Elimination, NULL, threadNum);
	pthread_mutex_init(&mutexI, NULL);

	// �����߳�
	for (int i = threadNum - 1; i >= 0; i--)//���򴴽�����֤���һ���߳����п�������ִ�������
	{
		thread_param_t[i].t_id = i;
		thread_param_t[i].mat = mat;
		thread_param_t[i].n = n;
		thread_param_t[i].threadNum = threadNum;
		if (i != 0)pthread_create(&threads[i], NULL, threadFunc_dynamic, (void*)(&thread_param_t[i]));
	}
	threadFunc_dynamic((void*)(&thread_param_t[0]));//���߳���Ϊ�߳�0���м���

	// �ȴ������߳�ִ�����?
	for (int i = 1; i < threadNum; i++)
	{
		pthread_join(threads[i], NULL);
	}

	// ����
	pthread_barrier_destroy(&barrier_Divsion);
	pthread_barrier_destroy(&barrier_Elimination);
	pthread_mutex_destroy(&mutexI);
}

// ��̬�߳��̺߳���
void* threadFunc_dynamic_thread(void* param)
{
	gaussPara* thread_param = (gaussPara*)param;
	int n = thread_param->n;
	int k = thread_param->k;
	int i = thread_param->t_id;
	float* mat = thread_param->mat;

	// ѭ����������ʹ���ش��¾���
	__m256 vaik = _mm256_set1_ps(mat[i * n + k]);
	int j = k + 1;
	for (; j < n; ++j)
	{
		//32�ֽڶ���
		if (!((ull(&mat[i * n + j])) & 0x1f)) break;
		mat[i * n + j] -= mat[k * n + j] * mat[i * n + k];
	}
	for (; j + 8 <= n; j += 8)
	{
		__m256 vakj = _mm256_load_ps(&mat[k * n + j]);
		__m256 vaij = _mm256_load_ps(&mat[i * n + j]);
		vaij = _mm256_fnmadd_ps(vakj, vaik, vaij);
		_mm256_store_ps(&mat[i * n + j], vaij);
	}
	for (; j < n; j++) mat[i * n + j] -= mat[k * n + j] * mat[i * n + k];
	mat[i * n + k] = 0.0f;
	pthread_exit(NULL);
	return NULL;
}

// ��̬�߳�
void pthreadDynamicThread(int n, float* mat, int threadNum=1)
{
	for (int k = 0; k < n; ++k)
	{
		__m256 vt = _mm256_set1_ps(mat[k * n + k]);
		int j = k + 1;
		for (; j < n; j++)
		{
			//32�ֽڶ���
			if (!((ull(&mat[k * n + j])) & 0x1f)) break;
			mat[k * n + j] /= mat[k * n + k];
		}
		for (; j + 8 <= n; j += 8)
		{
			__m256 va = _mm256_load_ps(&mat[k * n + j]);
			va = _mm256_div_ps(va, vt);
			_mm256_store_ps(&mat[k * n + j], va);
		}
		for (; j < n; j++) mat[k * n + j] /= mat[k * n + k];
		mat[k * n + k] = 1.0f;

		// �����߳�
		for (int i = k+1; i < n; i++)//���򴴽�����֤���һ���߳����п�������ִ�������
		{
			thread_param_t[i].t_id = i;
			thread_param_t[i].mat = mat;
			thread_param_t[i].n = n;
			thread_param_t[i].k = k;
			pthread_create(&threads[i], NULL, threadFunc_dynamic_thread, (void*)(&thread_param_t[i]));
		}

		// �ȴ������߳�ִ�����?
		for (int i = k + 1; i < n; i++)
		{
			pthread_join(threads[i], NULL);
		}
	}
}

// ѭ������+��ʹ��simd�̺߳���
void* threadFunc_nosimd(void* param)
{
	gaussPara* thread_param = (gaussPara*)param;
	int n = thread_param->n;
	int t_id = thread_param->t_id;
	int threadNum = thread_param->threadNum;
	float* mat = thread_param->mat;

	for (int k = 0; k < n; k++)
	{
		// ���������߳̽��г�������
		float pivot = mat[k * n + k];
		for (int j = k + 1 + t_id; j < n; j += threadNum)mat[k * n + j] /= pivot;
		
		pthread_barrier_wait(&barrier_Divsion);

		//����ֵ��ֵΪ1����������if (t_id == threadNum-1)����У����?���������̻߳�û��ɳ���?
		//֮�󲢲����õ�������ݣ������ڴ˿̸�ֵ�?1�����������?
		if (t_id == threadNum - 1) mat[k * n + k] = 1.0f;

		// ѭ����������ʹ���ش��¾���
		for (int i = k + 1 + t_id; i < n; i += threadNum)
		{
			float pivot = mat[i * n + k];
			for (int j = k + 1; j < n; ++j)
				mat[i * n + j] -= mat[k * n + j] * pivot;
			mat[i * n + k] = 0.0f;
		}

		// �����߳�׼��������һ��
		pthread_barrier_wait(&barrier_Elimination);
	}
	if (t_id != 0)pthread_exit(NULL);
	return NULL;
}

// ѭ������+��ʹ��simd
void pthreadNoSimd(int n, float* mat, int threadNum)
{
	// ��ʼ��
	pthread_barrier_init(&barrier_Divsion, NULL, threadNum);
	pthread_barrier_init(&barrier_Elimination, NULL, threadNum);

	// �����߳�
	for (int i = threadNum - 1; i >= 0; i--)//���򴴽�����֤���һ���߳����п�������ִ�������
	{
		thread_param_t[i].t_id = i;
		thread_param_t[i].mat = mat;
		thread_param_t[i].n = n;
		thread_param_t[i].threadNum = threadNum;
		if (i != 0)pthread_create(&threads[i], NULL, threadFunc_nosimd, (void*)(&thread_param_t[i]));
	}
	threadFunc_nosimd((void*)(&thread_param_t[0]));//���߳���Ϊ�߳�0���м���

	// �ȴ������߳�ִ�����?
	for (int i = 1; i < threadNum; i++)
	{
		pthread_join(threads[i], NULL);
	}

	// ����
	pthread_barrier_destroy(&barrier_Divsion);
	pthread_barrier_destroy(&barrier_Elimination);
}

#else
// ѭ������+neon�̺߳���
void* threadFunc_neon(void* param)
{
	float32x4_t vt, va, vaik, vakj, vaij;
	gaussPara* thread_param = (gaussPara*)param;
	int n = thread_param->n;
	int t_id = thread_param->t_id;
	int threadNum = thread_param->threadNum;
	float* mat = thread_param->mat;

	for (int k = 0; k < n; k++)
	{
		// ���������߳̽��г�������
		int step = threadNum * 4;
		vt = vdupq_n_f32(mat[k * n + k]);
		int j = k + 1 + t_id * 4;
		for (; j + 4 <= n; j += step)
		{
			va = vld1q_f32(&mat[k * n + j]);
			va = vmulq_f32(va, vt);
			vst1q_f32(&mat[k * n + j], va);
		}
		if (t_id == threadNum - 1)//�߳�threadNum-1�������������µ����ݣ���Ϊ���߳����п����ڴ�ʱ����
		{
			for (j = n - ((n - k - 1) % 4); j < n; j++)
				mat[k * n + j] /= mat[k * n + k];
		}

		pthread_barrier_wait(&barrier_Divsion);

		//����ֵ��ֵΪ1����������if (t_id == threadNum-1)����У����?���������̻߳�û��ɳ���?
		//֮�󲢲����õ�������ݣ������ڴ˿̸�ֵ�?1�����������?
		if (t_id == 0) mat[k * n + k] = 1.0f;

		// ѭ����������ʹ���ش��¾���
		for (int i = k + 1 + t_id; i < n; i += threadNum)
		{
			vaik = vdupq_n_f32(mat[i * n + k]);
			int j = k + 1;
			for (; j < n; ++j)
			{
				//16�ֽڶ���
				if (!((ull(&mat[i * n + j])) & 0xf)) break;
				mat[i * n + j] -= mat[k * n + j] * mat[i * n + k];
			}
			for (; j + 8 <= n; j += 8)
			{
				vakj = vld1q_f32(&mat[k * n + j]);
				vaij = vld1q_f32(&mat[i * n + j]);
				vaij = vmlsq_f32(vakj, vaik, vaij);
				vst1q_f32(&mat[i * n + j], vaij);
			}
			for (; j < n; j++) mat[i * n + j] -= mat[k * n + j] * mat[i * n + k];
			mat[i * n + k] = 0.0f;
		}

		// �����߳�׼��������һ��
		pthread_barrier_wait(&barrier_Elimination);
	}
	if (t_id != 0)pthread_exit(NULL);
	return NULL;
}

// ѭ������+neon
void pthreadNeon(int n, float* mat, int threadNum)
{
	// ��ʼ��
	pthread_barrier_init(&barrier_Divsion, NULL, threadNum);
	pthread_barrier_init(&barrier_Elimination, NULL, threadNum);

	// �����߳�
	for (int i = threadNum - 1; i >= 0; i--)//���򴴽�����֤���һ���߳����п�������ִ�������
	{
		thread_param_t[i].t_id = i;
		thread_param_t[i].mat = mat;
		thread_param_t[i].n = n;
		thread_param_t[i].threadNum = threadNum;
		if (i != 0)pthread_create(&threads[i], NULL, threadFunc_neon, (void*)(&thread_param_t[i]));
	}
	threadFunc_neon((void*)(&thread_param_t[0]));//���߳���Ϊ�߳�0���м���

	// �ȴ������߳�ִ�����?
	for (int i = 1; i < threadNum; i++)
	{
		pthread_join(threads[i], NULL);
	}

	// ����
	pthread_barrier_destroy(&barrier_Divsion);
	pthread_barrier_destroy(&barrier_Elimination);
}
#endif