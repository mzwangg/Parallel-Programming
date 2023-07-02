#define _CRT_SECURE_NO_WARNINGS
#include "config.h"
#include "EliminateClass.hpp"

int Nmat;//�������
int Nee;//����Ԫ������
int Ner;//��Ԫ������
int Ncol;;//תΪuint���Ӧ������
int comm_sz, my_rank;//MPI�ܽ���������ǰid
int start, task_num, end;//�黮�ֵ������Ϣ
int startArr[8];//ÿ������������Ŀ�ʼ������
uint* elierRow;//��ǰ����Ԫ��
Eliminater* globalEliminater;//ȫ����Ԫ����
Eliminatee* globalEliminatee;//ȫ�ֱ���Ԫ����
Eliminater* eliminater;//��Ԫ����
Eliminatee* eliminatee;//����Ԫ����

//��ȡ����
void dataInit(string DATA)
{
    //��ʼ����Ϣ
    stringstream ss(DATA);
    char c;
    ss >> Nmat >> c >> Nmat >> c >> Ner >> c >> Nee;
    if (DATA == "8_23045_18748_14325") Nmat = 23075;
    if (DATA == "9_37960_29304_14921") Nee = 14291;
    Ncol = (Nmat - (Nmat & 0xff) + 0x100) >> 5;
    startArr[0] = 0;
    for (int i = 1; i < comm_sz; i++) startArr[i] = startArr[i - 1] + (i - 1 < (Nee % comm_sz) ? Nee / comm_sz + 1 : Nee / comm_sz);
    task_num = my_rank < (Nee% comm_sz) ? Nee / comm_sz + 1 : Nee / comm_sz;

    //����������ռ�
    if (my_rank == 0)
    {
        eliminater = new Eliminater(Nmat, Ner);
        eliminatee = new Eliminatee(Nmat, Nee);
    }
    else {
        eliminatee = new Eliminatee(Nmat, task_num);
    }
    elierRow = (uint*)malloc(Ncol * sizeof(uint));//��һά����������Ԫʱ��������������

    if (my_rank == 0) {
        globalEliminater = new Eliminater(DATAPATH + DATA + Para1, Nmat, Ner);
        globalEliminatee = new Eliminatee(DATAPATH + DATA + Para2, Nmat, Nee);
    }
}

//���Ժ���
void timingAll(string nameArr[], double(*funArr[])(), int funNum)
{
    if (my_rank == 0) {
        for (int i = 0; i < funNum; ++i)cout << nameArr[i] << '\t';
        cout << endl;
    }
    
    for (int dataIndex = 0; dataIndex < dataSize; ++dataIndex)
    {
        string DATA = dataArr[dataIndex];
        dataInit(DATA);
        MPI_Barrier(MPI_COMM_WORLD);
        
        for (int funIndex = 0; funIndex < funNum; ++funIndex)
        {
            int step = 1;
            double totalTime = 0.0f;
            for (; step <= REPEAT_NUM; ++step)
            {
                if (my_rank == 0) {
                    eliminater->copy(*globalEliminater);
                    eliminatee->copy(*globalEliminatee);
                }
                MPI_Barrier(MPI_COMM_WORLD);//�������̵ȴ����ݳ�ʼ��
                
                //��������岿��
                double localTime=funArr[funIndex]();

                if (my_rank == 0) {
                    eliminatee->check(DATAPATH + DATA + Para3,step);
                    totalTime += localTime;
                }
                
                if (dataIndex == 7 || dataIndex == 8 || dataIndex == 9)break;//��7/8/9�����ݺ�ʱ̫����ֻ����һ�μ���
            }
            if (my_rank == 0)cout << totalTime / (double)step << '\t';
        }

        free(elierRow);
        delete eliminatee;
        if (my_rank == 0) {
            delete eliminater;
            delete globalEliminater;
            delete globalEliminatee;
            cout << '\n';
        }
    }
}

double serialGroebner()
{
    if (my_rank != 0)//��ʹ��0�Ž��̽��м���
        return 0;

    double start_time, end_time;
    start_time = MPI_Wtime();

    uint* elier, * eliee;
    for (int j = Nmat - 1; j >= 0; --j)
    { // ������Ԫ��
        elier = eliminater->m_elierList[j];
        if (!eliminater->m_isElierList[j])
        { // ��������ڶ�Ӧ��Ԫ���򽫱���Ԫ������
            int i;
            for (i = 0; i < Nee; ++i)
            {
                if (eliminatee->m_isElierList[i])
                    continue;
                uint* eliee = eliminatee->m_elieeList[i];
                if (getBit(eliee, j))
                {
                    eliminater->m_elierList[j] = eliee;
                    eliminater->m_isElierList[j] = true;
                    eliminatee->m_isElierList[i] = true;
                    elier = eliee;
                    break;
                }
            }
            if (i == Nee) continue;
        }
        for (int i = 0; i < Nee; ++i)
        { // ��������Ԫ��
            if (eliminatee->m_isElierList[i])
                continue;
            eliee = eliminatee->m_elieeList[i];
            if (getBit(eliee, j))
            { // �����ǰ����Ҫ��Ԫ
                for (int k = 0; k < Ncol; ++k)
                    eliee[k] ^= elier[k];
            }
        }
    }

    end_time = MPI_Wtime();
    return end_time - start_time;
}

double mpiGroebner()
 {
    double start_time, end_time;
    MPI_Barrier(MPI_COMM_WORLD);//�����߳�ͬʱ��ʼ�����ڼ�ʱ
    start_time = MPI_Wtime();

    //���Ƚ�ÿ������Ҫ��������ݷŵ�������ǰ��
    if (my_rank == 0) {// 0�Ž��̸�������ĳ�ʼ�ַ�����
        for (int p = 1; p < comm_sz; p++) {
            int count = p < (Nee% comm_sz) ? Nee / comm_sz + 1 : Nee / comm_sz;
            MPI_Send(eliminatee->m_elieeList[startArr[p]], count * Ncol, MPI_UNSIGNED, p, 0, MPI_COMM_WORLD);
            MPI_Send(eliminatee->m_isElierList, count, MPI_C_BOOL, p, 1, MPI_COMM_WORLD);
        }
    }else {// ��0�Ž��̸�������Ľ��չ���
        MPI_Recv(eliminatee->dataVector, task_num * Ncol, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(eliminatee->m_isElierList, task_num, MPI_C_BOOL, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    uint* elier, * eliee;
    for (int j = Nmat - 1; j >= 0; j--)
    { // ������Ԫ��
        if (my_rank == 0)
        {
            elier = eliminater->m_elierList[j];
            if (!eliminater->m_isElierList[j])
            { // ��������ڶ�Ӧ��Ԫ���򽫱���Ԫ������
                for (int i = 0; i < Nee; ++i)
                {
                    if (eliminatee->m_isElierList[i])
                        continue;
                    uint* eliee = eliminatee->m_elieeList[i];
                    if (getBit(eliee, j))
                    {
                        memcpy(eliminater->m_elierList[j], eliee, Ncol * sizeof(uint));
                        eliminater->m_isElierList[j] = true;
                        eliminatee->m_isElierList[i] = true;
                        elier = eliee;
                        break;
                    }
                }
            }
            memcpy(elierRow, elier, Ncol * sizeof(uint));
            for (int p = 1; p < comm_sz; p++) {
                int count = p < (Nee% comm_sz) ? Nee / comm_sz + 1 : Nee / comm_sz;
                MPI_Send(eliminatee->m_isElierList+startArr[p], count, MPI_C_BOOL, p, 2, MPI_COMM_WORLD);
            }
        }else {
            MPI_Recv(eliminatee->m_isElierList, task_num, MPI_C_BOOL, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        MPI_Bcast(elierRow, Ncol, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        for (int i = 0; i < task_num; ++i)
        { // ��������Ԫ��
            if (eliminatee->m_isElierList[i])
                continue;
            eliee = eliminatee->m_elieeList[i];
            if (getBit(eliee, j))
            { // �����ǰ����Ҫ��Ԫ
                for (int k = 0; k < Ncol; ++k)
                    eliee[k] ^= elierRow[k];
            }
        }

        //����Ԫ������ݵ�0�Ž�����
        if (my_rank == 0) {
            for (int p = 1; p < comm_sz; p++) {
                int count = p < (Nee% comm_sz) ? Nee / comm_sz + 1 : Nee / comm_sz;
                MPI_Recv(eliminatee->m_elieeList[startArr[p]], count * Ncol, MPI_UNSIGNED, p, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        else {
            MPI_Send(eliminatee->dataVector, task_num * Ncol, MPI_UNSIGNED, 0, 3, MPI_COMM_WORLD);
        }
    }
    

     MPI_Barrier(MPI_COMM_WORLD);//�ȴ�ȫ���߳����
     end_time = MPI_Wtime();
     return end_time - start_time;
 }

double mpiAvx2Groebner()
{
    double start_time, end_time;
    MPI_Barrier(MPI_COMM_WORLD);//�����߳�ͬʱ��ʼ�����ڼ�ʱ
    start_time = MPI_Wtime();

    //���Ƚ�ÿ������Ҫ��������ݷŵ�������ǰ��
    if (my_rank == 0) {// 0�Ž��̸�������ĳ�ʼ�ַ�����
        for (int p = 1; p < comm_sz; p++) {
            int count = p < (Nee% comm_sz) ? Nee / comm_sz + 1 : Nee / comm_sz;
            MPI_Send(eliminatee->m_elieeList[startArr[p]], count * Ncol, MPI_UNSIGNED, p, 0, MPI_COMM_WORLD);
            MPI_Send(eliminatee->m_isElierList, count, MPI_C_BOOL, p, 1, MPI_COMM_WORLD);
        }
    }
    else {// ��0�Ž��̸�������Ľ��չ���
        MPI_Recv(eliminatee->dataVector, task_num * Ncol, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(eliminatee->m_isElierList, task_num, MPI_C_BOOL, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    uint* elier, * eliee;
    __m256i velier, veliee;
    for (int j = Nmat - 1; j >= 0; j--)
    { // ������Ԫ��
        if (my_rank == 0)
        {
            elier = eliminater->m_elierList[j];
            if (!eliminater->m_isElierList[j])
            { // ��������ڶ�Ӧ��Ԫ���򽫱���Ԫ������
                for (int i = 0; i < Nee; ++i)
                {
                    if (eliminatee->m_isElierList[i])
                        continue;
                    uint* eliee = eliminatee->m_elieeList[i];
                    if (getBit(eliee, j))
                    {
                        memcpy(eliminater->m_elierList[j], eliee, Ncol * sizeof(uint));
                        eliminater->m_isElierList[j] = true;
                        eliminatee->m_isElierList[i] = true;
                        elier = eliee;
                        break;
                    }
                }
            }
            memcpy(elierRow, elier, Ncol * sizeof(uint));
            for (int p = 1; p < comm_sz; p++) {
                int count = p < (Nee% comm_sz) ? Nee / comm_sz + 1 : Nee / comm_sz;
                MPI_Send(eliminatee->m_isElierList + startArr[p], count, MPI_C_BOOL, p, 2, MPI_COMM_WORLD);
            }
        }
        else {
            MPI_Recv(eliminatee->m_isElierList, task_num, MPI_C_BOOL, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        MPI_Bcast(elierRow, Ncol, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        for (int i = 0; i < task_num; ++i)
        { // ��������Ԫ��
            if (eliminatee->m_isElierList[i])
                continue;
            eliee = eliminatee->m_elieeList[i];
            if (getBit(eliee, j))
            { // �����ǰ����Ҫ��Ԫ
                for (int k = (eliminatee->m_Ncol >> 5) - 8; k >= 0; k -= 8)
                {
                    velier = _mm256_load_si256((__m256i*)(elierRow + k));
                    veliee = _mm256_load_si256((__m256i*)(eliee + k));
                    _mm256_store_si256((__m256i*)(eliee + k), _mm256_xor_si256(velier, veliee));
                }
            }
        }

        //����Ԫ������ݵ�0�Ž�����
        if (my_rank == 0) {
            for (int p = 1; p < comm_sz; p++) {
                int count = p < (Nee% comm_sz) ? Nee / comm_sz + 1 : Nee / comm_sz;
                MPI_Recv(eliminatee->m_elieeList[startArr[p]], count * Ncol, MPI_UNSIGNED, p, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        else {
            MPI_Send(eliminatee->dataVector, task_num * Ncol, MPI_UNSIGNED, 0, 3, MPI_COMM_WORLD);
        }
    }


    MPI_Barrier(MPI_COMM_WORLD);//�ȴ�ȫ���߳����
    end_time = MPI_Wtime();
    return end_time - start_time;
}

 int main()
 {
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    //double(*funArr1[])() = { serialGroebner, mpiGroebner ,mpiAvx2Groebner };
    //string nameArr1[] = { "serialGroebner","mpiGroebner","mpiAvx2Groebner" };
    //timingAll(nameArr1, funArr1, 3);

    double(*funArr[])() = { mpiGroebner, };
    string nameArr[] = { "mpiGroebner" };
    timingAll(nameArr, funArr, 1);

    MPI_Finalize();
 }