#define _CRT_SECURE_NO_WARNINGS
#include "config.h"
#include "EliminateClass.hpp"

int Nmat;//矩阵阶数
int Nee;//被消元行行数
int Ner;//消元子行数
int Ncol;;//转为uint后对应的列数
int comm_sz, my_rank;//MPI总进程数及当前id
int start, task_num, end;//块划分的相关信息
int startArr[8];//每个进程所分配的开始的行数
uint* elierRow;//当前的消元子
Eliminater* globalEliminater;//全局消元子类
Eliminatee* globalEliminatee;//全局被消元行类
Eliminater* eliminater;//消元子类
Eliminatee* eliminatee;//被消元行类

//读取数据
void dataInit(string DATA)
{
    //初始化信息
    stringstream ss(DATA);
    char c;
    ss >> Nmat >> c >> Nmat >> c >> Ner >> c >> Nee;
    if (DATA == "8_23045_18748_14325") Nmat = 23075;
    if (DATA == "9_37960_29304_14921") Nee = 14291;
    Ncol = (Nmat - (Nmat & 0xff) + 0x100) >> 5;
    startArr[0] = 0;
    for (int i = 1; i < comm_sz; i++) startArr[i] = startArr[i - 1] + (i - 1 < (Nee % comm_sz) ? Nee / comm_sz + 1 : Nee / comm_sz);
    task_num = my_rank < (Nee% comm_sz) ? Nee / comm_sz + 1 : Nee / comm_sz;

    //给变量申请空间
    if (my_rank == 0)
    {
        eliminater = new Eliminater(Nmat, Ner);
        eliminatee = new Eliminatee(Nmat, Nee);
    }
    else {
        eliminatee = new Eliminatee(Nmat, task_num);
    }
    elierRow = (uint*)malloc(Ncol * sizeof(uint));//该一维向量用于消元时传递枢轴所在行

    if (my_rank == 0) {
        globalEliminater = new Eliminater(DATAPATH + DATA + Para1, Nmat, Ner);
        globalEliminatee = new Eliminatee(DATAPATH + DATA + Para2, Nmat, Nee);
    }
}

//测试函数
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
                MPI_Barrier(MPI_COMM_WORLD);//其他进程等待数据初始化
                
                //程序的主体部分
                double localTime=funArr[funIndex]();

                if (my_rank == 0) {
                    eliminatee->check(DATAPATH + DATA + Para3,step);
                    totalTime += localTime;
                }
                
                if (dataIndex == 7 || dataIndex == 8 || dataIndex == 9)break;//第7/8/9号数据耗时太长，只计算一次即可
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
    if (my_rank != 0)//仅使用0号进程进行计算
        return 0;

    double start_time, end_time;
    start_time = MPI_Wtime();

    uint* elier, * eliee;
    for (int j = Nmat - 1; j >= 0; --j)
    { // 遍历消元子
        elier = eliminater->m_elierList[j];
        if (!eliminater->m_isElierList[j])
        { // 如果不存在对应消元子则将被消元行升格
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
        { // 遍历被消元行
            if (eliminatee->m_isElierList[i])
                continue;
            eliee = eliminatee->m_elieeList[i];
            if (getBit(eliee, j))
            { // 如果当前行需要消元
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
    MPI_Barrier(MPI_COMM_WORLD);//所有线程同时开始，便于计时
    start_time = MPI_Wtime();

    //首先将每个进程要处理的数据放到数组最前面
    if (my_rank == 0) {// 0号进程负责任务的初始分发工作
        for (int p = 1; p < comm_sz; p++) {
            int count = p < (Nee% comm_sz) ? Nee / comm_sz + 1 : Nee / comm_sz;
            MPI_Send(eliminatee->m_elieeList[startArr[p]], count * Ncol, MPI_UNSIGNED, p, 0, MPI_COMM_WORLD);
            MPI_Send(eliminatee->m_isElierList, count, MPI_C_BOOL, p, 1, MPI_COMM_WORLD);
        }
    }else {// 非0号进程负责任务的接收工作
        MPI_Recv(eliminatee->dataVector, task_num * Ncol, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(eliminatee->m_isElierList, task_num, MPI_C_BOOL, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    uint* elier, * eliee;
    for (int j = Nmat - 1; j >= 0; j--)
    { // 遍历消元子
        if (my_rank == 0)
        {
            elier = eliminater->m_elierList[j];
            if (!eliminater->m_isElierList[j])
            { // 如果不存在对应消元子则将被消元行升格
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
        { // 遍历被消元行
            if (eliminatee->m_isElierList[i])
                continue;
            eliee = eliminatee->m_elieeList[i];
            if (getBit(eliee, j))
            { // 如果当前行需要消元
                for (int k = 0; k < Ncol; ++k)
                    eliee[k] ^= elierRow[k];
            }
        }

        //将消元结果传递到0号进程中
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
    

     MPI_Barrier(MPI_COMM_WORLD);//等待全部线程完成
     end_time = MPI_Wtime();
     return end_time - start_time;
 }

double mpiAvx2Groebner()
{
    double start_time, end_time;
    MPI_Barrier(MPI_COMM_WORLD);//所有线程同时开始，便于计时
    start_time = MPI_Wtime();

    //首先将每个进程要处理的数据放到数组最前面
    if (my_rank == 0) {// 0号进程负责任务的初始分发工作
        for (int p = 1; p < comm_sz; p++) {
            int count = p < (Nee% comm_sz) ? Nee / comm_sz + 1 : Nee / comm_sz;
            MPI_Send(eliminatee->m_elieeList[startArr[p]], count * Ncol, MPI_UNSIGNED, p, 0, MPI_COMM_WORLD);
            MPI_Send(eliminatee->m_isElierList, count, MPI_C_BOOL, p, 1, MPI_COMM_WORLD);
        }
    }
    else {// 非0号进程负责任务的接收工作
        MPI_Recv(eliminatee->dataVector, task_num * Ncol, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(eliminatee->m_isElierList, task_num, MPI_C_BOOL, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    uint* elier, * eliee;
    __m256i velier, veliee;
    for (int j = Nmat - 1; j >= 0; j--)
    { // 遍历消元子
        if (my_rank == 0)
        {
            elier = eliminater->m_elierList[j];
            if (!eliminater->m_isElierList[j])
            { // 如果不存在对应消元子则将被消元行升格
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
        { // 遍历被消元行
            if (eliminatee->m_isElierList[i])
                continue;
            eliee = eliminatee->m_elieeList[i];
            if (getBit(eliee, j))
            { // 如果当前行需要消元
                for (int k = (eliminatee->m_Ncol >> 5) - 8; k >= 0; k -= 8)
                {
                    velier = _mm256_load_si256((__m256i*)(elierRow + k));
                    veliee = _mm256_load_si256((__m256i*)(eliee + k));
                    _mm256_store_si256((__m256i*)(eliee + k), _mm256_xor_si256(velier, veliee));
                }
            }
        }

        //将消元结果传递到0号进程中
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


    MPI_Barrier(MPI_COMM_WORLD);//等待全部线程完成
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