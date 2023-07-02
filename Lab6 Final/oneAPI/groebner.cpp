#include <cstring>
#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <chrono>
#include <CL/sycl.hpp>
#include "systime.h"
using namespace std;

#define ull (unsigned long long)
#define uint unsigned int
//#define CHECK

#ifdef _WIN32
const string DATAPATH = "..\\data\\Groebner\\";
const string Para1 = "\\消元子.txt";
const string Para2 = "\\被消元行.txt";
const string Para3 = "\\消元结果.txt";
#elif __ARM_NEON
const string DATAPATH = "/home/data/Groebner/";
const string Para1 = "/1.txt";
const string Para2 = "/2.txt";
const string Para3 = "/3,txt";
#else
const string DATAPATH = "./data/Groebner\\";
const string Para1 = "/消元子.txt";
const string Para2 = "/被消元行.txt";
const string Para3 = "/消元结果,txt";
#endif

const float eps = 1e-9;
const int floatLimit = (1 << 24);
const int BLOCK_SIZE = 1024;
const int REPEAT_NUM = 6;
const int WARMUP = 2;
const int groupSize = 8;
const int gridSize = 16;
const int blockSize = 1024;
const string dataArr[] = { "1_130_22_8","2_254_106_53","3_562_170_53","4_1011_539_263","5_2362_1226_453","6_3799_2759_1953",\
                        "7_8399_6375_4535","8_23045_18748_14325","9_37960_29304_14921","10_43577_39477_54274","11_85401_5724_756" };


//*******************消元子和消元行类*********************88


//置1
void setBit(uint* bits, int k)
{
    bits[k >> 5] |= (0x80000000 >> (k & 0x1f));
}

//置0
void resetBit(uint* bits, int k)
{
    bits[k >> 5] &= ~(0x80000000 >> (k & 0x1f));
}

//访问
bool getBit(uint* bits, int k)
{
    return bits[k >> 5] & (0x80000000 >> (k & 0x1f));
}

void* myAlignedMalloc(size_t required_bytes, size_t alignment)
{
    int offset = alignment - 1 + sizeof(void*);
    void* p1 = (void*)malloc(required_bytes + offset);
    if (p1 == NULL)
        return NULL;
    void** p2 = (void**)(((size_t)p1 + offset) & ~(alignment - 1));
    p2[-1] = p1;
    return p2;
}

void myAlignedFree(void* p2)
{
    void* p1 = ((void**)p2)[-1];
    free(p1);
}

class Eliminater
{
public:
    int m_Nrow;
    int m_Ncol;
    bool* m_isElierList;
    uint* dataVector;
    uint** m_elierList;
    void copy(const Eliminater& para);
    Eliminater(string path, int Nmat, int Ner);
    Eliminater(int Nmat);
    Eliminater(const Eliminater& para);
    ~Eliminater();
};

class Eliminatee
{
public:
    int m_Nrow;
    int m_Ncol;
    bool* m_isElierList;
    uint* dataVector;
    uint** m_elieeList;
    void copy(const Eliminatee& para);
    Eliminatee(string path, int Nmat, int Nee);
    Eliminatee(int Nmat, int Nee);
    Eliminatee(const Eliminatee& para);
    Eliminatee(const Eliminatee& para, int my_rank);
    ~Eliminatee();
    void print();
    void check(string path);
    void clear();
};

Eliminater::Eliminater(string path, int Nmat, int Ner)
{
    m_Nrow = Nmat;
    m_Ncol = Nmat - (Nmat & 0xff) + 0x100;
    m_isElierList = (bool*)malloc(m_Nrow);
    memset(m_isElierList, 0, m_Nrow);
    m_elierList = (uint**)malloc(m_Nrow * sizeof(uint*));
    dataVector = (uint*)myAlignedMalloc(m_Nrow * (m_Ncol >> 3), 32);
    memset(dataVector, 0, m_Nrow * (m_Ncol >> 3));
    for (int i = 0; i < m_Nrow; ++i)
        m_elierList[i] = dataVector + i * (m_Ncol >> 5);

    ifstream ifs(path, ios::in);
    int temp, header;
    string line;
    uint* elier;

    for (int i = 0; i < Ner; ++i)
    {
        getline(ifs, line);
        stringstream ss(line);
        ss >> header;
        m_isElierList[header] = true;
        elier = m_elierList[header];
        setBit(elier, header);
        while (ss >> temp)
            setBit(elier, temp);
    }
    ifs.close();
}

Eliminater::Eliminater(int Nmat)
{
    m_Nrow = Nmat;
    m_Ncol = Nmat - (Nmat & 0xff) + 0x100;
    m_isElierList = (bool*)malloc(m_Nrow);
    m_elierList = (uint**)malloc(m_Nrow * sizeof(uint*));
    dataVector = (uint*)myAlignedMalloc(m_Nrow * (m_Ncol >> 3), 32);
    for (int i = 0; i < m_Nrow; ++i)
        m_elierList[i] = dataVector + i * (m_Ncol >> 5);
}

Eliminater::Eliminater(const Eliminater& para)
{
    m_Nrow = para.m_Nrow;
    m_Ncol = para.m_Ncol;
    m_isElierList = (bool*)malloc(m_Nrow);
    memcpy(m_isElierList, para.m_isElierList, m_Nrow);
    m_elierList = (uint**)malloc(m_Nrow * sizeof(uint*));
    dataVector = (uint*)myAlignedMalloc(m_Nrow * (m_Ncol >> 3), 32);
    memcpy(dataVector, para.dataVector, m_Nrow * (m_Ncol >> 3));
    for (int i = 0; i < m_Nrow; ++i)
        m_elierList[i] = dataVector + i * (m_Ncol >> 5);
}

void Eliminater::copy(const Eliminater& para)
{
    memcpy(m_isElierList, para.m_isElierList, m_Nrow);
    memcpy(dataVector, para.dataVector, m_Nrow * (m_Ncol >> 3));
}

Eliminater::~Eliminater()
{
    myAlignedFree(dataVector);
    free(m_elierList);
    free(m_isElierList);
}

Eliminatee::Eliminatee(string path, int Nmat, int Nee)
{
    m_Nrow = Nee;
    m_Ncol = Nmat - (Nmat & 0xff) + 0x100;
    m_isElierList = (bool*)malloc(m_Nrow);
    memset(m_isElierList, 0, m_Nrow);
    m_elieeList = (uint**)malloc(m_Nrow * sizeof(uint*));
    dataVector = (uint*)myAlignedMalloc(m_Nrow * (m_Ncol >> 3), 32);
    memset(dataVector, 0, m_Nrow * (m_Ncol >> 3));
    for (int i = 0; i < m_Nrow; ++i)
        m_elieeList[i] = dataVector + i * (m_Ncol >> 5);

    ifstream ifs(path, ios::in);
    int temp;
    string line;

    for (int i = 0; i < m_Nrow; ++i)
    {
        getline(ifs, line);
        stringstream ss(line);
        while (ss >> temp)
            setBit(m_elieeList[i], temp);
    }
    ifs.close();
}

Eliminatee::Eliminatee(int Nmat, int Nee)
{
    m_Nrow = Nee;
    m_Ncol = Nmat - (Nmat & 0xff) + 0x100;
    m_isElierList = (bool*)malloc(m_Nrow);
    m_elieeList = (uint**)malloc(m_Nrow * sizeof(uint*));
    dataVector = (uint*)myAlignedMalloc(m_Nrow * (m_Ncol >> 3), 32);
    for (int i = 0; i < m_Nrow; ++i)
        m_elieeList[i] = dataVector + i * (m_Ncol >> 5);
}

Eliminatee::Eliminatee(const Eliminatee& para)
{
    m_Nrow = para.m_Nrow;
    m_Ncol = para.m_Ncol;
    m_isElierList = (bool*)malloc(m_Nrow);
    memcpy(m_isElierList, para.m_isElierList, m_Nrow);
    m_elieeList = (uint**)malloc(m_Nrow * sizeof(uint*));
    dataVector = (uint*)myAlignedMalloc(m_Nrow * (m_Ncol >> 3), 32);
    memcpy(dataVector, para.dataVector, m_Nrow * (m_Ncol >> 3));
    for (int i = 0; i < m_Nrow; ++i)
        m_elieeList[i] = dataVector + i * (m_Ncol >> 5);
}

void Eliminatee::copy(const Eliminatee& para)
{
    memcpy(m_isElierList, para.m_isElierList, m_Nrow);
    memcpy(dataVector, para.dataVector, m_Nrow * (m_Ncol >> 3));
}

Eliminatee::~Eliminatee()
{
    myAlignedFree(dataVector);
    free(m_elieeList);
    free(m_isElierList);
}

void Eliminatee::print()
{
    for (int i = 0; i < m_Nrow; ++i)
    {
        uint* eliee = m_elieeList[i];
        for (int j = m_Ncol << 2; j >= 0; --j)
        {
            if (getBit(eliee, j)) printf("%d ", j);
        }
        printf("\n");
    }
}

void Eliminatee::check(string path)
{
    //判断是否需要检查
#ifndef CHECK
    return;
#endif

    ifstream ifs(path, ios::in);
    string line;
    int temp;
    uint* eliee;

    for (int i = 0; i < m_Nrow; ++i)
    {
        eliee = m_elieeList[i];
        getline(ifs, line);
        stringstream ss(line);
        for (int index = m_Ncol - 1; index >= 0; --index)
        {
            if (!getBit(eliee, index)) continue;
            if (!(ss >> temp))
            {
                printf("\ncheck result: Wrong!\n");
                printf("In row %d, it's empty row,but get index = %d", i, index);
                exit(0);
            }
            if (index != temp)
            {
                printf("\ncheck result: Wrong!\n");
                printf("In row %d, ans = %d, index = %d\n", i, temp, index);
                exit(0);
            }
        }
    }
    ifs.close();

    cout << "check result : Right!\n";
}

//***********************全局变量*******************************

int Nmat;//矩阵阶数
int Nee;//被消元行行数
int Ner;//消元子行数
int Ncol;//转为uint后对应的列数

Eliminater* globalEliminater;//全局消元子类
Eliminatee* globalEliminatee;//全局被消元行类
Eliminater* eliminater;//消元子类
Eliminatee* eliminatee;//被消元行类

//***********************基础函数*******************************

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

    eliminater = new Eliminater(Nmat);
    eliminatee = new Eliminatee(Nmat, Nee);
    globalEliminater = new Eliminater(DATAPATH + DATA + Para1, Nmat, Ner);
    globalEliminatee = new Eliminatee(DATAPATH + DATA + Para2, Nmat, Nee);
}

double serial1()
{
    auto start = chrono::high_resolution_clock::now();

    uint* elier, * eliee;
    for (int i = 0; i < Nee; ++i)
    {
        eliee = eliminatee->m_elieeList[i];
        for (int j = Nmat - 1; j >= 0; --j)
        {
            if (!getBit(eliee, j)) continue;
            if (eliminater->m_isElierList[j])
            {
                elier = eliminater->m_elierList[j];
                for (int k = 0; k < Ncol; ++k) eliee[k] ^= elier[k];
            }
            else
            {
                eliminater->m_elierList[j] = eliee;
                eliminater->m_isElierList[j] = true;
                eliminatee->m_isElierList[i] = true;
                break;
            }
        }
    }

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start) / 1000000.0;
    return duration.count();
}

double serial2()
{
    struct timeval start, end;
    gettimeofday(&start, NULL);

    uint* elier, * eliee;
    for (int j = Nmat - 1; j >= 0; --j) { // 遍历消元子
        elier = eliminater->m_elierList[j];
        if (!eliminater->m_isElierList[j]) { // 如果不存在对应消元子则将被消元行升格
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

    gettimeofday(&end, NULL);
    return (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
}

double serial3()
{
    struct timeval start, end;
    gettimeofday(&start, NULL);

    int endj = Nmat - 1;
    uint* elier, * eliee;
    for (int j = endj; j >= 0; j = endj) { // 遍历消元子

        //以groupSize为最大个数找到可以直接连续消元的行数
        for (int i = 0; i < groupSize; ++i, --endj)
            if (endj < 0 || !eliminater->m_isElierList[endj])
                break;

        if (j == endj) {// 如果不存在对应消元子则将被消元行升格
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
                    break;
                }
            }
            if (i == Nee) --endj;
            continue;
        }

        for (int jj = j; jj > endj; --jj) {
            elier = eliminater->m_elierList[jj];
            for (int i = 0; i < Nee; ++i) { // 遍历被消元行
                if (eliminatee->m_isElierList[i])
                    continue;
                eliee = eliminatee->m_elieeList[i];
                if (getBit(eliee, jj))
                { // 如果当前行需要消元
                    for (int k = 0; k < Ncol; ++k)
                        eliee[k] ^= elier[k];
                }
            }
        }
    }

    gettimeofday(&end, NULL);
    return (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
}

//***********************GPU*******************************

double eliminate_kernel(uint* elieeData, uint* eliers, bool* isElierList, sycl::queue& q)
{
    int m_Nmat = Nmat;//矩阵阶数
    int m_Nee = Nee;//被消元行行数
    int m_Ner = Ner;//消元子行数
    int m_Ncol = Ncol;//转为uint后对应的列数

    struct timeval start, end;
    gettimeofday(&start, NULL);

    //将消元子传递到GPU内存中
    memcpy(elieeData, eliminatee->dataVector, m_Nee * Ncol * sizeof(uint));

    int endj = m_Nmat - 1;
    uint* elier, * eliee;
    for (int j = endj; j >= 0; j = endj) { // 遍历消元子

        //以groupSize为最大个数找到可以直接连续消元的行数
        for (int i = 0; i < groupSize; ++i, --endj)
            if (endj < 0 || !eliminater->m_isElierList[endj])
                break;

        if (j == endj) {// 如果不存在对应消元子则将被消元行升格
            int i;
            for (i = 0; i < m_Nee; ++i)
            {
                if (eliminatee->m_isElierList[i])
                    continue;
                uint* eliee = eliminatee->m_elieeList[i];
                if (getBit(eliee, j))
                {
                    memcpy(eliminater->m_elierList[j], eliee, m_Ncol * sizeof(uint));
                    eliminater->m_isElierList[j] = true;
                    eliminatee->m_isElierList[i] = true;
                    break;
                }
            }
            if (i == m_Nee) --endj;
            continue;
        }

        //将消元子和是否升格的信息传递到GPU中
        memcpy(eliers, eliminater->dataVector + (endj + 1) * Ncol, (j - endj) * Ncol * sizeof(uint));
        memcpy(isElierList, eliminatee->m_isElierList, Ncol);

        auto e = q.submit([&](sycl::handler& h) {
            h.parallel_for(sycl::nd_range<2>({ gridSize,blockSize }, {1, blockSize }), [=](sycl::nd_item<2> item) {
                for (int jj = j; jj > endj; --jj) {
                    uint* elier = eliers + (jj - endj - 1) * m_Ncol;
                    for (int i = item.get_global_id()[0]; i < m_Nee; i += gridSize) {
                        if (isElierList[i])
                            continue;
                        uint* eliee = elieeData + i * m_Ncol;
                        if (eliee[jj >> 5] & (0x80000000 >> (jj & 0x1f)))
                        { // 如果当前行需要消元
                            for (int k = item.get_local_id()[1]; k < m_Ncol; k += blockSize)//通过循环划分划分数据
                                eliee[k] ^= elier[k];
                        }
                    }
                    item.barrier();
                }
            });
        });
        e.wait();

        if (!eliminater->m_isElierList[endj]) {//判断是否将被消元行复制回CPU中
            memcpy(eliminatee->dataVector, elieeData, m_Nee * m_Ncol * sizeof(uint));
        }
    }

    gettimeofday(&end, NULL);
    return (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
}

void gemm(int start, int end,int iterations, sycl::queue& q) {
    for (int dataIndex = start; dataIndex <= end; ++dataIndex) {
        string DATA = dataArr[dataIndex - 1];
        dataInit(DATA);

        //申请GPU内存
        auto isElierList = sycl::malloc_device<bool>(Nee, q);
        auto elieeData = sycl::malloc_device<uint>(Nee * Ncol, q);
        auto eliers = sycl::malloc_device<uint>(groupSize * Ncol, q);

        //存储时间
        double duration_cpu = 0.0f;
        double duration_gpu = 0.0f;

        // 时间测量
        //cpu串行高斯消元
        for (int run = 0; run < iterations + WARMUP; run++) {
            eliminater->copy(*globalEliminater);
            eliminatee->copy(*globalEliminatee);
            float duration = serial3();
            if (run >= WARMUP) duration_cpu += duration;//先进行热身，不计入耗时
        }
        duration_cpu = duration_cpu / iterations;
        eliminatee->check(DATAPATH + DATA + Para3);

        //GPU并行高斯消元
        for (int run = 0; run < iterations + WARMUP; run++) {
            eliminater->copy(*globalEliminater);
            eliminatee->copy(*globalEliminatee);
            float duration = eliminate_kernel(elieeData, eliers, isElierList, q);
            if (run >= WARMUP) duration_gpu += duration;//先进行热身，不计入耗时
        }
        duration_gpu = duration_gpu / iterations;
        memcpy(eliminatee->dataVector, elieeData, Nee * Ncol * sizeof(uint));
        eliminatee->check(DATAPATH + DATA + Para3);

        std::printf("%d\t%lf\t%lf\n", dataIndex, duration_cpu, duration_gpu);

        sycl::free(isElierList, q);
        sycl::free(elieeData, q);
        sycl::free(eliers, q);
        delete eliminatee;
        delete eliminater;
        delete globalEliminater;
        delete globalEliminatee;
    }
}

int main() {
    freopen("out.xls", "w", stdout);
    auto propList = cl::sycl::property_list{ cl::sycl::property::queue::enable_profiling() };
    sycl::queue my_gpu_queue(cl::sycl::cpu_selector_v, propList);

    //打印设备信息
    //cout << "Select device: " << my_gpu_queue.get_device().get_info<sycl::info::device::name>() << "\n";
    cout << "dataIndex\tcpu_time\tgpu_time\n";

    //for (int N = 64; N <= 1024; N += 64)
    gemm(11, 11, REPEAT_NUM, my_gpu_queue);

    return 0;
}