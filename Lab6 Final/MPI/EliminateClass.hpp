#pragma once
#include "config.h"
void setBit(uint* bits, int k);
void resetBit(uint* bits, int k);
bool getBit(uint* bits, int k);
void* myAlignedMalloc(size_t required_bytes, size_t alignment);
void myAlignedFree(void* p2);

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
    Eliminater(int Nmat, int Nee);
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
    Eliminatee(int Nmat,int Nee);
    Eliminatee(const Eliminatee& para);
    Eliminatee(const Eliminatee& para,int my_rank);
    ~Eliminatee();
    void print();
    void check(string path,int step);
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
        m_elierList[i] = dataVector + i * (m_Ncol>>5);

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

Eliminater::Eliminater(int Nmat, int Nee)
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

void Eliminatee::check(string path,int step)
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
        for (int index = m_Ncol -1; index >= 0; --index)
        {
            if (!getBit(eliee, index)) continue;
            if (!(ss >> temp))
            {
                printf("\ncheck result: Wrong!\n");
                printf("In row %d, it's empty row,but get index = %d", i, index);
                MPI_Finalize();
                exit(0);
            }
            if (index != temp)
            {
                printf("\ncheck result: Wrong!\n");
                printf("In row %d, ans = %d, index = %d\n", i, temp, index);
                MPI_Finalize();
                exit(0);
            }
        }
    }
    ifs.close();

    if (step == REPEAT_NUM)
        cout << "check result : Right!\n";
}

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