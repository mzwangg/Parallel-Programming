#pragma once
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
    uint** m_elierList;
    Eliminater(string path, int Nmat, int Ner);
    Eliminater(const Eliminater& para);
    ~Eliminater();
};

class Eliminatee
{
public:
    int m_Nrow;
    int m_Ncol;
    bool* m_isElierList;
    uint** m_elieeList;
    Eliminatee(string path, int Nmat, int Nee);
    Eliminatee(const Eliminatee& para);
    ~Eliminatee();
    void print();
    void check(string path);
    void clear();
};

Eliminater::Eliminater(string path, int Nmat, int Ner)
{
    m_Nrow = Nmat;
    m_Ncol = Nmat - (Nmat & 0xff) + 0x100;
    m_elierList = (uint**)malloc(m_Nrow * sizeof(uint*));
    for (int i = 0; i < m_Nrow; ++i)
        m_elierList[i] = nullptr;

    ifstream ifs(path, ios::in);
    int temp, header;
    string line;
    uint* elier;

    for (int i = 0; i < Ner; ++i)
    {
        getline(ifs, line);
        stringstream ss(line);
        ss >> header;
        m_elierList[header] = (uint*)myAlignedMalloc(m_Ncol>>3, 32);
        elier = m_elierList[header];
        memset(elier, 0, m_Ncol>>3);
        setBit(elier, header);
        while (ss >> temp)
            setBit(elier, temp);
    }
    ifs.close();
}

Eliminater::Eliminater(const Eliminater& para)
{
    m_Nrow = para.m_Nrow;
    m_Ncol = para.m_Ncol;
    m_elierList = (uint**)malloc(m_Nrow * sizeof(uint*));
    for (int i = 0; i < m_Nrow; ++i)
    {
        if (para.m_elierList[i] == nullptr)m_elierList[i] = nullptr;
        else
        {
            m_elierList[i]= (uint*)myAlignedMalloc(m_Ncol >> 3, 32);
            memcpy(m_elierList[i], para.m_elierList[i], m_Ncol >> 3);
        }
    }
}

Eliminater::~Eliminater()
{
    for (int i = 0; i < m_Nrow; ++i)
    {
        if (!m_elierList[i])continue;
        myAlignedFree(m_elierList[i]);
    }
    free(m_elierList);
}

Eliminatee::Eliminatee(string path, int Nmat, int Nee)
{
    m_Nrow = Nee;
    m_Ncol = Nmat - (Nmat & 0xff) + 0x100;
    m_elieeList = (uint**)malloc(m_Nrow * sizeof(uint*));
    m_isElierList = (bool*)malloc(m_Nrow);
    for (int i = 0; i < m_Nrow; ++i)
        m_elieeList[i] = nullptr;
    memset(m_isElierList, 0, m_Nrow);

    ifstream ifs(path, ios::in);
    int temp;
    string line;

    for (int i = 0; i < m_Nrow; ++i)
    {
        m_elieeList[i] = (uint*)myAlignedMalloc(m_Ncol >> 3, 32);
        memset(m_elieeList[i], 0, m_Ncol>>3);

        getline(ifs, line);
        stringstream ss(line);

        while (ss >> temp)
            setBit(m_elieeList[i], temp);
    }
    ifs.close();
}

Eliminatee::Eliminatee(const Eliminatee& para)
{
    m_Nrow = para.m_Nrow;
    m_Ncol = para.m_Ncol;
    m_elieeList = (uint**)malloc(m_Nrow * sizeof(uint*));
    m_isElierList = (bool*)malloc(m_Nrow);
    memcpy(m_isElierList, para.m_isElierList, m_Nrow);
    for (int i = 0; i < m_Nrow; ++i)
    {
        if (para.m_elieeList[i] == nullptr)m_elieeList[i] = nullptr;
        else
        {
            m_elieeList[i] = (uint*)myAlignedMalloc(m_Ncol >> 3, 32);
            memcpy(m_elieeList[i], para.m_elieeList[i], m_Ncol >> 3);
        }
    }
}

Eliminatee::~Eliminatee()
{
    for (int i = 0; i < m_Nrow; ++i)
    {
        if (m_isElierList[i])continue;
        myAlignedFree(m_elieeList[i]);
    }
    free(m_elieeList);
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
    ifstream ifs(path, ios::in);
    string line;
    int temp, index;
    uint* eliee;

    for (int i = 0; i < m_Nrow; ++i)
    {
        index = m_Ncol;
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
}

//ÖÃ1
void setBit(uint* bits, int k)
{
    bits[k >> 5] |= (0x80000000 >> (k & 0x1f));
}

//ÖÃ0
void resetBit(uint* bits, int k)
{
    bits[k >> 5] &= ~(0x80000000 >> (k & 0x1f));
}

//·ÃÎÊ
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