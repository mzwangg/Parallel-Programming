#include "config.h"

void setBit(uint* bits, int k);
void resetBit(uint* bits, int k);
bool getBit(uint* bits, int k);

class Eliminater
{
public:
    int m_Nrow;
    int m_Ncol;
    uint** m_elierList;
    Eliminater(string path, int Nmat, int Ner);
};

class Eliminatee
{
public:
    int m_Nrow;
    int m_Ncol;
    uint** m_elieeList;
    Eliminatee(string path, int Nmat, int Nee);
    void print();
    void check(string path);
};

Eliminater::Eliminater(string path, int Nmat, int Ner)
{
    m_Nrow = Nmat;
    m_Ncol = (Nmat & 0xffffff00) + 256;
    m_elierList = new uint* [m_Nrow];
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
        m_elierList[header] = (uint*)aligned_alloc(256, m_Ncol);//每一个消元子都申请一块256位对齐的空间
        //m_elierList[header] = new uint[m_Ncol>>5];
        elier = m_elierList[header];
        memset(elier, 0, m_Ncol);
        setBit(elier, header);
        while (ss >> temp)
            setBit(elier, temp);
    }
    ifs.close();
}

Eliminatee::Eliminatee(string path, int Nmat, int Nee)
{
    m_Nrow = Nee;
    //m_Ncol = Nmat+(256-(Nmat&0xff));//被消元行的列数需要与256位对齐
    m_Ncol = (Nmat & 0xffffff00) + 256;
    m_elieeList = new uint* [m_Nrow];

    ifstream ifs(path, ios::in);
    int temp;
    string line;

    for (int i = 0; i < m_Nrow; ++i)
    {
        m_elieeList[i] = (uint*)aligned_alloc(256, m_Ncol);//每一个被消元行都申请一块256位对齐的空间
        //m_elieeList[i] = new uint[m_Ncol>>5];
        memset(m_elieeList[i], 0, m_Ncol);

        getline(ifs, line);
        stringstream ss(line);

        while (ss >> temp)
            setBit(m_elieeList[i], temp);
    }
    ifs.close();
}

void Eliminatee::print()
{
    for (int i = 0; i < m_Nrow; ++i)
    {
        uint* eliee = m_elieeList[i];
        for (int j = m_Ncol; j >= 0; --j)
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
        while (ss >> temp)
        {
            while (!getBit(eliee, --index));
            if (index!=temp)
            {
                printf("\ncheck result: Wrong!\n");
                printf("In row %d, temp = %d, index = %d\n",i,temp,index);
                return;
            }
        }
    }
    ifs.close();
    printf("\ncheck result: Right!\n");
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
