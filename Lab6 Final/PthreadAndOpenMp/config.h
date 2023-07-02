#pragma once
#include <cstring>
#include <iostream>
#include <random>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <pthread.h>
#include "systime.h"
using namespace std;
#pragma comment(lib, "pthreadVC2.lib")


#define ull (unsigned long long)
#define uint unsigned int


#ifndef __ARM_NEON
#include <immintrin.h> //MMX
#include <xmmintrin.h> //SSE
#include <emmintrin.h> //SSE2
#include <pmmintrin.h> //SSE3
#include <tmmintrin.h> //SSSE3
#include <smmintrin.h> //SSE4.1
#include <nmmintrin.h> //SSSE4.2
#include <immintrin.h> //AVX+AVX2AVX-512
#else
#include <arm_neon.h>
#endif


#ifdef _WIN32
const string DATAPATH = ".\\data\\Groebner\\";
const string Para1 = "\\消元子.txt";
const string Para2 = "\\被消元行.txt";
const string Para3 = "\\消元结果.txt";
#elif __ARM_NEON
const string DATAPATH = "/home/data/Groebner/";
const string Para1 = "/1.txt";
const string Para2 = "/2.txt";
const string Para3 = "/3,txt";
#else
const string DATAPATH = "/home/data/Groebner/";
const string Para1 = "/消元子.txt";
const string Para2 = "/被消元行.txt";
const string Para3 = "/消元结果,txt";
#endif


#ifdef TIMING
#define MAXN 2048
#define REPEAT_NUM 5
#define MAX_PRINT_N 0
#else
#define CHECK
#define MAXN 256
#define REPEAT_NUM 2
#define MAX_PRINT_N 20
#endif

struct gaussPara
{
    int n;
    int k;
    int t_id;
    int threadNum;
    float* mat;
};

struct groebnerPara
{
    void* eliminater;
    void* eliminatee;
    int t_id;
    int threadNum;
};

const float eps = 1e-9;
const int floatLimit = (1 << 24);
const string dataArr[] = { "1_130_22_8","2_254_106_53","3_562_170_53","4_1011_539_263","5_2362_1226_453","6_3799_2759_1953",\
                        "7_8399_6375_4535","8_23045_18748_14325","9_37960_29304_14921","10_43577_39477_54274","11_85401_5724_756" };
const int threadNumArr[] = { 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,20,24,28,32 };
const int dataSize = 11;