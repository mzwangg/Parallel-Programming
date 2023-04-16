#pragma once
#include <cstring>
#include <iostream>
#include <random>
#include <algorithm>
#include <math.h>
#include <sys/time.h>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>


#ifdef __ARM_NEON
#include <arm_neon.h>

#elif __x86_64
#include <immintrin.h> //MMX
#include <xmmintrin.h> //SSE
#include <emmintrin.h> //SSE2
#include <pmmintrin.h> //SSE3
#include <tmmintrin.h> //SSSE3
#include <smmintrin.h> //SSE4.1
#include <nmmintrin.h> //SSSE4.2
#include <immintrin.h> //AVX、AVX2、AVX-512
#endif


using namespace std;


#define ull (unsigned long long)
#define uint unsigned int


#ifndef BASH_TIMING
#define CHECK
#define N 1024
#define MAX_PRINT_N 20
#define STEP 1
#endif

#ifndef DATA
#define DATA "1_130_22_8"
#endif


const float eps = 1e-9;


