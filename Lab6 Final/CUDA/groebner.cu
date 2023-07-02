#include "EliminateClass.h"

int Nmat;//矩阵阶数
int Nee;//被消元行行数
int Ner;//消元子行数
int Ncol;//转为uint后对应的列数

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

   eliminater = new Eliminater(Nmat);
   eliminatee = new Eliminatee(Nmat, Nee);
   globalEliminater = new Eliminater(DATAPATH + DATA + Para1, Nmat, Ner);
   globalEliminatee = new Eliminatee(DATAPATH + DATA + Para2, Nmat, Nee);
}

void timing(double(*func)(int), int dataIndex, int blockSize = 1024)
{
   string DATA = dataArr[dataIndex - 1];
   dataInit(DATA);

   double totalTime = 0.0f;
   //先进行WARMUP次热身，不计入耗时，使得测量时间更加精确
   for (int run = 0; run < REPEAT_NUM + WARMUP; run++) {
       eliminater->copy(*globalEliminater);
       eliminatee->copy(*globalEliminatee);
       float duration = func(blockSize);
       if (run >= WARMUP) totalTime += duration;
   }
   cout << totalTime / REPEAT_NUM << '\t';

   //检验答案
   eliminatee->check(DATAPATH + DATA + Para3);

   //释放内存
   delete eliminatee;
   delete eliminater;
   delete globalEliminater;
   delete globalEliminatee;
}

void timingAll(string nameArr[], double(*funcArr[])(int), int funcNum, int start, int end, int blockSize = 1024)
{
   //输出表头
   cout << "Problem\t";
   for (int i = 0; i < funcNum; ++i)
       cout << nameArr[i] << '\t';
   cout << endl;

   for (int dataIndex = start; dataIndex <= end; ++dataIndex) {
       cout << dataIndex << '\t';
       string DATA = dataArr[dataIndex - 1];
       dataInit(DATA);

       for (int funcIndex = 0; funcIndex < funcNum; ++funcIndex) {
           double totalTime = 0.0f;
           //先进行WARMUP次热身，不计入耗时，使得测量时间更加精确
           for (int run = 0; run < REPEAT_NUM + WARMUP; run++) {
               eliminater->copy(*globalEliminater);
               eliminatee->copy(*globalEliminatee);
               float duration = funcArr[funcIndex](blockSize);
               if (run >= WARMUP) totalTime += duration;
           }
           cout << totalTime / REPEAT_NUM << '\t';
           //检验答案
           eliminatee->check(DATAPATH + DATA + Para3);
       }
       cout << endl;

       delete eliminatee;
       delete eliminater;
       delete globalEliminater;
       delete globalEliminatee;
   }
}

double serial1(int blockSize = 0)
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

double serial2(int blockSize = 0)
{
   auto start = chrono::high_resolution_clock::now();

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

   auto stop = chrono::high_resolution_clock::now();
   auto duration = chrono::duration_cast<chrono::microseconds>(stop - start) / 1000000.0;
   return duration.count();
}


double serial3(int blockSize = 0)
{
   auto start = chrono::high_resolution_clock::now();

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

   auto stop = chrono::high_resolution_clock::now();
   auto duration = chrono::duration_cast<chrono::microseconds>(stop - start) / 1000000.0;
   return duration.count();
}

__global__ void eliminate_kernel(uint* elieeData, uint* eliers, bool* isElierList, int j, int endj, int Nee, int Ncol)
{
   uint* elier, * eliee;
   for (int jj = j; jj > endj; --jj) {
       elier = eliers + (jj - endj - 1) * Ncol;
       for (int i = blockIdx.x; i < Nee; i += gridDim.x) {
           if (isElierList[i])
               continue;
           eliee = elieeData + i * Ncol;
           if (eliee[jj >> 5] & (0x80000000 >> (jj & 0x1f)))
           { // 如果当前行需要消元
               for (int k = threadIdx.x; k < Ncol; k += blockDim.x)//通过循环划分划分数据
                   eliee[k] ^= elier[k];
           }
       }
       __syncthreads();//块内同步
   }
}

double cuda(int blockSize) {
   //开始计时
   float elapsedTime;
   cudaEvent_t start, stop;
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   cudaEventRecord(start, 0);

   //得到多处理器数量，从而确定最佳网格维度
   int deviceId;
   int numberOfSMs;
   cudaGetDevice(&deviceId);
   cudaDeviceGetAttribute(&numberOfSMs, cudaDevAttrMultiProcessorCount, deviceId);

   // 申请GPU内存
   bool* isElierList;
   uint* elieeData, * eliers;
   cudaMalloc((void**)&isElierList, Nee);
   cudaMalloc((void**)&elieeData, Nee * Ncol * sizeof(uint));
   cudaMalloc((void**)&eliers, groupSize * Ncol * sizeof(uint));

   //将数据从主机内存复制到GPU内存
   cudaMemcpy(elieeData, eliminatee->dataVector, Nee * Ncol * sizeof(uint), cudaMemcpyHostToDevice);

   //执行核函数
   int endj = Nmat - 1;
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
                   memcpy(eliminater->m_elierList[j], eliee, Ncol * sizeof(uint));
                   eliminater->m_isElierList[j] = true;
                   eliminatee->m_isElierList[i] = true;
                   break;
               }
           }
           if (i == Nee) --endj;
           continue;
       }

       //将消元子和是否升格的信息传递到GPU中
       cudaMemcpy(eliers, eliminater->dataVector + (endj + 1) * Ncol, (j - endj) * Ncol * sizeof(uint), cudaMemcpyHostToDevice);
       cudaMemcpy(isElierList, eliminatee->m_isElierList, Ncol , cudaMemcpyHostToDevice);

       //进行消元
       eliminate_kernel << <numberOfSMs, blockSize >> > (elieeData, eliers, isElierList, j, endj, Nee, Ncol);

       //等待所有进程消元完成
       cudaDeviceSynchronize();


       if (!eliminater->m_isElierList[endj]) {//判断是否将被消元行复制回CPU中
           cudaMemcpy(eliminatee->dataVector, elieeData, Nee * Ncol * sizeof(uint), cudaMemcpyDeviceToHost);
       }
   }


   // 将数据从GPU内存复制到主机内存
   cudaMemcpy(eliminatee->dataVector, elieeData, Nee * Ncol * sizeof(uint), cudaMemcpyDeviceToHost);

   //结束计时
   cudaEventRecord(stop, 0);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&elapsedTime, start, stop);

   //释放内存
   cudaFree(isElierList);
   cudaFree(elieeData);
   cudaFree(eliers);

   return elapsedTime / 1000.0;
}

int main()
{
   //timing(serial1, 2, 128);

   string namearr[] = { "serial1","serial3","cuda"};
   double(*funcArr[])(int) = { serial1,serial3,cuda };
   timingAll(namearr, funcArr, 3, 11, 11);
}