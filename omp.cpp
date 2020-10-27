#include <iostream>
#include <omp.h>
#include <stdlib.h>
int main()
{
    int sum = 0;
    int a[10] = {1,2,3,4,5,6,7,8,9,10};
    int coreNum = omp_get_num_procs();//获得处理器个数
    int* sumArray = new int[coreNum];//对应处理器个数，先生成一个数组
    for (int i=0;i<coreNum;i++)//将数组各元素初始化为0
        sumArray[i] = 0;
#pragma omp parallel for
    for (int i=0;i<10;i++)
    {
        int k = omp_get_thread_num();//获得每个线程的ID
        sumArray[k] = sumArray[k]+a[i];
        printf("\t%d\n", k);
    }
    for (int i = 0;i<coreNum;i++)
        sum = sum + sumArray[i];
   std::cout<<"sum: "<<sum<<std::endl;
   getchar();
   return 0;
}
