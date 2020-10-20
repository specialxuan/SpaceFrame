#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <mpi.h>
#include <omp.h>
#include <time.h>

#pragma warning(disable:4996)
#define _USE_MATH_DEFINES

/*全局变量声明*/
int TNN;//节点总数
int NFIN;//固定节点数
int NFRN;//可动节点数
int NOR;//杆件数
int WD;//危险杆件判断数
double* sigma;//杆件应力
int* SF;//杆件校核结果
double Adstress;//材料许用应力
double* XCN;//节点x方向的坐标
double* YCN;//节点y方向的坐标
int* BNR;//杆件始端节点号
int* ENR;//杆件末端节点号
double* TSR;//杆件的抗拉刚度
int* NWL;//具有载荷的节点号
double* XCL;//节点载荷x方向的分量
double* YCL;//节点载荷y方向的分量
double** LCS;//杆件的长度，倾角余弦，倾角正弦
double* DON;//节点位移分量
double* IFR;//杆件内力
double* ROU;//杆件密度
double* AREA;//杆件横截面积
double* I;//惯性矩
double* pp;//pp存放载荷向量
double** kk;//总刚度阵
double* Fcr;//压杆临界力

/*子函数声明*/
//从文本文件中读取数据
void PHRead();

//组集总刚度阵
double** PHBuildTotalStif();
//求杆件的长度，倾角余弦和倾角正弦
void PHLCosSin(int k);
//求单根杆件的单元刚度阵
void PHBuildUnitStif(int k, double us[2][2]);
//节点对号函数
int* PHI0J0(int k);
//组集载荷向量
double* PHBuildLoadVector();
//改进的平方根法求解方程组
void PHCholesky(double** a, double* b, double* x, int n);
//求解杆件内力
void PHRodForce();
//格式输出函数
void PHaaa();
//结果输出函数 
void PHPrint();
//内存释放
void PHFree();
//杆件自重函数
void PHselfweight();
//压杆稳定性校核
void PHstability();

int main()
{
	time_t start, end;


	//printf("%d\n", omp_get_num_procs());
	int k = 0;
	char c = 0;
	//存放用户输入的字符
	char value = 0;
	start = clock();
	printf("欢迎使用平面桁架结构求解器！\n");
	printf("\n\n*****\n");
	PHRead();//从数据文件读入数据
	kk = PHBuildTotalStif();//组集总刚度阵并将其指针赋给kk
	//printf("是否考虑杆件自重？(Y/N):");
	////用户根据信息提示选择输入字符
	//scanf("%c", &c);
	//if (c == 'Y' || c == 'y')
	//{
	//	PHselfweight();//执行自重函数
	//}

	pp = PHBuildLoadVector();//组集载荷向量并将其指针赋给pp
	PHCholesky(kk, pp, DON, 2 * NFRN);//改进的平方根法求节点位移，结果存放在DON中
	PHRodForce();//计算杆件内力
	PHPrint();//结构参数以及计算结果的输出	
	end = clock();
	printf("time=%f\n", ((double)end - start) / CLK_TCK);
}

void PHRead()
{
	FILE* fp;//定义文件指针
	char c = 0;//存放临时的字符型数据
	int i, j;//循环控制变量
	fp = fopen("structure_data.txt", "r");//为读取数据打开文件
	fseek(fp, 22L, 0);//将fp所指位置从初始位置向后移动22个字节
	fscanf(fp, "%d", &TNN);//读取fp所指位置的整形数据，存放在TNN中
	fseek(fp, 23L, 1);//将fp所指位置从当前位置向后移动23个字节
	fscanf(fp, "%d", &NFIN);//读取fp所指位置的整形数据，存放在NFIN中
	fseek(fp, 17L, 1);//将fp所指位置从当前位置向后移动17个字节
	fscanf(fp, "%d", &NOR);//读取fp所指位置的整形数据，存放在NOR中
	fseek(fp, 9L, 1);//将fp所指位置从当前位置向后移动9个字节
	fscanf(fp, "%lf", &Adstress);//读取fp所指位置的整形数据，存放在Adstress中
	fseek(fp, 2L, 1);//将fp所指位置从当前位置向后移动2个字节
	NFRN = TNN - NFIN;//计算可动节点数
	XCN = (double*)calloc(TNN, sizeof(double));//为XCN分配TNN个长度等于double变量的内存空间，下同
	memset(XCN, 0, TNN * sizeof(double));//内存空间初始化，下同
	YCN = (double*)calloc(TNN, sizeof(double));
	memset(YCN, 0, TNN * sizeof(double));
	BNR = (int*)calloc(NOR, sizeof(int));
	memset(BNR, 0, NOR * sizeof(int));
	ENR = (int*)calloc(NOR, sizeof(int));
	memset(ENR, 0, NOR * sizeof(int));
	TSR = (double*)calloc(NOR, sizeof(double));
	memset(TSR, 0, NOR * sizeof(double));
	NWL = (int*)calloc(TNN, sizeof(int));
	memset(NWL, 0, TNN * sizeof(int));
	XCL = (double*)calloc(TNN, sizeof(double));
	memset(XCL, 0, TNN * sizeof(double));
	YCL = (double*)calloc(TNN, sizeof(double));
	memset(YCL, 0, TNN * sizeof(double));
	DON = (double*)calloc(2 * NFRN, sizeof(double));
	memset(DON, 0, 2 * NFRN * sizeof(double));
	ROU = (double*)calloc(NOR, sizeof(double));
	memset(ROU, 0, NOR * sizeof(double));
	AREA = (double*)calloc(NOR, sizeof(double));
	memset(AREA, 0, NOR * sizeof(double));
	I = (double*)calloc(NOR, sizeof(double));
	memset(I, 0, NOR * sizeof(double));

	for (i = 0; i < 11; i++)//分别读取8组数据放在11个变量中
	{
		fseek(fp, 4L, 1);//将fp所指位置从当前位置向后移动4个字节
		j = 0;
		do
		{
			switch (i)//用switch语句控制对每组数据的读取
			{
			case 0:fscanf(fp, "%lf", &XCN[j]); break;
			case 1:fscanf(fp, "%lf", &YCN[j]); break;
			case 2:fscanf(fp, "%d", &BNR[j]); break;
			case 3:fscanf(fp, "%d", &ENR[j]); break;
			case 4:fscanf(fp, "%lf", &TSR[j]); break;
			case 5:fscanf(fp, "%d", &NWL[j]); break;
			case 6:fscanf(fp, "%lf", &XCL[j]); break;
			case 7:fscanf(fp, "%lf", &YCL[j]); break;
			case 8:fscanf(fp, "%lf", &ROU[j]); break;
			case 9:fscanf(fp, "%lf", &AREA[j]); break;
			case 10:fscanf(fp, "%lf", &I[j]); break;
			}
			fscanf(fp, "%c", &c);//读取每个数据后的逗号或换行符
			j++;//数据指标自加
		} while (c != '\n');//若读取的数据后面不是换行符则继续读取
	}
}
double** PHBuildTotalStif()
{
	double** kk, us[2][2];//kk为总刚度阵，us为单元刚度阵
	int i, j, k, m, n, * p;//i,j,m,n为循环控制变量，k为杆件程序编号，p为存放杆端节点对号位置的数组的指针
	kk = (double**)calloc(2 * NFRN, sizeof(double*));//以下三行语句为kk申请二维存储空间
	for (i = 0; i < 2 * NFRN; i++)
		*(kk + i) = (double*)calloc(2 * NFRN, sizeof(double));
	for (i = 0; i < 2 * NFRN; i++)//以下三行语句对kk指向的总刚度阵清零
		for (j = 0; j < 2 * NFRN; j++)
			kk[i][j] = 0;
	LCS = (double**)calloc(3, sizeof(double*));//以下三行语句对存放杆件几何参数的LCS申请二维存储空间
	for (i = 0; i < 3; i++)
		*(LCS + i) = (double*)calloc(NOR, sizeof(double));
	for (k = 0; k < NOR; k++)//kk为杆件程序编号，对每一根杆件循环，组装总刚度阵
	{
		PHLCosSin(k + 1);//计算杆件的几何参数
		PHBuildUnitStif(k, us);//计算程序编号为k的杆件的单元刚度阵分块，并存放在us中
		p = PHI0J0(k + 1);//计算程序编号为k的杆件端点对号位置，并存放在p指向的数组中
		for (i = 0; i < 2; i++)
		{
			if (p[i] >= 0)//符合条件说明为可动节点并叠加，否则不叠加
			{
				for (m = 0; m < 2; m++)
					for (n = 0; n < 2; n++)
						kk[p[i] + m][p[i] + n] += us[m][n];//对us的四个元素按相应位置进行叠加
			}
		}
		if (p[0] >= 0 && p[1] >= 0)//符合条件说明两端点均为可动节点并进行叠加，否则不叠加
		{
			for (i = 0; i < 2; i++)
			{
				for (m = 0; m < 2; m++)
					for (n = 0; n < 2; n++)
						kk[p[i] + m][p[1 - i] + n] -= us[m][n];//对us中的四个元素的相反数按相应位置进行叠加
			}
		}
	}

	return kk;
}
void PHLCosSin(int k)//k为杆件实际编号
{
	int i, j;
	k--;//k自减，即为杆件程序编号
	i = BNR[k] - 1;//i存放杆件的始端节点对应程序中的数组指标
	j = ENR[k] - 1;//j存放杆件的末端节点对应程序中的数组指标
	LCS[1][k] = XCN[j] - XCN[i];//杆件始末端节点横坐标之差
	LCS[2][k] = YCN[j] - YCN[i];//杆件始末端节点纵坐标之差
	LCS[0][k] = sqrt(LCS[1][k] * LCS[1][k] + LCS[2][k] * LCS[2][k]);//求杆件长度
	LCS[1][k] /= LCS[0][k];//求杆件余弦值
	LCS[2][k] /= LCS[0][k];//求杆件正弦值
}
void PHBuildUnitStif(int k, double us[2][2])//k为杆件程序编号，us为单元刚度矩阵
{
	int i, j;//i,j为循环控制变量
	double rd;//rd存放抗拉刚度系数
	rd = TSR[k] / LCS[0][k];//计算抗拉刚度系数
	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
			us[i][j] = LCS[i + 1][k] * LCS[j + 1][k] * rd;//计算us中各元素值并赋值
}
int* PHI0J0(int k)
{
	int bl, br, ij[2];
	bl = BNR[k - 1];//bl存放杆件的始端节点号
	br = ENR[k - 1];//br存放杆件的末端节点号
	ij[0] = 2 * (bl - NFIN - 1);//将始端节点在总刚度阵中对应的位置编号存放在ij数组中
	ij[1] = 2 * (br - NFIN - 1);//将末端节点在总刚度阵中对应的位置编号存放在ij数组中
	return ij;//返回ij数组指针
}
void PHselfweight()//自重函数
{
	int bl, br, k, * p;
	double w = 0, g = 9.81;
	for (k = 0; k < NOR; k++)
	{
		bl = BNR[k];//bl存放杆件的始端节点号
		br = ENR[k];//br存放杆件的末端节点号
		PHLCosSin(k + 1);//调用函数，求解杆长
		w = ROU[k] * AREA[k] * LCS[0][k] * g;//计算杆件自重
		if ((bl - NFRN - 1) >= 0)//左端点非固定
		{
			YCL[bl - 1] -= 0.5 * w / 1000;//在载荷向量中叠加自重
		}
		if ((br - NFRN - 1) >= 0)//右端点非固定
		{
			YCL[br - 1] -= 0.5 * w / 1000;//在载荷向量中叠加自重
		}
	}
}
double* PHBuildLoadVector()
{
	int i;//i为循环变量
	pp = (double*)calloc(2 * NFRN, sizeof(double));//为pp分配2*NFRN个长度等于double的内存空间
	for (i = 0; i < 2 * NFRN; i++)//载荷变量清零
	{
		pp[i] = 0;
	}
	i = 0;//循环变量清零
	for (i = NFIN; i < TNN; i++)
	{
		pp[2 * (i - NFIN)] = XCL[i];
		pp[2 * (i - NFIN) + 1] = YCL[i];
	}
	return pp;//返回载荷变量数组指针
}
void PHCholesky(double** A, double* b, double* x, int n)
//A为对称系数阵，b为常数向量，x为未知数向量，n为维数
{
	int i, j, k;//循环控制变量
	double s, ** L, * D, * y;//s为中间变量，L,D为分解矩阵，y为中间向量
	L = (double**)calloc(n, sizeof(double*));//为L申请n个长度为double的内存空间
	for (i = 0; i < n; i++)
		*(L + i) = (double*)calloc(n, sizeof(double));
	D = (double*)calloc(n, sizeof(double));//为D申请n个长度为double 的内存空间
	y = (double*)calloc(n, sizeof(double));//为y申请n个长度为double的内存空间
	for (i = 0; i < n; i++)
		L[i][i] = 1;//L初始化
	/*将A分解为LDL(t)*/
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			s = 0;
			for (k = 0; k < j; k++)
			{
				s = s + L[i][k] * D[k] * L[j][k];
			}
			L[i][j] = (A[i][j] - s) / D[j];
		}
		s = 0;
		for (k = 0; k < i; k++)
			s = s + L[i][k] * L[i][k] * D[k];
		D[i] = A[i][i] - s;
	}
	/*由Ly=b求解y*/
	y[0] = b[0];
	for (i = 1; i < n; i++)
	{
		s = 0;
		for (k = 0; k < i; k++)
		{
			s = s + L[i][k] * y[k];
		}
		y[i] = b[i] - s;
	}
	/*由DL(T)x=y求解x*/
	x[n - 1] = y[n - 1] / D[n - 1];
	for (i = n - 2; i >= 0; i--)
	{
		s = 0;
		for (k = i + 1; k < n; k++)
		{
			s = s + L[k][i] * x[k];
		}
		x[i] = y[i] / D[i] - s;
	}
}
void PHRodForce()//计算杆件内力
{
	int i, j, k, * p, WD = 0;//i,j为普通变量，k为杆件程序号
	//p为存放杆端节点对号位置的数组的指针；
	double d1[2], d2[2], rd;//d1,d2分别为杆件始末端的位移分量，rd为抗拉刚度系数
	IFR = (double*)calloc(NOR, sizeof(double));//为IFR分配NOR个长度等于double变量的内存空间
	sigma = (double*)calloc(NOR, sizeof(double));
	SF = (int*)calloc(NOR, sizeof(int));
	for (k = 0; k < NOR; k++)//对所有杆件进行循环
	{
		p = PHI0J0(k + 1);//计算程序编号为k的杆件端点对号位置，并存放在p指向的数组中
		i = p[0];
		j = p[1];//对i,j进行赋值
		if (i < 0)//i<0则节点为固定节点
		{
			d1[0] = d1[1] = 0;//固定节点位移分量为0
		}
		else
		{
			d1[0] = DON[i];
			d1[1] = DON[i + 1];//始端节点为可动节点时将DON中数据赋给d1
		}
		d2[0] = DON[j];
		d2[1] = DON[j + 1];//将DON中末端节点位移数据赋给d2(根据编号规则末端节点必为可动节点）
		rd = TSR[k] / LCS[0][k];//计算杆件的抗拉刚度
		IFR[k] = rd * (LCS[1][k] * (d2[0] - d1[0]) + LCS[2][k] * (d2[1] - d1[1]));//计算杆件内力
		sigma[k] = IFR[k] / AREA[k];//计算杆件应力
		Fcr = (double*)calloc(NOR, sizeof(double));//为Fcr调用NOR个长度为double的内存空间
		for (i = 0; i < NOR; i++)
		{

			PHLCosSin(i + 1);
			Fcr[i] = 3.1415 * 3.1415 * TSR[i] / AREA[i] * I[i] / (LCS[0][i] * LCS[0][i]) / 1000.0;//计算压杆临界力，单位为kN
		}
		for (i = 0; i < NOR; i++)
		{
			if (IFR[i] < 0)
			{

				if (fabs(IFR[i]) > Fcr[i] || fabs(sigma[i]) > Adstress)//压杆稳定性与应力校核
				{
					SF[i] = i + 1;
				}
			}
			else
			{
				if (fabs(sigma[i]) > Adstress)//拉杆内力校核
				{
					SF[i] = i + 1;
				}
			}
		}
	}
}
void PHFree()
{
	free(LCS);//以下均为申请内存空间的内存释放
	free(XCN);
	free(YCN);
	free(BNR);
	free(ENR);
	free(TSR);
	free(NWL);
	free(XCL);
	free(YCL);
	free(DON);
	free(I);
	free(sigma);
	free(Fcr);
	free(ROU);
	free(AREA);
	free(IFR);
	free(SF);
}
void PHaaa()
{
	printf("*********************************\n");
}
void PHPrint()
{
	int i, j, * p = NWL, n = 0;//指针p指向NWL首地址
	printf("\t\t\t\t平面桁架结构计算\n");
	PHaaa();
	printf("\t节点总数=%d\t固定节点总数=%d\t可动节点总数=%d\t杆件数=%d\n", TNN, NFIN, NFRN, NOR);
	PHaaa();
	printf("\t节点号\t\tx坐标\ty坐标\n");
	for (i = 1; i <= TNN; i++)
		printf("\t%d\t\t%5.4f\t\t%5.4f\n", i, XCN[i - 1], YCN[i - 1]);
	PHaaa();
	printf("\t杆件号\t\t始端节点号\t末端节点号\t抗拉刚度EA\n");
	for (i = 1; i <= NOR; i++)
		printf("\t%d\t\t%d\t\t%d\t\t%5.4f\n", i, BNR[i - 1], ENR[i - 1], TSR[i - 1]);
	PHaaa();
	printf("\t荷载节点号\tx分量\t\ty分量\n");

	for (i = NFIN; i < TNN; i++)
	{
		printf("\t%d\t\t%5.4f\t\t%5.4f\n", i + 1, pp[2 * (i - NFIN)], pp[2 * (i - NFIN) + 1]);
	}
	PHaaa();
	printf("计算结果\t结点位移输出:\n");
	PHaaa();
	printf("\t节点号\t\t位移x分量\t位移y分量\n");
	for (i = 0; i < NFRN; i++)
	{
		printf("\t%d\t\t%5.7f\t%5.7f\n", NFIN + i + 1, DON[2 * i], DON[2 * i + 1]);
	}
	PHaaa();
	printf("计算结果\t杆件内力输出：\n");
	PHaaa();
	printf("\t杆件号\t\t截面系数EA\t\t\t杆长l\t\t杆件内力\t\t杆件应力\n");
	for (i = 1; i <= NOR; i++)
	{
		printf("\t%d\t\t%5.4f\t\t\t%5.4f\t\t%5.7f\t\t%5.7f\t\n", i, TSR[i - 1], LCS[0][i - 1], IFR[i - 1], sigma[i - 1]);
	}
	PHaaa();

	printf("计算结果\t杆件压杆稳定及应力校核结果:\n");
	for (i = 1; i <= NOR; i++)
	{
		if (SF[i - 1] != 0)
		{
			n++;
			printf("\t\t%d\t\n", SF[i - 1]);
		}
	}
	if (n == 0)
	{
		printf("\n\t\t所有杆件均安全\n");
	}
	else
		printf("\t\t以上为危险杆件号！\n");

	PHaaa();

	printf("\n\t\t\t\t感谢您的使用!\n\n\n\n\n\n");
}