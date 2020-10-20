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

/*ȫ�ֱ�������*/
int TNN;//�ڵ�����
int NFIN;//�̶��ڵ���
int NFRN;//�ɶ��ڵ���
int NOR;//�˼���
int WD;//Σ�ո˼��ж���
double* sigma;//�˼�Ӧ��
int* SF;//�˼�У�˽��
double Adstress;//��������Ӧ��
double* XCN;//�ڵ�x���������
double* YCN;//�ڵ�y���������
int* BNR;//�˼�ʼ�˽ڵ��
int* ENR;//�˼�ĩ�˽ڵ��
double* TSR;//�˼��Ŀ����ն�
int* NWL;//�����غɵĽڵ��
double* XCL;//�ڵ��غ�x����ķ���
double* YCL;//�ڵ��غ�y����ķ���
double** LCS;//�˼��ĳ��ȣ�������ң��������
double* DON;//�ڵ�λ�Ʒ���
double* IFR;//�˼�����
double* ROU;//�˼��ܶ�
double* AREA;//�˼�������
double* I;//���Ծ�
double* pp;//pp����غ�����
double** kk;//�ܸն���
double* Fcr;//ѹ���ٽ���

/*�Ӻ�������*/
//���ı��ļ��ж�ȡ����
void PHRead();

//�鼯�ܸն���
double** PHBuildTotalStif();
//��˼��ĳ��ȣ�������Һ��������
void PHLCosSin(int k);
//�󵥸��˼��ĵ�Ԫ�ն���
void PHBuildUnitStif(int k, double us[2][2]);
//�ڵ�Ժź���
int* PHI0J0(int k);
//�鼯�غ�����
double* PHBuildLoadVector();
//�Ľ���ƽ��������ⷽ����
void PHCholesky(double** a, double* b, double* x, int n);
//���˼�����
void PHRodForce();
//��ʽ�������
void PHaaa();
//���������� 
void PHPrint();
//�ڴ��ͷ�
void PHFree();
//�˼����غ���
void PHselfweight();
//ѹ���ȶ���У��
void PHstability();

int main()
{
	time_t start, end;


	//printf("%d\n", omp_get_num_procs());
	int k = 0;
	char c = 0;
	//����û�������ַ�
	char value = 0;
	start = clock();
	printf("��ӭʹ��ƽ����ܽṹ�������\n");
	printf("\n\n*****\n");
	PHRead();//�������ļ���������
	kk = PHBuildTotalStif();//�鼯�ܸն��󲢽���ָ�븳��kk
	//printf("�Ƿ��Ǹ˼����أ�(Y/N):");
	////�û�������Ϣ��ʾѡ�������ַ�
	//scanf("%c", &c);
	//if (c == 'Y' || c == 'y')
	//{
	//	PHselfweight();//ִ�����غ���
	//}

	pp = PHBuildLoadVector();//�鼯�غ�����������ָ�븳��pp
	PHCholesky(kk, pp, DON, 2 * NFRN);//�Ľ���ƽ��������ڵ�λ�ƣ���������DON��
	PHRodForce();//����˼�����
	PHPrint();//�ṹ�����Լ������������	
	end = clock();
	printf("time=%f\n", ((double)end - start) / CLK_TCK);
}

void PHRead()
{
	FILE* fp;//�����ļ�ָ��
	char c = 0;//�����ʱ���ַ�������
	int i, j;//ѭ�����Ʊ���
	fp = fopen("structure_data.txt", "r");//Ϊ��ȡ���ݴ��ļ�
	fseek(fp, 22L, 0);//��fp��ָλ�ôӳ�ʼλ������ƶ�22���ֽ�
	fscanf(fp, "%d", &TNN);//��ȡfp��ָλ�õ��������ݣ������TNN��
	fseek(fp, 23L, 1);//��fp��ָλ�ôӵ�ǰλ������ƶ�23���ֽ�
	fscanf(fp, "%d", &NFIN);//��ȡfp��ָλ�õ��������ݣ������NFIN��
	fseek(fp, 17L, 1);//��fp��ָλ�ôӵ�ǰλ������ƶ�17���ֽ�
	fscanf(fp, "%d", &NOR);//��ȡfp��ָλ�õ��������ݣ������NOR��
	fseek(fp, 9L, 1);//��fp��ָλ�ôӵ�ǰλ������ƶ�9���ֽ�
	fscanf(fp, "%lf", &Adstress);//��ȡfp��ָλ�õ��������ݣ������Adstress��
	fseek(fp, 2L, 1);//��fp��ָλ�ôӵ�ǰλ������ƶ�2���ֽ�
	NFRN = TNN - NFIN;//����ɶ��ڵ���
	XCN = (double*)calloc(TNN, sizeof(double));//ΪXCN����TNN�����ȵ���double�������ڴ�ռ䣬��ͬ
	memset(XCN, 0, TNN * sizeof(double));//�ڴ�ռ��ʼ������ͬ
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

	for (i = 0; i < 11; i++)//�ֱ��ȡ8�����ݷ���11��������
	{
		fseek(fp, 4L, 1);//��fp��ָλ�ôӵ�ǰλ������ƶ�4���ֽ�
		j = 0;
		do
		{
			switch (i)//��switch�����ƶ�ÿ�����ݵĶ�ȡ
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
			fscanf(fp, "%c", &c);//��ȡÿ�����ݺ�Ķ��Ż��з�
			j++;//����ָ���Լ�
		} while (c != '\n');//����ȡ�����ݺ��治�ǻ��з��������ȡ
	}
}
double** PHBuildTotalStif()
{
	double** kk, us[2][2];//kkΪ�ܸն���usΪ��Ԫ�ն���
	int i, j, k, m, n, * p;//i,j,m,nΪѭ�����Ʊ�����kΪ�˼������ţ�pΪ��Ÿ˶˽ڵ�Ժ�λ�õ������ָ��
	kk = (double**)calloc(2 * NFRN, sizeof(double*));//�����������Ϊkk�����ά�洢�ռ�
	for (i = 0; i < 2 * NFRN; i++)
		*(kk + i) = (double*)calloc(2 * NFRN, sizeof(double));
	for (i = 0; i < 2 * NFRN; i++)//������������kkָ����ܸն�������
		for (j = 0; j < 2 * NFRN; j++)
			kk[i][j] = 0;
	LCS = (double**)calloc(3, sizeof(double*));//�����������Դ�Ÿ˼����β�����LCS�����ά�洢�ռ�
	for (i = 0; i < 3; i++)
		*(LCS + i) = (double*)calloc(NOR, sizeof(double));
	for (k = 0; k < NOR; k++)//kkΪ�˼������ţ���ÿһ���˼�ѭ������װ�ܸն���
	{
		PHLCosSin(k + 1);//����˼��ļ��β���
		PHBuildUnitStif(k, us);//���������Ϊk�ĸ˼��ĵ�Ԫ�ն���ֿ飬�������us��
		p = PHI0J0(k + 1);//���������Ϊk�ĸ˼��˵�Ժ�λ�ã��������pָ���������
		for (i = 0; i < 2; i++)
		{
			if (p[i] >= 0)//��������˵��Ϊ�ɶ��ڵ㲢���ӣ����򲻵���
			{
				for (m = 0; m < 2; m++)
					for (n = 0; n < 2; n++)
						kk[p[i] + m][p[i] + n] += us[m][n];//��us���ĸ�Ԫ�ذ���Ӧλ�ý��е���
			}
		}
		if (p[0] >= 0 && p[1] >= 0)//��������˵�����˵��Ϊ�ɶ��ڵ㲢���е��ӣ����򲻵���
		{
			for (i = 0; i < 2; i++)
			{
				for (m = 0; m < 2; m++)
					for (n = 0; n < 2; n++)
						kk[p[i] + m][p[1 - i] + n] -= us[m][n];//��us�е��ĸ�Ԫ�ص��෴������Ӧλ�ý��е���
			}
		}
	}

	return kk;
}
void PHLCosSin(int k)//kΪ�˼�ʵ�ʱ��
{
	int i, j;
	k--;//k�Լ�����Ϊ�˼�������
	i = BNR[k] - 1;//i��Ÿ˼���ʼ�˽ڵ��Ӧ�����е�����ָ��
	j = ENR[k] - 1;//j��Ÿ˼���ĩ�˽ڵ��Ӧ�����е�����ָ��
	LCS[1][k] = XCN[j] - XCN[i];//�˼�ʼĩ�˽ڵ������֮��
	LCS[2][k] = YCN[j] - YCN[i];//�˼�ʼĩ�˽ڵ�������֮��
	LCS[0][k] = sqrt(LCS[1][k] * LCS[1][k] + LCS[2][k] * LCS[2][k]);//��˼�����
	LCS[1][k] /= LCS[0][k];//��˼�����ֵ
	LCS[2][k] /= LCS[0][k];//��˼�����ֵ
}
void PHBuildUnitStif(int k, double us[2][2])//kΪ�˼������ţ�usΪ��Ԫ�նȾ���
{
	int i, j;//i,jΪѭ�����Ʊ���
	double rd;//rd��ſ����ն�ϵ��
	rd = TSR[k] / LCS[0][k];//���㿹���ն�ϵ��
	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
			us[i][j] = LCS[i + 1][k] * LCS[j + 1][k] * rd;//����us�и�Ԫ��ֵ����ֵ
}
int* PHI0J0(int k)
{
	int bl, br, ij[2];
	bl = BNR[k - 1];//bl��Ÿ˼���ʼ�˽ڵ��
	br = ENR[k - 1];//br��Ÿ˼���ĩ�˽ڵ��
	ij[0] = 2 * (bl - NFIN - 1);//��ʼ�˽ڵ����ܸն����ж�Ӧ��λ�ñ�Ŵ����ij������
	ij[1] = 2 * (br - NFIN - 1);//��ĩ�˽ڵ����ܸն����ж�Ӧ��λ�ñ�Ŵ����ij������
	return ij;//����ij����ָ��
}
void PHselfweight()//���غ���
{
	int bl, br, k, * p;
	double w = 0, g = 9.81;
	for (k = 0; k < NOR; k++)
	{
		bl = BNR[k];//bl��Ÿ˼���ʼ�˽ڵ��
		br = ENR[k];//br��Ÿ˼���ĩ�˽ڵ��
		PHLCosSin(k + 1);//���ú��������˳�
		w = ROU[k] * AREA[k] * LCS[0][k] * g;//����˼�����
		if ((bl - NFRN - 1) >= 0)//��˵�ǹ̶�
		{
			YCL[bl - 1] -= 0.5 * w / 1000;//���غ������е�������
		}
		if ((br - NFRN - 1) >= 0)//�Ҷ˵�ǹ̶�
		{
			YCL[br - 1] -= 0.5 * w / 1000;//���غ������е�������
		}
	}
}
double* PHBuildLoadVector()
{
	int i;//iΪѭ������
	pp = (double*)calloc(2 * NFRN, sizeof(double));//Ϊpp����2*NFRN�����ȵ���double���ڴ�ռ�
	for (i = 0; i < 2 * NFRN; i++)//�غɱ�������
	{
		pp[i] = 0;
	}
	i = 0;//ѭ����������
	for (i = NFIN; i < TNN; i++)
	{
		pp[2 * (i - NFIN)] = XCL[i];
		pp[2 * (i - NFIN) + 1] = YCL[i];
	}
	return pp;//�����غɱ�������ָ��
}
void PHCholesky(double** A, double* b, double* x, int n)
//AΪ�Գ�ϵ����bΪ����������xΪδ֪��������nΪά��
{
	int i, j, k;//ѭ�����Ʊ���
	double s, ** L, * D, * y;//sΪ�м������L,DΪ�ֽ����yΪ�м�����
	L = (double**)calloc(n, sizeof(double*));//ΪL����n������Ϊdouble���ڴ�ռ�
	for (i = 0; i < n; i++)
		*(L + i) = (double*)calloc(n, sizeof(double));
	D = (double*)calloc(n, sizeof(double));//ΪD����n������Ϊdouble ���ڴ�ռ�
	y = (double*)calloc(n, sizeof(double));//Ϊy����n������Ϊdouble���ڴ�ռ�
	for (i = 0; i < n; i++)
		L[i][i] = 1;//L��ʼ��
	/*��A�ֽ�ΪLDL(t)*/
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
	/*��Ly=b���y*/
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
	/*��DL(T)x=y���x*/
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
void PHRodForce()//����˼�����
{
	int i, j, k, * p, WD = 0;//i,jΪ��ͨ������kΪ�˼������
	//pΪ��Ÿ˶˽ڵ�Ժ�λ�õ������ָ�룻
	double d1[2], d2[2], rd;//d1,d2�ֱ�Ϊ�˼�ʼĩ�˵�λ�Ʒ�����rdΪ�����ն�ϵ��
	IFR = (double*)calloc(NOR, sizeof(double));//ΪIFR����NOR�����ȵ���double�������ڴ�ռ�
	sigma = (double*)calloc(NOR, sizeof(double));
	SF = (int*)calloc(NOR, sizeof(int));
	for (k = 0; k < NOR; k++)//�����и˼�����ѭ��
	{
		p = PHI0J0(k + 1);//���������Ϊk�ĸ˼��˵�Ժ�λ�ã��������pָ���������
		i = p[0];
		j = p[1];//��i,j���и�ֵ
		if (i < 0)//i<0��ڵ�Ϊ�̶��ڵ�
		{
			d1[0] = d1[1] = 0;//�̶��ڵ�λ�Ʒ���Ϊ0
		}
		else
		{
			d1[0] = DON[i];
			d1[1] = DON[i + 1];//ʼ�˽ڵ�Ϊ�ɶ��ڵ�ʱ��DON�����ݸ���d1
		}
		d2[0] = DON[j];
		d2[1] = DON[j + 1];//��DON��ĩ�˽ڵ�λ�����ݸ���d2(���ݱ�Ź���ĩ�˽ڵ��Ϊ�ɶ��ڵ㣩
		rd = TSR[k] / LCS[0][k];//����˼��Ŀ����ն�
		IFR[k] = rd * (LCS[1][k] * (d2[0] - d1[0]) + LCS[2][k] * (d2[1] - d1[1]));//����˼�����
		sigma[k] = IFR[k] / AREA[k];//����˼�Ӧ��
		Fcr = (double*)calloc(NOR, sizeof(double));//ΪFcr����NOR������Ϊdouble���ڴ�ռ�
		for (i = 0; i < NOR; i++)
		{

			PHLCosSin(i + 1);
			Fcr[i] = 3.1415 * 3.1415 * TSR[i] / AREA[i] * I[i] / (LCS[0][i] * LCS[0][i]) / 1000.0;//����ѹ���ٽ�������λΪkN
		}
		for (i = 0; i < NOR; i++)
		{
			if (IFR[i] < 0)
			{

				if (fabs(IFR[i]) > Fcr[i] || fabs(sigma[i]) > Adstress)//ѹ���ȶ�����Ӧ��У��
				{
					SF[i] = i + 1;
				}
			}
			else
			{
				if (fabs(sigma[i]) > Adstress)//��������У��
				{
					SF[i] = i + 1;
				}
			}
		}
	}
}
void PHFree()
{
	free(LCS);//���¾�Ϊ�����ڴ�ռ���ڴ��ͷ�
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
	int i, j, * p = NWL, n = 0;//ָ��pָ��NWL�׵�ַ
	printf("\t\t\t\tƽ����ܽṹ����\n");
	PHaaa();
	printf("\t�ڵ�����=%d\t�̶��ڵ�����=%d\t�ɶ��ڵ�����=%d\t�˼���=%d\n", TNN, NFIN, NFRN, NOR);
	PHaaa();
	printf("\t�ڵ��\t\tx����\ty����\n");
	for (i = 1; i <= TNN; i++)
		printf("\t%d\t\t%5.4f\t\t%5.4f\n", i, XCN[i - 1], YCN[i - 1]);
	PHaaa();
	printf("\t�˼���\t\tʼ�˽ڵ��\tĩ�˽ڵ��\t�����ն�EA\n");
	for (i = 1; i <= NOR; i++)
		printf("\t%d\t\t%d\t\t%d\t\t%5.4f\n", i, BNR[i - 1], ENR[i - 1], TSR[i - 1]);
	PHaaa();
	printf("\t���ؽڵ��\tx����\t\ty����\n");

	for (i = NFIN; i < TNN; i++)
	{
		printf("\t%d\t\t%5.4f\t\t%5.4f\n", i + 1, pp[2 * (i - NFIN)], pp[2 * (i - NFIN) + 1]);
	}
	PHaaa();
	printf("������\t���λ�����:\n");
	PHaaa();
	printf("\t�ڵ��\t\tλ��x����\tλ��y����\n");
	for (i = 0; i < NFRN; i++)
	{
		printf("\t%d\t\t%5.7f\t%5.7f\n", NFIN + i + 1, DON[2 * i], DON[2 * i + 1]);
	}
	PHaaa();
	printf("������\t�˼����������\n");
	PHaaa();
	printf("\t�˼���\t\t����ϵ��EA\t\t\t�˳�l\t\t�˼�����\t\t�˼�Ӧ��\n");
	for (i = 1; i <= NOR; i++)
	{
		printf("\t%d\t\t%5.4f\t\t\t%5.4f\t\t%5.7f\t\t%5.7f\t\n", i, TSR[i - 1], LCS[0][i - 1], IFR[i - 1], sigma[i - 1]);
	}
	PHaaa();

	printf("������\t�˼�ѹ���ȶ���Ӧ��У�˽��:\n");
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
		printf("\n\t\t���и˼�����ȫ\n");
	}
	else
		printf("\t\t����ΪΣ�ո˼��ţ�\n");

	PHaaa();

	printf("\n\t\t\t\t��л����ʹ��!\n\n\n\n\n\n");
}