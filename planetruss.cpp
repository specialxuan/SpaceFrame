#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>
#include <Windows.h>

int TNN;  // total number of nodes
int NFIN; // number of fixed nodes
int NFRN; // number of free nodes
int NOR;  // total number of rods

double* XCN; // X coordinates
double* YCN; // Y coordinates
int* BNR;	 // begining numbers of rods
int* ENR;	 // ending numbers of rods
double* TSR; // tensile stiffness of rods

int* NWL;	 // numbers of nods with load
double* XCL; // X axis component load
double* YCL; // Y axis component load

double* LCS; // length cosine sine
double* DON;  // displacement of nodes
double* IFR;  // internal forces of rods

double* ROU;	 // density of rods
double* AREA;	 // sectional area
double g = 9.81; // acceleration of gravity

// read data from .txt
bool PTRead();
// build total stiff matrix
double* PTBuildTotalStiff();
// calculate the length, cosine and sine of rods
bool PTLCosSin(int);
// build unit stiffness matrix
bool PTBuildUnitStiff(int, double*);
// match the number of displacement
int* PTI0J0(int);
// calculate self weight
bool PTselfweight();
// build the vector of load
double* PTBuildLoadVector();
// solve matrix equation
bool PTCholesky(double*, double*, double*, int);
// calculate the force inside rods
bool PTRodForce();
// print "-------------------"
bool PTPrintLine();
// print "*******************"
bool PTPrintLine2();
// print result
bool PTPrint();
// print error
bool PTPrintError(int);

bool solve_conjugate_gradient(double* A, double* b, double* x, int N);

int main()
{
	double* kk = 0, * pp = 0;
	clock_t start1 = 0, end1 = 0;

    int coreNum = omp_get_num_procs();//获得处理器个数
    printf("%d\n", coreNum);

#pragma omp parallel
    {
        printf("Hello OpenMP!ThreadID=%d\n",omp_get_thread_num());
    }

    printf("Welcome to use the calculator of plane truss!\nPress any key to start");
	char value = getchar();
    start1 = clock();
    DWORD start,end;  
    start= GetTickCount();  

	PTPrintLine(); //"---------------------------------"
	if (PTRead())
	{
		PTPrintError(1);
		printf("\nPress any key to exit\n");
		value = getchar();

		return 1;
	}
	else
		printf("Data input succeeded!\n");
	if ((kk = PTBuildTotalStiff()) == 0)
	{
		PTPrintError(2);
		printf("\nPress any key to exit\n");
		value = getchar();

		return 1;
	}
	else
		printf("Building total stiffness matrix succeeded!\n");
	if ((pp = PTBuildLoadVector()) == 0)
	{
		PTPrintError(3);
		printf("\nPress any key to exit\n");
		value = getchar();

		return 1;
	}
	else
		printf("Building load vector succeeded!\n");
	if (solve_conjugate_gradient(kk, pp, DON, 2 * NFRN))
	{
		PTPrintError(4);
		printf("\nPress any key to exit\n");
		value = getchar();

		return 1;
	}
	else
		printf("Solving equation succeeded!\n");

	free(kk);
	free(pp);

	if (PTRodForce())
	{
		PTPrintError(5);
		printf("\nPress any key to exit\n");
		value = getchar();

		return 1;
	}
	else
		printf("Calculating internal force succeeded!\n");
	//PTPrint();

	end1 = clock();
	printf("time = %f\n", (double)(end1 - start1) / CLOCKS_PER_SEC);
    end= GetTickCount();  
    printf("realtime=%f\n",(double)(end-start)/1000); 

	printf("Press any key to exit\n");
	value = getchar();




	return 0;
}

bool PTRead()
{
	FILE* fp;
	fp = fopen("structure_data_1.txt", "r");

	fseek(fp, 22L, 0);
	if (fscanf(fp, "%5d", &TNN) == EOF)
	{
		PTPrintError(9);
		return 1;
	}
	fseek(fp, 24L, 1);
	if (fscanf(fp, "%5d", &NFIN) == EOF)
	{
		PTPrintError(9);
		return 1;
	}
	fseek(fp, 17L, 1);
	if (fscanf(fp, "%5d", &NOR) == EOF)
	{
		PTPrintError(9);
		return 1;
	}

	NFRN = TNN - NFIN;

	XCN = (double*)malloc(TNN * sizeof(double));
	if (XCN != NULL)
		memset(XCN, 0, TNN * sizeof(double));
	YCN = (double*)malloc(TNN * sizeof(double));
	if (YCN != NULL)
		memset(YCN, 0, TNN * sizeof(double));
	BNR = (int*)malloc(NOR * sizeof(int));
	if (BNR != NULL)
		memset(BNR, 0, NOR * sizeof(int));
	ENR = (int*)malloc(NOR * sizeof(int));
	if (ENR != NULL)
		memset(ENR, 0, NOR * sizeof(int));
	TSR = (double*)malloc(NOR * sizeof(double));
	if (TSR != NULL)
		memset(TSR, 0, NOR * sizeof(double));
	NWL = (int*)malloc(TNN * sizeof(int));
	if (NWL != NULL)
		memset(NWL, 0, TNN * sizeof(int));
	XCL = (double*)malloc(TNN * sizeof(double));
	if (XCL != NULL)
		memset(XCL, 0, TNN * sizeof(double));
	YCL = (double*)malloc(TNN * sizeof(double));
	if (YCL != NULL)
		memset(YCL, 0, TNN * sizeof(double));
	DON = (double*)malloc(2 * NFRN * sizeof(double));
	if (DON != NULL)
		memset(DON, 0, 2 * NFRN * sizeof(double));
	IFR = (double*)malloc(NOR * sizeof(double));
	if (IFR != NULL)
		memset(IFR, 0, NOR * sizeof(double));
	ROU = (double*)malloc(NOR * sizeof(double));
	if (ROU != NULL)
		memset(ROU, 0, NOR * sizeof(double));
	AREA = (double*)malloc(NOR * sizeof(double));
	if (AREA != NULL)
		memset(AREA, 0, NOR * sizeof(double));

	fseek(fp, 2L, 1);
	for (int i = 0; i < 10; i++)
	{
		fseek(fp, 4L, 1);
		char c = 0;
		int j = 0;

		do
		{
			switch (i)
			{
			case 0:
				if (fscanf(fp, "%lf", XCN + j) == EOF)
				{
					PTPrintError(9);
					return 1;
				}
				break;
			case 1:
				if (fscanf(fp, "%lf", YCN + j) == EOF)
				{
					PTPrintError(9);
					return 1;
				}
				break;
			case 2:
				if (fscanf(fp, "%5d", BNR + j) == EOF)
				{
					PTPrintError(9);
					return 1;
				}
				break;
			case 3:
				if (fscanf(fp, "%5d", ENR + j) == EOF)
				{
					PTPrintError(9);
					return 1;
				}
				break;
			case 4:
				if (fscanf(fp, "%lf", TSR + j) == EOF)
				{
					PTPrintError(9);
					return 1;
				}
				break;
			case 5:
				if (fscanf(fp, "%5d", NWL + j) == EOF)
				{
					PTPrintError(9);
					return 1;
				}
				break;
			case 6:
				if (fscanf(fp, "%lf", XCL + j) == EOF)
				{
					PTPrintError(9);
					return 1;
				}
				break;
			case 7:
				if (fscanf(fp, "%lf", YCL + j) == EOF)
				{
					PTPrintError(9);
					return 1;
				}
				break;
			case 8:
				if (fscanf(fp, "%lf", ROU + j) == EOF)
				{
					PTPrintError(9);
					return 1;
				}
				break;
			case 9:
				if (fscanf(fp, "%lf", AREA + j) == EOF)
				{
					PTPrintError(9);
					return 1;
				}
				break;
			default:
				PTPrintError(9);
				return 1;
				break;
			}
			fscanf(fp, "%c", &c);
			j++;
		} while (c != '\n');
	}

	return 0;
}

double* PTBuildTotalStiff()
{
	double* kk = 0, * us = 0;
	const int nods = 2 * NFRN;
	int* p = 0;

	kk = (double*)malloc(nods * nods * sizeof(double));
	if (kk != NULL)
		memset(kk, 0, nods * nods * sizeof(double));

	us = (double*)malloc(2 * 2 * sizeof(double));
	if (us != NULL)
		memset(us, 0, 2 * 2 * sizeof(double));

	LCS = (double*)malloc(3 * NOR * sizeof(double));
	if (LCS != NULL)
		memset(LCS, 0, 3 * NOR * sizeof(double));

	for (int num = 0; num < NOR; num++)
	{
		if (PTLCosSin(num))
		{
			PTPrintError(6);
			return 0;
		}
		if (PTBuildUnitStiff(num, us))
		{
			PTPrintError(7);
			return 0;
		}
		if ((p = PTI0J0(num)) == 0)
		{
			PTPrintError(8);
			return 0;
		}
		for (int i = 0; i < 2; i++)
		{
			if (p[i] >= 0)
			{
                
				for (int m = 0; m < 2; m++)
					for (int n = 0; n < 2; n++)
					{
						kk[(p[i] + m) * nods + (p[i] + n)] += us[m * 2 + n];
					}
			}
		}
		if (p[0] >= 0 && p[1] >= 0)
		{
        
			for (int i = 0; i < 2; i++)
			{
				for (int m = 0; m < 2; m++)
					for (int n = 0; n < 2; n++)
					{
						kk[(p[i] + m) * nods + (p[1 - i] + n)] -= us[m * 2 + n];
					}
			}
		}
		free(p);
	}
	free(us);

	return kk;
}

bool PTLCosSin(int num)
{
	int i = BNR[num] - 1, j = ENR[num] - 1;
	LCS[1 * NOR + num] = XCN[j] - XCN[i];
	LCS[2 * NOR + num] = YCN[j] - YCN[i];
	LCS[0 * NOR + num] = sqrt(LCS[1 * NOR + num] * LCS[1 * NOR + num] + LCS[2 * NOR + num] * LCS[2 * NOR + num]);
	LCS[1 * NOR + num] /= LCS[0 * NOR + num];
	LCS[2 * NOR + num] /= LCS[0 * NOR + num];

	return 0;
}

bool PTBuildUnitStiff(int num, double* us)
{
	double rd = TSR[num] / LCS[0 * NOR + num];
   
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			us[i * 2 + j] = LCS[(i + 1) * NOR + num] * LCS[(j + 1) * NOR + num] * rd;
		}
	}

	return 0;
}

int* PTI0J0(int num)
{
	int* ij = 0;
	ij = (int*)malloc(2 * sizeof(int));

	ij[0] = 2 * (BNR[num] - NFIN - 1);
	ij[1] = 2 * (ENR[num] - NFIN - 1);

	return ij;
}

bool PTselfweight()
{
	double w = 0;
	int bl = 0, br = 0;

    #pragma omp parallel for
	for (int k = 0; k < NOR; k++)
	{
		bl = BNR[k];
		br = ENR[k];
		w = ROU[k] * AREA[k] * LCS[0 * NOR + k] * g;
		YCL[bl - 1] -= 0.5 * w;
		YCL[br - 1] -= 0.5 * w;
	}

	return 0;
}

double* PTBuildLoadVector()
{
	PTselfweight();

	double* pp = 0;
	pp = (double*)malloc(2 * NFRN * sizeof(double));
	if (pp != NULL)
		memset(pp, 0, 2 * NFRN * sizeof(double));
    #pragma omp parallel for
	for (int i = NFIN; i <= TNN; i++)
	{
		pp[2 * (i - NFIN)] = XCL[i];
		pp[2 * (i - NFIN) + 1] = YCL[i];
	}

	return pp;
}

bool PTCholesky(double* A, double* b, double* x, int n)
{
	if (A == NULL)
	{
		PTPrintError(4);
		return 1;
	}
	else if (b == NULL)
	{
		PTPrintError(4);
		return 1;
	}
	else if (x == NULL)
	{
		PTPrintError(4);
		return 1;
	}
	else if (n == 0)
	{
		PTPrintError(4);
		return 1;
	}

	double sum = 0, * L = 0, * D = 0, * y = 0;

	L = (double*)malloc(n * n * sizeof(double));
	if (L != NULL)
		memset(L, 0, n * n * sizeof(double));
	for (int i = 0; i < n; i++)
	{
		L[i * n + i] = 1;
	}
	D = (double*)malloc(n * sizeof(double));
	memset(D, 0, n * sizeof(double));
	y = (double*)malloc(n * sizeof(double));
	if (y != NULL)
		memset(y, 0, n * sizeof(double));

	// factorize matrix A
	D[0] = A[0 * n + 0];
	for (int i = 1; i < n; i++)
	{
		for (int j = 0; j < i; j++)
		{
			sum = 0;
			for (int k = 0; k < j; k++)
			{
				sum = sum + L[i * n + k] * D[k] * L[j * n + k];
			}
			L[i * n + j] = (A[i * n + j] - sum) / D[j];
		}
		sum = 0;
		for (int k = 0; k < i; k++)
		{
			sum = sum + L[i * n + k] * L[i * n + k] * D[k];
		}
		D[i] = A[i * n + i] - sum;
	}

	// solve y by Ly=b
	y[0] = b[0];
	for (int i = 1; i < n; i++)
	{
		sum = 0;
		for (int j = 0; j < i; j++)
		{
			sum = sum + L[i * n + j] * y[j];
		}
		y[i] = b[i] - sum;
	}

	// solve x by DL(t)=y
	x[n - 1] = y[n - 1] / D[n - 1];
	for (int i = n - 2; i >= 0; i--)
	{
		sum = 0;
		for (int j = i + 1; j < n; j++)
		{
			sum = sum + L[j * n + i] * x[j];
		}
		x[i] = y[i] / D[i] - sum;
	}

	free(L);
	free(D);
	free(y);

	return 0;
}

bool PTRodForce()
{
	int i = 0, j = 0, * p = 0;
	double d1[2] = { 0 }, d2[2] = { 0 }, rd = 0;

	for (int k = 0; k < NOR; k++)
	{
		p = PTI0J0(k);
		i = p[0];
		j = p[1];
		if (i < 0)
		{
			d1[0] = d1[1] = 0;
		}
		else
		{
			d1[0] = DON[i];
			d1[1] = DON[i + 1];
		}
		if (j < 0)
		{
			d2[0] = d2[1] = 0;
		}
		else
		{
			d2[0] = DON[j];
			d2[1] = DON[j + 1];
		}
		rd = TSR[k] / LCS[0 * NOR + k];
		if (IFR + k != NULL)
			*(IFR + k) = rd * (LCS[1 * NOR + k] * (d2[0] - d1[0]) + LCS[2 * NOR + k] * (d2[1] - d1[1]));
	}

	return 0;
}

bool PTPrintLine()
{
	printf("-------------------------------------------------------------------------\n");
	return 0;
}

bool PTPrintLine2()
{
	printf("*************************************************************************\n");
	return 0;
}

bool PTPrint()
{
	PTPrintLine();
	printf("| \t\t\tCalculation of Plane Truss\t\t\t|\n");
	PTPrintLine();
	printf("| TNN = %9d | NFIN = %8d | NFRN = %8d | NOR = %9d |\n", TNN, NFIN, NFRN, NOR);
	PTPrintLine();
	printf("|          Number |        X        |        Y        |                 |\n");
	for (int i = 0; i < TNN; i++)
	{
		printf("| %15d | %15.4f | %15.4f |                 |\n", i, XCN[i], YCN[i]);
	}
	PTPrintLine();
	printf("|          Number |       Beginning |          Ending |              EA |\n");
	for (int i = 0; i < NOR; i++)
	{
		printf("| %15d | %15d | %15d | %15.4f |\n", i, BNR[i], ENR[i], TSR[i]);
	}
	PTPrintLine();
	printf("| Nodes with load |   X component   |   Y component   |                 |\n");
	for (int i = 0; i < TNN; i++)
	{
		printf("| %15d | %15.4f | %15.4f |                 |\n", i + 1, XCL[i], YCL[i]);
	}
	PTPrintLine();
	printf("| \t\t\tResult\tDisplacement of nodes\t\t\t|\n");
	PTPrintLine();
	printf("|          Number | X displacement  | Y displacement  |                 |\n");
	for (int i = 0; i < NFRN; i++)
	{
		printf("| %15d | %15.7f | %15.7f |                 |\n", NFIN + i, DON[2 * i], DON[2 * i + 1]);
	}
	PTPrintLine();
	printf("| \t\t\tResult\tInternal force of rods\t\t\t|\n");
	PTPrintLine();
	printf("|          Number | Section factor  | Length of rods  | Internal force  |\n");
	for (int i = 0; i < NOR; i++)
	{
		printf("| %15d | %15.4f | %15.4f | %15.7f |\n", i + 1, TSR[i], LCS[0 * NOR + i], IFR[i]);
	}
	PTPrintLine();
	printf("Thank you for your usage!\n");

	free(XCN);
	free(YCN);
	free(BNR);
	free(ENR);
	free(TSR);
	free(NWL);
	free(XCL);
	free(YCL);
	free(DON);
	free(IFR);
	free(ROU);
	free(AREA);

	return 0;
}

bool PTPrintError(int error)
{
	printf("ERROR:\t");
	switch (error)
	{
	case 1:
		printf("Data input failed!\n");
		break;
	case 2:
		printf("Building total stiffness matrix failed!\n");
		break;
	case 3:
		printf("Building load vector failed!\n");
		break;
	case 4:
		printf("Solving equation failed!\n");
		break;
	case 5:
		printf("Calculating internal force failed!\n");
		break;
	case 6:
		printf("Calculating length, cosine and sine failed!\n");
		break;
	case 7:
		printf("Building unit stiffness matrix failed!\n");
		break;
	case 8:
		printf("Matching displacement failed!\n");
		break;
	case 9:
		printf("Data format is wrong!\n");
	default:
		break;
	}
	printf("There is at least one error in your file, please check it and try it one more time.\n");

	return 0;
}

bool solve_conjugate_gradient(double* A, double* b, double* x, int N)
{
	double* r, * p , * z, tol = 1e-16;
	double gamma, gamma_new, alpha, beta;

	r = (double*)malloc(N * sizeof(double));
	p = (double*)malloc(N * sizeof(double));
	z = (double*)malloc(N * sizeof(double));

	// x = [0 ... 0]
	// r = b - A * x
	// p = r
	// gamma = r' * r
	gamma = 0.0;
    #pragma omp parallel for reduction(+:gamma)
	for (int i = 0; i < N; ++i) {
		x[i] = 0.0;
		r[i] = b[i];
		p[i] = r[i];
		gamma += r[i] * r[i];
	}
    
	for (int n = 0; n < N; ++n) 
    {
		// z = A * p
        #pragma omp parallel for
		for (int i = 0; i < N; ++i) 
        {
			z[i] = 0.0;
			for (int j = 0; j < N; ++j)
				z[i] += A[i * N + j] * p[j];
        }

		// alpha = gamma / (p' * z)
		alpha = 0.0;
        #pragma omp parallel for reduction(+:alpha)
		for (int i = 0; i < N; ++i)
        {
            // printf("%d\n", i);
            alpha += p[i] * z[i];
        }
		alpha = gamma / alpha;

		// x = x + alpha * p
		// r = r - alpha * z
		// gamma_new = r' * r
		gamma_new = 0.0;
        #pragma omp parallel for reduction(+:gamma_new)
		for (int i = 0; i < N; ++i) 
        {
			x[i] += alpha * p[i];
			r[i] -= alpha * z[i];
			gamma_new += r[i] * r[i];
		}

		if (sqrt(gamma_new) < tol)
			break;

		beta = gamma_new / gamma;

		// p = r + (gamma_new / gamma) * p;
        #pragma omp parallel for
		for (int i = 0; i < N; ++i)
        {
			p[i] = r[i] + beta * p[i];
        }

		// gamma = gamma_new
		gamma = gamma_new;
	}

	free(r);
	free(p);
	free(z);

	return 0;
}