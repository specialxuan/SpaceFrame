#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <mpi.h>
#include <omp.h>
#include <time.h>
void add(int*, char***);
int main1(int argc, char* argv[])
{
	time_t start, end;
	int* p_argc = &argc;
	char*** p_argv = &argv;
	start = clock();
	add(p_argc, p_argv);
	end = clock();
	printf("time=%f\n", ((double)end - start) / CLK_TCK);
	return 0;
}

void add(int *p1, char ***p2)
{
	double global_solutions=0;
	int id, p;
	double A = 0;
	MPI_Init(p1, p2);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
    //#pragma omp parallel for
    for (int i = id; i < 100000000; i+=p)
	{
		for (int j = 0; j < 20; j++)
		{
		A ++;
		}
	}
	MPI_Reduce(&A, &global_solutions, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	printf("Process %d is done\n", id);
	fflush(stdout);
	MPI_Finalize();
	if (id==0) printf("%f", global_solutions);
}