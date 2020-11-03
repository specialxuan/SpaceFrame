#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

int main()
{
    FILE *fp = 0;
    fp = fopen("structure_data_1.txt", "w");
    
    int n = 0;
    scanf("%d", &n);

    fprintf(fp, "Total number of nodes:%d\n", 2 * n + 1);
    fprintf(fp, "Number of fixed nodes:%d\n", n + 1);
    fprintf(fp, "Number of rods:%d\n", 2 * n);
    
    fprintf(fp, "XCN");
    for (int i = 0; i < n + 1; i++)
    {
        fprintf(fp, ",%d", 2 * i);
    }
    for (int i = n + 1; i < 2 * n + 1; i++)
    {
        fprintf(fp, ",%d", 2 * (i - n - 1) + 1);
    }
    fprintf(fp, "\n");

    fprintf(fp, "YCN");
    for (int i = 0; i < n + 1; i++)
    {
        fprintf(fp, ",0");
    }
    for (int i = n + 1; i < 2 * n + 1; i++)
    {
        fprintf(fp, ",1");
    }
    fprintf(fp, "\n");

    fprintf(fp, "BNR");
    for (int i = 0; i < n; i++)
    {
        fprintf(fp, ",%d,%d", i + 1, i + 2);
    }
    fprintf(fp, "\n");

    fprintf(fp, "ENR");
    for (int i = 0; i < n; i++)
    {
        fprintf(fp, ",%d,%d", n + i + 2, n + i + 2);
    }
    fprintf(fp, "\n");

    fprintf(fp, "TSR");
    for (int i = 0; i < 2 * n; i += 2)
    {
        int tsr = 100 * rand() % (3000000) + 2000000;
        fprintf(fp, ",%d,%d", tsr, tsr);
    }
    fprintf(fp, "\n");

    fprintf(fp, "NWL");
    for (int i = 0; i < n + 1; i++)
    {
        fprintf(fp, ",%d", 0);
    }
    for (int i = n + 1; i < 2 * n + 1; i++)
    {
        fprintf(fp, ",%d", 1);
    }
    fprintf(fp, "\n");

    fprintf(fp, "XCL");
    for (int i = 0; i < 2 * n + 1; i++)
    {
        fprintf(fp, ",%d", 0);
    }
    fprintf(fp, "\n");

    fprintf(fp, "YCL");
    for (int i = 0; i < n + 1; i++)
    {
        fprintf(fp, ",%d", 0);
    }
    for (int i = n + 1; i < 2 * n + 1; i++)
    {
        fprintf(fp, ",%d", -10000);
    }
    fprintf(fp, "\n");

    fprintf(fp, "ROU");
    for (int i = 0; i < 2 * n; i++)
    {
        fprintf(fp, ",%d", 0);
    }
    fprintf(fp, "\n");    

    fprintf(fp, "ARA");
    for (int i = 0; i < 2 * n; i++)
    {
        fprintf(fp, ",%d", 1);
    }
    fprintf(fp, "\n");

    return 0;
}