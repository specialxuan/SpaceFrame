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
    fp = fopen("sf_test.csv", "w");
    
    int n = 0;
    scanf("%d", &n);

    fprintf(fp, "Stress Test\n");
    fprintf(fp, "TNN,%d,\n", 2 * n);
    fprintf(fp, "NFIN,%d,\n", n);
    fprintf(fp, "NOR,%d,\n", n);
    fprintf(fp, "NOL,%d,\n", n);
    fprintf(fp, "NOS,%d,\n", n);
    
    fprintf(fp, "XCN,");
    for (int i = 0; i < 2 * n; i++)
    {
        fprintf(fp, "0,");
    }
    fprintf(fp, "\n");

    fprintf(fp, "YCN,");
    for (int i = 0; i < n; i++)
    {
        fprintf(fp, "%d,", i);
    }
    for (int i = n; i < 2 * n; i++)
    {
        fprintf(fp, "%d,", i - n);
    }
    fprintf(fp, "\n");

    fprintf(fp, "ZCN,");
    for (int i = 0; i < n; i++)
    {
        fprintf(fp, "0,");
    }
    for (int i = n; i < 2 * n; i++)
    {
        fprintf(fp, "1,");
    }
    fprintf(fp, "\n");


    fprintf(fp, "BNR,");
    for (int i = 0; i < n; i++)
    {
        fprintf(fp, "%d,", i + 1);
    }
    fprintf(fp, "\n");

    fprintf(fp, "ENR,");
    for (int i = 0; i < n; i++)
    {
        fprintf(fp, "%d,", n + i + 1);
    }
    fprintf(fp, "\n");

    fprintf(fp, "ELASTIC,");
    for (int i = 0; i < n; i++)
    {
        fprintf(fp, "%d,", 210000000 + 100000 * (rand() % 1000));
    }
    fprintf(fp, "\n");

    fprintf(fp, "SHEAR,");
    for (int i = 0; i < n; i++)
    {
        fprintf(fp, "%d,", 80769000);
    }
    fprintf(fp, "\n");

    fprintf(fp, "AREA,");
    for (int i = 0; i < n; i++)
    {
        fprintf(fp, "%f,", 0.007854);
    }
    fprintf(fp, "\n");

    fprintf(fp, "IMY,");
    for (int i = 0; i < n; i++)
    {
        fprintf(fp, "%.10f,", 0.0000040001 + 0.0000000001 * (rand() % 10000));
    }
    fprintf(fp, "\n");

    fprintf(fp, "IMZ,");
    for (int i = 0; i < n; i++)
    {
        fprintf(fp, "%.10f,", 0.0000040001 + 0.0000000001 * (rand() % 10000));
    }
    fprintf(fp, "\n");

    fprintf(fp, "THETA,");
    for (int i = 0; i < n; i++)
    {
        fprintf(fp, "%d,", 0);
    }
    fprintf(fp, "\n");

    fprintf(fp, "NRL,");
    for (int i = 0; i < n; i++)
    {
        fprintf(fp, "%d,", i + 1);
    }
    fprintf(fp, "\n");

    fprintf(fp, "PLI,");
    for (int i = 0; i < n; i++)
    {
        fprintf(fp, "%d,", 0);
    }
    fprintf(fp, "\n");

    fprintf(fp, "KOL,");
    for (int i = 0; i < n; i++)
    {
        fprintf(fp, "%d,", 1);
    }
    fprintf(fp, "\n");

    fprintf(fp, "VOL,");
    for (int i = 0; i < n; i++)
    {
        fprintf(fp, "%d,", 1);
    }
    fprintf(fp, "\n");    

    fprintf(fp, "DLB,");
    for (int i = 0; i < n; i++)
    {
        fprintf(fp, "%d,", 1);
    }
    fprintf(fp, "\n");

    fprintf(fp, "NRS,");
    for (int i = 0; i < n; i++)
    {
        fprintf(fp, "%d,", i + 1);
    }
    fprintf(fp, "\n");

    fprintf(fp, "DSB,");
    for (int i = 0; i < n; i++)
    {
        fprintf(fp, "%d,", 0);
    }
    fprintf(fp, "\nEND,");

    return 0;
}