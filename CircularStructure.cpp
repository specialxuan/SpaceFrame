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

    int n = 0, m = 0, l = 0;
    scanf("%d %d %d", &m, &n, &l);

    fprintf(fp, "Stress Test, degree of freedom is %d,\n", ((m + 1) * (n + 1) * (l + 1) - (m + 1)* (n + 1)) * 6);
    fprintf(fp, "TNN,%d,\n", (m + 1) * (n + 1) * (l + 1));
    fprintf(fp, "NFIN,%d,\n", (m + 1)* (n + 1));
    int nor = ((2 * m + 1) * (2 * n + 1) - m * n) * l;
    fprintf(fp, "NOR,%d,\n", nor);
    fprintf(fp, "NOL,%d,\n", (m + 1) * (n + 1));
    fprintf(fp, "NOS,%d,\n", (m + 1) * (n + 1));
    
    fprintf(fp, "XCN,");
    for (int i = 0; i < l + 1; i++)
        for (int j = 0; j < n + 1; j++)
            for (int k = 0; k < m + 1; k++)
                fprintf(fp, "%d,", k);
    fprintf(fp, "\n");    

    fprintf(fp, "YCN,");
    for (int i = 0; i < l + 1; i++)
        for (int j = 0; j < n + 1; j++)
            for (int k = 0; k < m + 1; k++)
                fprintf(fp, "%d,", j);
    fprintf(fp, "\n");

    fprintf(fp, "ZCN,");
    for (int i = 0; i < l + 1; i++)
        for (int j = 0; j < n + 1; j++)
            for (int k = 0; k < m + 1; k++)
                fprintf(fp, "%d,", i);
    fprintf(fp, "\n");

    fprintf(fp, "BNR,");
    for (int i = 0; i < l; i++)
    {
        for (int j = 0; j < (m + 1) * (n + 1); j++)
            fprintf(fp, "%d,", j + 1 + i * (m + 1) * (n + 1));
        for (int j = 0; j < (m + 1) * n; j++)
            fprintf(fp, "%d,", j + 1 + (i + 1) * (m + 1) * (n + 1));
        for (int j = 0; j < n + 1; j++)
            for (int k = 0; k < m; k++)
                fprintf(fp, "%d,", k + 1 + j * (m + 1) + (i + 1) * (m + 1) * (n + 1));
    }
    fprintf(fp, "\n");

    fprintf(fp, "ENR,");
    for (int i = 0; i < l; i++)
    {
        for (int j = 0; j < (m + 1) * (n + 1); j++)
            fprintf(fp, "%d,", j + 1 + (i + 1) * (m + 1) * (n + 1));
        for (int j = 0; j < (m + 1) * n; j++)
            fprintf(fp, "%d,", j + 1 + m + 1 + (i + 1) * (m + 1) * (n + 1));
        for (int j = 0; j < n + 1; j++)
            for (int k = 0; k < m; k++)
                fprintf(fp, "%d,", k + 2 + j * (m + 1) + (i + 1) * (m + 1) * (n + 1));
    }
    fprintf(fp, "\n");

    fprintf(fp, "ELASTIC,");
    for (int i = 0; i < nor; i++)
    {
        fprintf(fp, "%d,", 210000000 + 100000 * (rand() % 1000));
    }
    fprintf(fp, "\n");

    fprintf(fp, "SHEAR,");
    for (int i = 0; i < nor; i++)
    {
        fprintf(fp, "%d,", 80769000);
    }
    fprintf(fp, "\n");

    fprintf(fp, "AREA,");
    for (int i = 0; i < nor; i++)
    {
        fprintf(fp, "%f,", 0.007854);
    }
    fprintf(fp, "\n");

    fprintf(fp, "IMY,");
    for (int i = 0; i < nor; i++)
    {
        fprintf(fp, "%.10f,", 0.0000040001 + 0.0000000001 * (rand() % 10000));
    }
    fprintf(fp, "\n");

    fprintf(fp, "IMZ,");
    for (int i = 0; i < nor; i++)
    {
        fprintf(fp, "%.10f,", 0.0000040001 + 0.0000000001 * (rand() % 10000));
    }
    fprintf(fp, "\n");

    fprintf(fp, "THETA,");
    for (int i = 0; i < nor; i++)
    {
        fprintf(fp, "%d,", 0);
    }
    fprintf(fp, "\n");

    fprintf(fp, "NRL,");
    for (int i = 0; i < (m + 1) * (n + 1); i++)
    {
        fprintf(fp, "%d,", i + 1 + ((2 * m + 1) * (2 * n + 1) - m * n) * (l - 1));
    }
    fprintf(fp, "\n");

    fprintf(fp, "PLI,");
    for (int i = 0; i < (m + 1) * (n + 1); i++)
    {
        fprintf(fp, "%d,", 0);
    }
    fprintf(fp, "\n");

    fprintf(fp, "KOL,");
    for (int i = 0; i < (m + 1) * (n + 1); i++)
    {
        fprintf(fp, "%d,", 3);
    }
    fprintf(fp, "\n");

    fprintf(fp, "VOL,");
    for (int i = 0; i < (m + 1) * (n + 1); i++)
    {
        fprintf(fp, "%d,", 1000 + rand() % 1000);
    }
    fprintf(fp, "\n");    

    fprintf(fp, "DLB,");
    for (int i = 0; i < (m + 1) * (n + 1); i++)
    {
        fprintf(fp, "%d,", 1);
    }
    fprintf(fp, "\n");

    fprintf(fp, "NRS,");
    for (int i = 0; i < (m + 1) * (n + 1); i++)
    {
        fprintf(fp, "%d,", i + 1);
    }
    fprintf(fp, "\n");

    fprintf(fp, "DSB,");
    for (int i = 0; i < (m + 1) * (n + 1); i++)
    {
        fprintf(fp, "%f,", 0.5);
    }
    fprintf(fp, "\nEND,");

    return 0;
}