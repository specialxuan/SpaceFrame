#include "SpaceFrame.h"

int main()
{
    FILE *fp = fopen("source&result/testresult.csv", "w");
    SpaceFrame Frame;
    DWORD start, end;
    fprintf(fp, " DOF      , 1e-15    , 1e-10    , 1e-5     ,\n");

    for (int i = 3; i <= 15; i++)
    {
        // for (int j = 2; j < 5; j++)
        // {
        //     for (int k = 1; k < 5; k++)
        //     {
        //         Frame.sfCircularStructure(i, j, k);
        //         fprintf(fp, "%9d ,", ((i + 1) * (j + 1) * (k + 1) - (i + 1) * (j + 1)) * 6);

        //         start = GetTickCount();
        //         Frame.sfInput();
        //         Frame.sfCalculate(true, false);
        //         Frame.sfOutput();
        //         end = GetTickCount();
        //         fprintf(fp, "%9.2f ,", (double)(end - start) / 1000);

        //         start = GetTickCount();
        //         Frame.sfInput();
        //         Frame.sfCalculate(true, false);
        //         Frame.sfOutput();
        //         end = GetTickCount();
        //         fprintf(fp, "%9.2f ,\n", (double)(end - start) / 1000);
        //     }
        // }
        Frame.sfCircularStructure(i, i, i);
        fprintf(fp, "%9d ,", ((i + 1) * (i + 1) * (i + 1) - (i + 1) * (i + 1)) * 6);

        start = GetTickCount();
        Frame.sfInput();
        Frame.sfCalculate(true, false, 1e-15);
        Frame.sfOutput();
        end = GetTickCount();
        fprintf(fp, "%9.2f ,", (double)(end - start) / 1000);

        start = GetTickCount();
        Frame.sfInput();
        Frame.sfCalculate(true, false, 1e-10);
        Frame.sfOutput();
        end = GetTickCount();
        fprintf(fp, "%9.2f ,", (double)(end - start) / 1000);

        start = GetTickCount();
        Frame.sfInput();
        Frame.sfCalculate(true, false, 1e-5);
        Frame.sfOutput();
        end = GetTickCount();
        fprintf(fp, "%9.2f ,\n", (double)(end - start) / 1000);
    }

    fclose(fp);
    fp = NULL;
    return 0;
}