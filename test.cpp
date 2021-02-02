#include "SpaceFrame.h"

int main()
{
    FILE *fp = fopen("source&result/testresult.csv", "w");
    SpaceFrame Frame;
    DWORD start, end;
    fprintf(fp, " DOF      , PB off   , PB on    ,\n");

    for (int i = 3; i < 16; i++)
    {
        Frame.sfCircularStructure(i, i, i);
        fprintf(fp, "%9d ,", ((i + 1) * (i + 1) * (i + 1) - (i + 1) * (i + 1)) * 6);

        start = GetTickCount();
        Frame.sfInput();
        Frame.sfCalculate(true, false);
        Frame.sfOutput();
        end = GetTickCount();
        fprintf(fp, "%9.2f ,", (double)(end - start) / 1000);

        start = GetTickCount();
        Frame.sfInput();
        Frame.sfCalculate();
        Frame.sfOutput();
        end = GetTickCount();
        fprintf(fp, "%9.2f ,\n", (double)(end - start) / 1000);
    }

    fclose(fp);
    fp = NULL;
    return 0;
}