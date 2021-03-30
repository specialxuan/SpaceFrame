#include "SpaceFrame v4.4.1.h"
#include "SpaceFrame.h"

int main()
{
    ofstream fout("source&result/testresult.csv", ios::out);
    SpaceFrame Frame;
    SpaceFrame_v4 Frame_v4;
    DWORD start, end, start1, end1, start2, end2;
    fout << " DOF      , 1e-15    , 1e-10    , 1e-5     ,\n";

    for (int i = 3; i <= 6; i++)
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
        fout << setw(10) << ((i + 1) * (i + 1) * (i + 1) - (i + 1) * (i + 1)) * 6 << ",";

        start1 = GetTickCount();
        Frame.sfInput();
        Frame.sfCalculate(true, false, 1e-15);
        Frame.sfOutput();
        end1 = GetTickCount();
        fout << setw(10) << (double)(end1 - start1) / 1000 << ",";
        cout << setw(10) << (double)(end1 - start1) / 1000 << ",";

        start2 = GetTickCount();
        Frame_v4.sfInput();
        Frame_v4.sfCalculate(true, false, 1e-15);
        Frame_v4.sfOutput();
        end2 = GetTickCount();
        fout << setw(10) << (double)(end2 - start2) / 1000 << ",";
        cout << setw(10) << (double)(end2 - start2) / 1000 << ",";

        // start = GetTickCount();
        // Frame.sfInput();
        // Frame.sfCalculate(true, false, 1e-15);
        // Frame.sfOutput();
        // end = GetTickCount();
        fout << setw(10) << (double)((end1 - start1) - (end2 - start2)) / 1000 << ",\n";
    }

    fout.close();

    return 0;
}