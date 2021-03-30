#include "SpaceFrame.h"
#include "SpaceFrame v4.4.1.h"

int main()
{
    clock_t start1 = 0, end1 = 0;
    DWORD start, end;
    start1 = clock();
    start = GetTickCount();

    SpaceFrame Frame;
    // Frame.sfCircularStructure(2, 2, 1);
    Frame.sfInput();
    Frame.sfCalculate(true, true, 1e-15);
    Frame.sfOutput(true);

    SpaceFrame Frame2(Frame);
    // Frame2.sfInput();
    // Frame2.sfCalculate(true, true, 1e-15);
    // Frame2.sfOutput(true);
    end1 = clock();
    cout << (double)(end1 - start1) / CLOCKS_PER_SEC << endl;
    end = GetTickCount();
    cout << (double)(end - start) / 1000 << endl;

    // system("pause");

    return 0;
}
