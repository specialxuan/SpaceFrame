#include "SpaceFrame.h"

int main()
{
    clock_t start1 = 0, end1 = 0;
    DWORD start, end;
    start1 = clock();
    start = GetTickCount();

    SpaceFrame Frame;
    Frame.sfCircularStructure(5, 5, 5);
    Frame.sfInput();
    // Frame.sfCalculate();
    // Frame.sfOutput();
    // Frame.~SpaceFrame();
    // Frame.sfInput();
    Frame.sfCalculate(true, true, 1e-15);
    Frame.sfOutput();
    // SpaceFrame Frame2(Frame);
    // Frame2.sfInput();
    // Frame2.sfCalculate();
    // Frame2.sfOutput();
    end1 = clock();
    printf("time = %f\n", (double)(end1 - start1) / CLOCKS_PER_SEC);
    end = GetTickCount();
    printf("realtime=%f\n", (double)(end - start) / 1000);

    system("pause");

    return 0;
}
