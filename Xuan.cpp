#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define EPS 1e-15

int TNN;  //total number of nodes
int NFIN; //number of fixed nodes
int NFRN; //number of free nodes
int NOR;  //number of rods
int NOL;  //number of loads
int NOS;  //number of sections

double *XCN; //X coordinate of nodes
double *YCN; //Y coordinate of nodes
double *ZCN; //Z coordinate of nodes

int *BNR;        //the beginning node number of rods
int *ENR;        //the end node number of rods
double *ELASTIC; //elastic modulus
double *SHEAR;   //shear modulus
double *AREA;    //area
double *IMY;     //inertia moment of Y axis
double *IMZ;     //inertia moment of Z axis
double *THETA;   //theta the deflection angle of main inertia axis

int *NRL;    //the number of rods with load
int *PLI;    //the plane of the load's in
int *KOL;    //the kind of load
double *VOL; //the value of load
double *DLB; //the distance between load and the beginning node

int *NRS;    //the number of rod with section
double *DSB; //the distance between section and the beginning node

double *LCS; //the length, sine and cosine of rods
double *DON; //the displacement of nodes
double *IFS; //the internal force in the section
double *RFE; //the reaction force of the end node

int *iv;      //the location of diagonal element
int nsi;      //upper limit
int maxibdw;  //half bandwidth
int *peribdw; //bandwidth per line

//read data from .csv
bool sfInput();
//build total stiffness matrix
bool sfBuildTotalStiff(double *);
//calculate the length sine and cosine of rods
bool sfLCosSin();
//build unit stiffness matrix
bool sfBuildUnitStiff(int, int, double *);
//build local stiffness matrix
bool sfBuildLocalStiff(int, int, double *);
//build transpose matrix
bool sfBuildTrans(int, double *);
//build load vector
bool sfBuildLoadVector(double *);
//calculate reaction force
bool sfReactionForce(int, double *, double *);
//solve equation of matrix by cholesky
bool sfCholesky(double *, double *, double *, int);
//solve equation of matrix by conjugate gradient
bool solve_conjugate_gradient(double *A, double *b, double *x, int N);
//calculate internal force of rods
bool sfInternalForce(int, int, double);
//calculate internal force of cantilever beam
bool sfCtlInternalForce(int, double, double *);
//calculate internal force of displacement
bool sfDisplacementForce(int, double *);
//print"----------------------------------------"
bool sfPrintLine();
//print"****************************************"
bool sfPrintLine2();
//output data
bool sfOutput();
//print error
bool sfPrintError(int);
bool dovidw();

int main()
{
    double *ts = 0, *lv = 0; //declare total stiffness and load vector

    printf("Welcome to use the calculator of space frame!\nPress any key to start");
    char value = getchar(); //pause

    sfPrintLine(); //"------------------------------"
    if (sfInput()) //input data
    {
        sfPrintError(1);
        printf("\nPress any key to exit\n");
        value = getchar();

        return 1;
    }
    else
        printf("Data input succeeded!\n");

    int dof = 6 * NFRN;
    // ts = (double *)malloc(dof * dof * sizeof(double)); //allocate memory for total stiffness matrix
    // memset(ts, 0, dof * dof * sizeof(double));
    lv = (double *)malloc(dof * sizeof(double)); //allocate memory for load vector
    memset(lv, 0, dof * sizeof(double));

    if (sfBuildTotalStiff(ts)) //build total stiffness matrix
    {
        sfPrintError(2);
        printf("\nPress any key to exit\n");
        value = getchar();

        return 1;
    }
    else
        printf("Building total stiffness matrix succeeded!\n");

    if (sfBuildLoadVector(lv)) //build load stiffness vector
    {
        sfPrintError(3);
        printf("\nPress any key to exit\n");
        value = getchar();

        return 1;
    }
    else
        printf("Building load vector succeeded!\n");

    printf("nsi\t%d\n", nsi);
    printf("maxibdw\t%d\n", maxibdw);
    for (int i = 0; i < 6 * NFRN; i++)
    {
        printf("\t%d\n", iv[i]);
    }

    // if (sfCholesky(ts, lv, DON, 6 * NFRN)) //solve matrix equation
    // {
    //     sfPrintError(4);
    //     printf("\nPress any key to exit\n");
    //     value = getchar();

    //     return 1;
    // }
    // else
    //     printf("Solving equation succeeded!\n");

    // if (solve_conjugate_gradient(ts, lv, DON, 6 * NFRN)) //solve matrix equation
    // {
    //     sfPrintError(4);
    //     printf("\nPress any key to exit\n");
    //     value = getchar();

    //     return 1;
    // }
    // else
    //     printf("Solving equation succeeded!\n");

    // sfPrintLine2();
    // for (int i = 0; i < 6 * NFRN; i++)
    // {
    //     printf("%15.7f", DON[i]);
    //     printf("\n");
    // }
    // sfPrintLine2();

    for (int i = 0; i < NOS; i++)
        if (sfInternalForce(6 * i, NRS[i], DSB[i])) //calculate the internal force of each rods
        {
            sfPrintError(5);
            printf("\nPress any key to exit\n");
            value = getchar();

            return 1;
        }
    sfOutput(); //output data.

    sfPrintLine2();
    for (int i = 0; i < 6 * NOS; i++)
    {
        printf("\t%f\n", IFS[i]);
    }
    sfPrintLine2();

    printf("Press any key to exit\n");
    value = getchar(); //pause

    return 0;
}

bool sfInput()
{
    TNN = 7;
    NFIN = 1;
    NFRN = TNN - NFIN;
    NOR = 6;
    NOL = 5;
    NOS = 3;

    XCN = (double *)malloc(TNN * sizeof(double));
    memset(XCN, 0, TNN * sizeof(double));
    YCN = (double *)malloc(TNN * sizeof(double));
    memset(YCN, 0, TNN * sizeof(double));
    ZCN = (double *)malloc(TNN * sizeof(double));
    memset(ZCN, 0, TNN * sizeof(double));
    BNR = (int *)malloc(NOR * sizeof(int));
    memset(BNR, 0, NOR * sizeof(int));
    ENR = (int *)malloc(NOR * sizeof(int));
    memset(ENR, 0, NOR * sizeof(int));
    ELASTIC = (double *)malloc(NOR * sizeof(double));
    memset(ELASTIC, 0, NOR * sizeof(double));
    SHEAR = (double *)malloc(NOR * sizeof(double));
    memset(SHEAR, 0, NOR * sizeof(double));
    AREA = (double *)malloc(NOR * sizeof(double));
    memset(AREA, 0, NOR * sizeof(double));
    IMY = (double *)malloc(NOR * sizeof(double));
    memset(IMY, 0, NOR * sizeof(double));
    IMZ = (double *)malloc(NOR * sizeof(double));
    memset(IMY, 0, NOR * sizeof(double));
    THETA = (double *)malloc(NOR * sizeof(double));
    memset(THETA, 0, NOR * sizeof(double));
    NRL = (int *)malloc(NOL * sizeof(int));
    memset(NRL, 0, NOL * sizeof(int));
    PLI = (int *)malloc(NOL * sizeof(int));
    memset(PLI, 0, NOL * sizeof(int));
    KOL = (int *)malloc(NOL * sizeof(int));
    memset(KOL, 0, NOL * sizeof(int));
    VOL = (double *)malloc(NOL * sizeof(double));
    memset(VOL, 0, NOL * sizeof(double));
    DLB = (double *)malloc(NOL * sizeof(double));
    memset(DLB, 0, NOL * sizeof(double));
    NRS = (int *)malloc(NOS * sizeof(int));
    memset(NRS, 0, NOS * sizeof(int));
    DSB = (double *)malloc(NOS * sizeof(double));
    memset(DSB, 0, NOS * sizeof(double));
    DON = (double *)malloc(6 * NFRN * sizeof(double));
    memset(DON, 0, 6 * NFRN * sizeof(double));
    IFS = (double *)malloc(3 * NOS * sizeof(double));
    memset(IFS, 0, 3 * NOS * sizeof(double));
    RFE = (double *)malloc(6 * NOR * sizeof(double));
    memset(RFE, 0, 6 * NOR * sizeof(double));

    XCN[0] = 1;
    XCN[1] = 2;
    XCN[2] = 3;
    XCN[3] = 4;
    XCN[4] = 5;
    XCN[5] = 6;
    XCN[6] = 7;

    YCN[0] = 0;
    YCN[1] = 0;
    YCN[2] = 0;
    YCN[3] = 0;
    YCN[4] = 0;
    YCN[5] = 0;
    YCN[6] = 0;

    ZCN[0] = 0;
    ZCN[1] = 0;
    ZCN[2] = 0;
    ZCN[3] = 0;
    ZCN[4] = 0;
    ZCN[5] = 0;
    ZCN[6] = 0;

    BNR[0] = 1;
    BNR[1] = 2;
    BNR[2] = 3;
    BNR[3] = 4;
    BNR[4] = 5;
    BNR[5] = 6;

    ENR[0] = 2;
    ENR[1] = 3;
    ENR[2] = 4;
    ENR[3] = 5;
    ENR[4] = 6;
    ENR[5] = 7;

    ELASTIC[0] = 210000000;
    ELASTIC[1] = 210000000;
    ELASTIC[2] = 210000000;
    ELASTIC[3] = 210000000;
    ELASTIC[4] = 210000000;
    ELASTIC[5] = 210000000;

    SHEAR[0] = 80769000;
    SHEAR[1] = 80769000;
    SHEAR[2] = 80769000;
    SHEAR[3] = 80769000;
    SHEAR[4] = 80769000;
    SHEAR[5] = 80769000;

    AREA[0] = 0.007854;
    AREA[1] = 0.007854;
    AREA[2] = 0.007854;
    AREA[3] = 0.007854;
    AREA[4] = 0.007854;
    AREA[5] = 0.007854;

    IMY[0] = 0.0000049087;
    IMY[1] = 0.0000049087;
    IMY[2] = 0.0000049087;
    IMY[3] = 0.0000049087;
    IMY[4] = 0.0000049087;
    IMY[5] = 0.0000049087;

    IMZ[0] = 0.0000049087;
    IMZ[1] = 0.0000049087;
    IMZ[2] = 0.0000049087;
    IMZ[3] = 0.0000049087;
    IMZ[4] = 0.0000049087;
    IMZ[5] = 0.0000049087;

    THETA[0] = 0;
    THETA[1] = 0;
    THETA[2] = 0;
    THETA[3] = 0;
    THETA[4] = 0;
    THETA[5] = 0;

    NRS[0] = 1;
    NRS[1] = 2;
    NRS[2] = 3;

    DSB[0] = 3.0;
    DSB[1] = 1.5;
    DSB[2] = 2.598;

    NRL[0] = 1;
    NRL[1] = 1;
    NRL[2] = 1;
    NRL[3] = 1;
    NRL[4] = 2;

    PLI[0] = 0;
    PLI[1] = 1;
    PLI[2] = 0;
    PLI[3] = 1;
    PLI[4] = 0;

    KOL[0] = 2;
    KOL[1] = 1;
    KOL[2] = 1;
    KOL[3] = 6;
    KOL[4] = 1;

    VOL[0] = -0.8;
    VOL[1] = 4.0;
    VOL[2] = -1.0;
    VOL[3] = -3.0;
    VOL[4] = 2.0;

    DLB[0] = 3.0;
    DLB[1] = 3.0;
    DLB[2] = 6.0;
    DLB[3] = 6.0;
    DLB[4] = 3.0;

    return 0;
}

bool sfBuildTotalStiff(double *ts) //ts is total stiffness matrix
{
    dovidw();

    double us[36] = {0};            //unit stiffness matrix
    int p[2] = {0}, dof = 6 * NFRN; //p is a temperary vector for i0j0, dof is the degree of freedom of nods

    ts = (double *)malloc(nsi * sizeof(double)); //allocate memory for total stiffness matrix
    memset(ts, 0, nsi * sizeof(double));

    LCS = (double *)malloc(4 * NOR * sizeof(double)); //allocate memory for rods' parameter
    memset(LCS, 0, 4 * NOR * sizeof(double));

    if (sfLCosSin()) //calculate the length, cosine and sine of all rods
    {
        sfPrintError(6);
        return 1;
    }

    for (int k = 0; k < NOR; k++)
    {
        p[0] = 6 * (BNR[k] - NFIN - 1); // tag: match the displacement with nods
        p[1] = 6 * (ENR[k] - NFIN - 1);

        for (int i = 0; i < 2; i++)
        {
            if (p[i] >= 0) //determine free node
            {
                if (sfBuildUnitStiff(k, i + 1, us)) //build unit stiffness matrix
                {
                    sfPrintError(7);
                    return 1;
                }
                for (int m = 0; m < 6; m++)
                    for (int n = 0; n <= m; n++)
                    {
                        ts[iv[(p[i] + m)] + (p[i] + n) - (p[i] + m) - 1] += us[m * 6 + n]; //superpose
                        //printf("%d, %f\n", iv[(p[i] + m)] + (p[i] + n) - (p[i] + m), ts[iv[(p[i] + m)] + (p[i] + n) - (p[i] + m)]);
                    }
            }
        }
        if (p[0] >= 0 && p[1] >= 0)
        {
            if (sfBuildUnitStiff(k, 3 + 1, us)) //build unit stiffness matrix
            {
                sfPrintError(7);
                return 1;
            }
            for (int m = 0; m < 6; m++)
            {
                for (int n = 0; n < 6; n++)
                {
                    ts[iv[(p[1] + m)] + (p[0] + n) - (p[1] + m) - 1] += us[m * 6 + n]; //superpose
                }
            }
        }
    }

    for (int mm = 0; mm < 6 * NFRN; mm++)
    {
        for (int nn = 0; nn < 6 * NFRN; nn++)
        {
            if (mm == nn)
            {
                printf("%.1f,", ts[iv[mm] - 1]);
            }
            else if (nn > mm)
            {
                if ((iv[nn] - nn + mm) > iv[nn - 1])
                {
                    printf("%.2f,", ts[iv[nn] - nn + mm - 1]);
                }
                else
                {
                    printf("%d,", 0);
                }
            }
            else if (mm > nn)
            {
                if ((iv[mm] - mm + nn) > iv[mm - 1])
                    printf("%.1f,", ts[iv[mm] - mm + nn - 1]);
                else
                    printf("%d,", 0);
            }
        }
        printf("\n");
    }
    // for (int mm = 0; mm < 6 * NFRN;mm++)
    // {
    //     printf("%f \n", ts[iv[mm]]);
    // }
    // printf("%f\n", ts[iv[5]]);
    return 0;
}

bool sfLCosSin()
{
    for (int k = 0; k < NOR; k++)
    {
        int i = BNR[k] - 1, j = ENR[k] - 1; //index of beginning and end nodes of rods
        LCS[1 * NOR + k] = XCN[j] - XCN[i];
        LCS[2 * NOR + k] = YCN[j] - YCN[i];
        LCS[3 * NOR + k] = ZCN[j] - ZCN[i];
        LCS[0 * NOR + k] = sqrt(LCS[1 * NOR + k] * LCS[1 * NOR + k] + LCS[2 * NOR + k] * LCS[2 * NOR + k] + LCS[3 * NOR + k] * LCS[3 * NOR + k]);
        if (LCS[0 * NOR + k] < EPS) //if the length of rod is too small, then return error
        {
            sfPrintError(8);
            return 1;
        }
        LCS[1 * NOR + k] = LCS[1 * NOR + k] / LCS[0 * NOR + k];
        LCS[2 * NOR + k] = LCS[2 * NOR + k] / LCS[0 * NOR + k];
        LCS[3 * NOR + k] = LCS[3 * NOR + k] / LCS[0 * NOR + k];
    }

    return 0;
}

bool sfBuildUnitStiff(int k, int flag, double *us) //k is the number of rods, flag is the index of matrix parts, us is the unit stiffness matrix
{
    double rd[36] = {0}, t[36] = {0}, c[36] = {0}, tmp = 0; //rd is local stiffness matrix, t is transpose matrix, c is a temperary matrix
    memset(us, 0, 36 * sizeof(double));

    if (sfBuildLocalStiff(k, flag, rd)) //build local stiffness matrix
    {
        sfPrintError(9);
        return 1;
    }
    if (sfBuildTrans(k, t)) //build transpose matrix
    {
        sfPrintError(10);
        return 1;
    }

    for (int i = 0; i < 6; i++) //transpose matrix times local stiffness matrix, store the result in c
    {
        for (int m = 0; m < 6; m++)
        {
            tmp = t[i * 6 + m];
            for (int j = 0; j < 6; j++)
            {
                c[i * 6 + j] += tmp * rd[m * 6 + j];
            }
        }
    }

    for (int i = 0; i < 6; i++) //c times the transposition of transpose matrix, store the result in unit stiff
    {
        for (int j = 0; j < 6; j++)
        {
            for (int m = 0; m < 6; m++)
            {
                us[i * 6 + j] += c[i * 6 + m] * t[j * 6 + m];
            }
        }
    }

    return 0;
}

bool sfBuildLocalStiff(int k, int flag, double *rd) //k is the number of rods, flag is the number of matrix
{
    double a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0, l = LCS[0 * NOR + k];

    a = ELASTIC[k] * AREA[k] / l;         //EA/1
    b = SHEAR[k] * (IMY[k] + IMZ[k]) / l; //GJ(p)/1
    c = 4 * ELASTIC[k] * IMY[k] / l;      //4EJ(y)/1
    d = c / 2 * 3 / l;                    //6EJ(z)/l/l
    e = 2 * d / l;                        //12EJ(y)/l/l/l
    f = 4 * ELASTIC[k] * IMZ[k] / l;      //4EJ(z)/l
    g = f / 2 * 3 / l;                    //6EJ(Z)/l/l
    h = 2 * g / l;                        //12EJ(z)/l/l/l

    switch (flag)
    {
    case 1: //k11
        rd[0 * 6 + 0] = a;
        rd[1 * 6 + 1] = h;
        rd[1 * 6 + 5] = rd[5 * 6 + 1] = g;
        rd[2 * 6 + 2] = e;
        rd[2 * 6 + 4] = rd[4 * 6 + 2] = -d;
        rd[3 * 6 + 3] = b;
        rd[4 * 6 + 4] = c;
        rd[5 * 6 + 5] = f;
        break;
    case 2: //k22
        rd[0 * 6 + 0] = a;
        rd[1 * 6 + 1] = h;
        rd[1 * 6 + 5] = rd[5 * 6 + 1] = -g;
        rd[2 * 6 + 2] = e;
        rd[2 * 6 + 4] = rd[4 * 6 + 2] = d;
        rd[3 * 6 + 3] = b;
        rd[4 * 6 + 4] = c;
        rd[5 * 6 + 5] = f;
        break;
    case 3: //k12
        rd[0 * 6 + 0] = -a;
        rd[1 * 6 + 1] = -h;
        rd[1 * 6 + 5] = g;
        rd[5 * 6 + 1] = -g;
        rd[2 * 6 + 2] = -e;
        rd[2 * 6 + 4] = -d;
        rd[4 * 6 + 2] = d;
        rd[3 * 6 + 3] = -b;
        rd[4 * 6 + 4] = c / 2;
        rd[5 * 6 + 5] = f / 2;
        break;
    case 4: //k21
        rd[0 * 6 + 0] = -a;
        rd[1 * 6 + 1] = -h;
        rd[1 * 6 + 5] = -g;
        rd[5 * 6 + 1] = g;
        rd[2 * 6 + 2] = -e;
        rd[2 * 6 + 4] = d;
        rd[4 * 6 + 2] = -d;
        rd[3 * 6 + 3] = -b;
        rd[4 * 6 + 4] = c / 2;
        rd[5 * 6 + 5] = f / 2;
        break;

    default:
        break;
    }

    return 0;
}

bool sfBuildTrans(int k, double *t) //k is the number of rods, t is transpose matrix
{
    double coa = 0, cob = 0, coc = 0, sic = 0, sit = 0, cot = 0, m = 0, n = 0; //co means cosine, si means sine, m and n is temperary variable

    memset(t, 0, 36 * sizeof(double));

    coa = LCS[1 * NOR + k]; //cosine alpha
    cob = LCS[2 * NOR + k]; //cosine beta
    coc = LCS[3 * NOR + k]; //cosine gama
    sit = sin(THETA[k]);    //sine theta
    cot = cos(THETA[k]);    //cosine theta

    if (fabs(coc - 1) < EPS) //vertical(z axis positive direction) rods' transpose matrix
    {
        t[2 * 6 + 0] = t[5 * 6 + 3] = 1;
        t[0 * 6 + 1] = t[3 * 6 + 4] = t[1 * 6 + 2] = t[4 * 6 + 5] = sit;
        t[1 * 6 + 1] = t[4 * 6 + 4] = cot;
        t[0 * 6 + 2] = t[3 * 6 + 5] = -cot;
    }
    else if (fabs(coc + 1) < EPS) //vertical(z axis negative direction) rods' transpose matrix
    {
        t[2 * 6 + 0] = t[5 * 6 + 3] = -1;
        t[0 * 6 + 1] = t[3 * 6 + 4] = sit;
        t[1 * 6 + 2] = t[4 * 6 + 5] = -sit;
        t[1 * 6 + 1] = t[4 * 6 + 4] = t[0 * 6 + 2] = t[3 * 6 + 5] = cot;
    }
    else
    {
        sic = sqrt(1 - coc * coc); //sine gama
        m = coa * coc;             //cosine alpha times cosine gama
        n = cob * coc;             //cosine beta times cosine gama

        t[0 * 6 + 0] = t[3 * 6 + 3] = coa;
        t[1 * 6 + 0] = t[4 * 6 + 3] = cob;
        t[2 * 6 + 0] = t[5 * 6 + 3] = coc;
        t[0 * 6 + 1] = t[3 * 6 + 4] = (cob * sit - m * cot) / sic;
        t[1 * 6 + 1] = t[4 * 6 + 4] = -(n * cot + coa * sit) / sic;
        t[2 * 6 + 1] = t[5 * 6 + 4] = cot * sic;
        t[0 * 2 + 2] = t[3 * 6 + 5] = (m * sit + cob * cot) / sic;
        t[1 * 6 + 2] = t[4 * 6 + 5] = (n * sit - coa * cot) / sic;
        t[2 * 6 + 2] = t[5 * 6 + 5] = -sit * sic;
    }

    return 0;
}

bool sfBuildLoadVector(double *lv) //lv is the load vector
{
    int rod = 0, p[2] = {0};          //rod is the number of rods, dof is the degree of freedom
    double rf[12] = {0}, t[36] = {0}; //rf is the reaction force matrix, t is the transpose matrix, p is a temperary vector for i0j0

    for (int i = 0; i < NOL; i++)
    {
        rod = NRL[i] - 1;                   //the number of rods with load
        memset(rf, 0, 12 * sizeof(double)); //zero clearing

        if (sfReactionForce(i, &rf[0 * 6], &rf[1 * 6])) //calculate reaction force
        {
            sfPrintError(11);
            return 1;
        }
        for (int j = 0; j < 6; j++) //add reaction force to RFE
        {
            RFE[6 * rod + j] += rf[1 * 6 + j];
        }
        if (sfBuildTrans(rod, t)) //build transpose matrix
        {
            sfPrintError(10);
            return 1;
        }

        p[0] = 6 * (BNR[rod] - NFIN - 1); // tag: match the displacement with nods
        p[1] = 6 * (ENR[rod] - NFIN - 1);

        for (int j = 0; j < 2; j++) //add reaction force to load vector
        {
            if (p[j] >= 0) //determine free node
                for (int m = 0; m < 6; m++)
                    for (int n = 0; n < 6; n++)
                        lv[p[j] + m] -= t[m * 6 + n] * rf[j * 6 + n];
        }
    }

    return 0;
}

bool sfReactionForce(int i, double *rfb, double *rfe) //i is the number of load, rfb and rfe is the reaction force at begining and end of rods
{
    double ra = 0, rb = 0, a = 0, b = 0, q = VOL[i], xq = DLB[i]; //ra, rb, a and b are middle variable
    int rod = NRL[i] - 1, pm = PLI[i], t = 0;                     //rod is the number of rods

    if (pm == 0) //load is in XY plane
    {
        t = -1; //The bending moment in the support-reaction equation is positive clockwise, convert it to positive to the coordinate axis
    }
    else if (pm == 1) //load is in XZ plane
    {
        t = 1; //The bending moment in the support-reaction equation is positive clockwise, convert it to positive to the coordinate axis
    }
    ra = DLB[i] / LCS[0 * NOR + rod]; //x(q) / L
    rb = 1 - ra;                      //1 - x(q) / L
    switch (KOL[i])
    {
    case 1: //vertical concentrating load
        a = rb * rb;
        rfb[pm + 1] = -q * rb * (1 + ra - 2 * ra * ra);
        rfe[pm + 1] = -q - rfb[pm + 1];
        rfb[5 - pm] = t * q * rb * ra * (LCS[0 * NOR + rod] - xq);
        rfe[5 - pm] = -t * q * ra * rb * xq;
        break;
    case 2: //vertical uniform load
        a = q * xq;
        b = a * xq / 12;
        rfb[pm + 1] = -a * (1 + 0.5 * ra * ra * ra - ra * ra);
        rfe[pm + 1] = -a - rfb[pm + 1];
        rfb[5 - pm] = t * b * (6 - 8 * ra + 3 * ra * ra);
        rfe[5 - pm] = -t * b * (4 * ra - 3 * ra * ra);
        break;
    case 3: //axial concentrating force when PLI == 0, torque when PLI ==1
        rfb[3 * pm] = -q * rb;
        rfe[3 * pm] = -q * ra;
        break;
    case 4: //axial uniform load
        a = q * xq;
        rfe[3 * pm] = -a * ra / 2;
        rfb[3 * pm] = -a - rfe[3 * pm];
        break;
    case 5: //vertical triangle distributed load
        a = q * xq / 2;
        b = -0.4 * ra * ra;
        rfb[pm + 1] = -2 * a * (0.5 - 0.75 * ra * ra + 0.4 * ra * ra * ra);
        rfe[pm + 1] = -a - rfb[pm + 1];
        rfb[5 - pm] = t * a * (2 / 3 + b - ra);
        rfe[5 - pm] = -t * a * (0.5 * ra + b);
        break;
    case 6: //concentrating bending moment
        rfb[2 - pm] = t * 6 * q * rb * ra / LCS[0 * NOR + rod];
        rfe[2 - pm] = -rfb[2 - pm];
        rfb[pm + 4] = t * q * rb * (-1 + 3 * ra);
        rfe[pm + 4] = t * q * ra * (2 - 3 * ra);
        break;
    case 7: //unifrom temperature rise
        rfb[0] = q * xq * ELASTIC[rod] * AREA[rod];
        rfe[0] = -rfb[0];
        break;
    case 8: //different temperature rise
        if (pm == 0)
        {
            a = IMZ[rod];
        }
        else if (pm == 1)
        {
            a = IMY[rod];
        }
        rfb[5 - pm] = t * q * 2 * ELASTIC[rod] * a * xq;
        rfe[5 - pm] = -rfb[5 - pm];
        break;
    default:
        break;
    }

    return 0;
}

bool sfCholesky(double *A, double *b, double *x, int n) //Ax=b, n=size(A)
{
    if (A == NULL)
    {
        sfPrintError(12);
        return 1;
    }
    else if (b == NULL)
    {
        sfPrintError(12);
        return 1;
    }
    else if (x == NULL)
    {
        sfPrintError(12);
        return 1;
    }
    else if (n == 0)
    {
        sfPrintError(12);
        return 1;
    }

    double sum = 0, *L = 0, *D = 0, *y = 0;

    L = (double *)malloc(n * n * sizeof(double));
    if (L != NULL)
        memset(L, 0, n * n * sizeof(double));
    for (int i = 0; i < n; i++)
    {
        L[i * n + i] = 1;
    }
    D = (double *)malloc(n * sizeof(double));
    memset(D, 0, n * sizeof(double));
    y = (double *)malloc(n * sizeof(double));
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

bool solve_conjugate_gradient(double *A, double *b, double *x, int N)
{
    double *r, *p, *a, *z;
    double gamma, gamma_new, alpha, beta;

    r = (double *)malloc(N * sizeof(double));
    memset(r, 0, sizeof(double));
    p = (double *)malloc(N * sizeof(double));
    memset(p, 0, sizeof(double));
    z = (double *)malloc(N * sizeof(double));
    memset(z, 0, sizeof(double));

    // x = [0 ... 0]
    // r = b - A * x
    // p = r
    // gamma = r' * r
    gamma = 0.0;
    for (int i = 0; i < N; ++i)
    {
        x[i] = 0.0;
        r[i] = b[i];
        p[i] = r[i];
        gamma += r[i] * r[i];
    }

    for (int n = 0; 1; ++n)
    {
        // z = A * p
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                if (i == j)
                {
                    a[6 * i + j] = A[iv[i] - 1];
                }
                else if (j > i)
                {
                    if ((iv[j] - j + i) > iv[j - 1])
                    {
                        a[i * 6 + j] = A[iv[j] - j + i - 1];
                    }
                    else
                    {
                        a[i * 6 + j] = 0;
                    }
                }
                else if (i > j)
                {
                    if ((iv[i] - i + j) > iv[i - 1])
                        a[i * 6 + j] = A[iv[i] - i + j - 1];
                    else
                        a[i * 6 + j] = 0;
                }
                z[i] += a[i*6+j] * p[j];
            }
        }

        // alpha = gamma / (p' * z)
        alpha = 0.0;
        for (int i = 0; i < N; ++i)
            alpha += p[i] * z[i];
        alpha = gamma / alpha;

        // x = x + alpha * p
        // r = r - alpha * z
        // gamma_new = r' * r
        gamma_new = 0.0;
        for (int i = 0; i < N; ++i)
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * z[i];
            gamma_new += r[i] * r[i];
        }

        if (sqrt(gamma_new) < EPS)
            break;

        beta = gamma_new / gamma;

        // p = r + (gamma_new / gamma) * p;
        for (int i = 0; i < N; ++i)
            p[i] = r[i] + beta * p[i];

        // gamma = gamma_new
        gamma = gamma_new;
    }

    free(r);
    free(p);
    free(z);

    return 0;
}

bool sfInternalForce(int m, int k, double xp) //m is the number of sections, k is the actual number of rods, xp is the distance between the section and the begining of rods
{
    int n = 6 * (k - 1); //n is the matching place of rods
    double tf[6] = {0};  //tf is temperary variable

    IFS[m] = RFE[n]; //calculate internal force cause by reaction force at the end of rods
    IFS[m + 1] = -RFE[n + 1];
    IFS[m + 2] = -RFE[n + 2];
    IFS[m + 3] = RFE[n + 3];
    IFS[m + 4] = -RFE[n + 4] + RFE[n + 2] * (LCS[0 * NOR + k - 1] - xp);
    IFS[m + 5] = RFE[n + 5] + RFE[n + 1] * (LCS[0 * NOR + k - 1] - xp);

    for (int i = 0; i < NOL; i++) //for every rods
    {
        if (NRL[i] == k) //if load is on rod k
        {
            memset(tf, 0, 6 * sizeof(double)); //zero clear tf
            if (sfCtlInternalForce(i, xp, tf)) // calculate internal force of cantilever beam
            {
                sfPrintError(13);
                return 1;
            }
            for (int j = 0; j < 6; j++) //add internal force of cantilever into IFR
            {
                IFS[m + j] += tf[j];
            }
        }
    }
    if (sfDisplacementForce(k, tf)) //calculate end force
    {
        sfPrintError(14);
        return 1;
    }

    IFS[m] -= tf[0]; //calculate section force cause by end force
    IFS[m + 1] += tf[1];
    IFS[m + 2] += tf[2];
    IFS[m + 3] -= tf[3];
    IFS[m + 4] += tf[4] + tf[2] * xp;
    IFS[m + 5] += -tf[5] + tf[1] * xp;

    return 0;
}

bool sfCtlInternalForce(int i, double xp, double *tf) //i is the number of load, xp is the distance between the section and the begining of rod, tf is internal force
{
    double xq = DLB[i], t = xq - xp, r = xp / xq, q = VOL[i]; //t and r are temperary variables
    int e = PLI[i];
    switch (KOL[i]) //calculate section force according to kind of loads
    {
    case 1:
        if (xp < xq)
        {
            tf[e + 1] = -q;
            tf[5 - e] = q * t;
        }
        break;
    case 2:
        if (xp < xq)
        {
            tf[e + 1] = -q * t;
            tf[5 - e] = 0.5 * q * t * t;
        }
        break;
    case 3:
        if (xp < xq)
        {
            tf[3 * e] = q;
        }
        break;
    case 4:
        if (xp < xq)
        {
            tf[3 * e] = q * t;
        }
        break;
    case 5:
        if (xp < xq)
        {
            tf[e + 1] = -q * (1 + r) * t / 2;
            tf[5 - e] = q * t * t * (2 + r) / 6;
        }
        break;
    case 6:
        if (xp < xq)
        {
            tf[e + 4] = (2 * e - 1) * q;
        }
        break;
    case 7: //temperature change don't generate internal force on cantilever beam
        break;
    case 8:
        break;

    default:
        break;
    }

    return 0;
}

bool sfDisplacementForce(int k, double *tref) //k is the actual number of rods, tref is the end force of rods
{
    int p[2] = {0};                                  //p is a temperary vector for i0j0
    double rd[36] = {0}, rdb[36] = {0}, t[36] = {0}; //rd

    memset(tref, 0, 6 * sizeof(double));

    if (sfBuildTrans(k - 1, t)) //calculate transpose matrix
    {
        sfPrintError(10);
        return 1;
    }

    p[0] = 6 * (BNR[k - 1] - NFIN - 1); // tag: match the displacement with nods
    p[1] = 6 * (ENR[k - 1] - NFIN - 1);

    for (int i = 0; i < 2; i++)
    {
        if (p[i] >= 0) //determine free node
        {
            if (sfBuildLocalStiff(k - 1, 2 * i + 1, rd)) //build unit stiffness matrix
            {
                sfPrintError(9);
                return 1;
            }

            memset(rdb, 0, 36 * sizeof(double)); //zero clean rdb

            for (int j = 0; j < 6; j++) //rd times transposition of transpose matrix
                for (int m = 0; m < 6; m++)
                    for (int n = 0; n < 6; n++)
                        rdb[j * 6 + m] += rd[j * 6 + n] * t[m * 6 + n];
            for (int j = 0; j < 6; j++) //rdb times DON
                for (int m = 0; m < 6; m++)
                    tref[j] += rdb[j * 6 + m] * DON[p[i] + m];
        }
        else //fixed node
        {
            for (int j = 0; j < 3; j++)
            {
                tref[j] += 0;
            }
        }
    }

    return 0;
}

bool sfPrintLine()
{
    printf("--------------------------------------------------------------------------\n");
    return 0;
}

bool sfPrintLine2()
{
    printf("**************************************************************************\n");
    return 0;
}

bool sfOutput()
{
    return 0;
}

bool sfPrintError(int error)
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
        printf("The length of a rod is too small!\n");
        break;
    case 9:
        printf("Building local stiffness matrix filed!\n");
        break;
    case 10:
        printf("Building transpose matrix failed!\n");
        break;
    case 11:
        printf("Calculating reaction force failed!\n");
        break;
    case 12:
        printf("There is something wrong in the equation!\n");
        break;
    case 13:
        printf("calculating internal force of cantilever beam failed!\n");
        break;
    case 14:
        printf("Calculating end force failed!\n");
        break;
    default:
        break;
    }
    printf("There is at least one error in your file, please check it and try it one more time.\n");

    return 0;
}
bool dovidw()
{
    int it = 0, mm = 0, i, j;
    maxibdw = 0;
    peribdw = (int *)malloc(TNN * sizeof(int));
    memset(peribdw, 0, NFRN * sizeof(int));
    iv = (int *)malloc(6 * NFRN * sizeof(int));
    memset(iv, 0, 6 * NFRN * sizeof(int));
    for (i = 0; i < NOR; i++)
    {
        if (BNR[i] > NFIN)
        {
            mm = ENR[i] - BNR[i];
            if (mm > peribdw[ENR[i] - 1])
            {
                peribdw[ENR[i] - 1] = mm;
                printf("%d %d\n", ENR[i], peribdw[ENR[i] - 1]);
            }
        }
    }
    for (i = 0; i < TNN; i++)
    {
        printf("%d %d\n", i, peribdw[i]);
    }
    for (i = NFIN; i < TNN; i++)
    {
        if (peribdw[i - NFIN + 1] > maxibdw)
        {
            maxibdw = peribdw[i - NFIN + 1];
        }
        for (j = 1; j <= 6; j++)
        {
            it = it + 1;
            if (it == 1)
            {
                iv[it - 1] = 6 * peribdw[i - NFIN + 1] + j;
            }
            else
            {
                iv[it - 1] = iv[it - 2] + 6 * peribdw[i - NFIN + 1] + j;
            }
        }
    }
    maxibdw = 6 * maxibdw + 5;
    nsi = iv[6 * NFRN - 1];
}
bool VibmRead(double *a, int mm, int nn, double *ts)
{
    if (mm == nn)
    {
        a[6 * mm + nn] = ts[iv[mm] - 1];
    }
    else if (nn > mm)
    {
        if ((iv[nn] - nn + mm) > iv[nn - 1])
        {
            a[mm * 6 + nn] = ts[iv[nn] - nn + mm - 1];
        }
        else
        {
            a[mm * 6 + nn] = 0;
        }
    }
    else if (mm > nn)
    {
        if ((iv[mm] - mm + nn) > iv[mm - 1])
            a[mm * nn] = ts[iv[mm] - mm + nn - 1];
        else
            a[mm + nn] = 0;
    }
    return 0;
}
