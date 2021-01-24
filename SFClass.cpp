#include <Windows.h>
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

class SpaceFrame
{
private:
    double EPS;
    double MAXTS;
    double MAXLV;

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
    double *SIGMA_1;
    double *SIGMA_2;
    int *DANGER;

    int *IV;     //the location of diagonal element
    int NSI;     //upper limit
    int MAXIBDW; //half bandwidth

    double *TS; //total stiffness

    int MAXTNN;   //the max total number of nodes with specify displacement
    int TNNSD[6]; //the total number of nodes with specify displacement.
                  //for example 4,0,0,0,1,0 means that there are four nodes with delta x and one node with theta y
    int *NNSD;    //the number of nodes with specify displacement
    double *VSD;  //the value of specify displacement

    double *ROU; //the density of rods
    double g;    //acceleration of gravity

    //build total stiffness matrix
    bool sfBuildTotalStiff() //ts is total stiffness matrix
    {
        dovidw();

        double us[36] = {0}; //unit stiffness matrix
        int p[2] = {0};      //p is a temperary vector for i0j0, dof is the degree of freedom of nods

        TS = (double *)malloc(NSI * sizeof(double)); //allocate memory for total stiffness matrix
        memset(TS, 0, NSI * sizeof(double));
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
                            TS[IV[(p[i] + m)] + (p[i] + n) - (p[i] + m) - 1] += us[m * 6 + n]; //superpose
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
                        TS[IV[(p[1] + m)] + (p[0] + n) - (p[1] + m) - 1] += us[m * 6 + n]; //superpose
                    }
                }
            }
        }

        for (int i = 0; i < NSI; i++)
        {
            // printf("%f\n", TS[i]);
            if (fabs(TS[i]) > MAXTS)
                MAXTS = TS[i];
        }

        // sfPrintLine2();
        // for (int i = 0; i < dof; i++)
        // {
        //     for (int j = 0; j < dof; j++)
        //     {
        //         printf("%15.2f", ts[i * dof + j]);
        //     }
        //     printf("\n");
        // }
        // sfPrintLine2();

        return 0;
    }
    //prepare to build variable bandwidth matrix
    bool dovidw()
    {
        int it = 0, mm = 0, i, j;
        int *peribdw; //bandwidth per line
        MAXIBDW = 0;
        peribdw = (int *)malloc(TNN * sizeof(int));
        memset(peribdw, 0, TNN * sizeof(int));
        IV = (int *)malloc(6 * NFRN * sizeof(int));
        memset(IV, 0, 6 * NFRN * sizeof(int));
        for (i = 0; i < NOR; i++)
        {
            if (BNR[i] > NFIN)
            {
                mm = ENR[i] - BNR[i];
                if (mm > peribdw[ENR[i] - 1])
                {
                    peribdw[ENR[i] - 1] = mm;
                    // printf("%d %d\n", ENR[i], peribdw[ENR[i] - 1]);
                }
            }
        }
        // for (i = 0; i < TNN; i++)
        // {
        //     printf("%d %d\n", i, peribdw[i]);
        // }
        for (i = NFIN; i < TNN; i++)
        {
            if (peribdw[i] > MAXIBDW)
            {
                MAXIBDW = peribdw[i];
            }
            for (j = 1; j <= 6; j++)
            {
                it = it + 1;
                if (it == 1)
                {
                    IV[it - 1] = 6 * peribdw[i] + j;
                }
                else
                {
                    IV[it - 1] = IV[it - 2] + 6 * peribdw[i] + j;
                }
            }
        }
        MAXIBDW = 6 * MAXIBDW + 5;
        NSI = IV[6 * NFRN - 1];
        free(peribdw);
        return 0;
    }
    //calculate the length sine and cosine of rods
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
    //build unit stiffness matrix
    bool sfBuildUnitStiff(int k, int flag, double *us) //k is the number of rods, flag is the index of matrix parts, us is the unit stiffness matrix
    {
        if (k < 0)
        {
            sfPrintError(16);
            return 0;
        }
        if (flag < 1 || flag > 4)
        {
            sfPrintError(16);
            return 0;
        }
        if (us == NULL)
        {
            sfPrintError(16);
            return 0;
        }

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
    //build local stiffness matrix
    bool sfBuildLocalStiff(int k, int flag, double *rd) //k is the number of rods, flag is the number of matrix
    {
        if (k < 0)
        {
            sfPrintError(17);
            return 0;
        }
        if (flag < 0 || flag > 4)
        {
            sfPrintError(17);
            return 0;
        }
        if (rd == NULL)
        {
            sfPrintError(17);
            return 0;
        }

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
    //build transpose matrix
    bool sfBuildTrans(int k, double *t) //k is the number of rods, t is transpose matrix
    {
        if (k < 0)
        {
            sfPrintError(18);
            return 0;
        }
        if (t == NULL)
        {
            sfPrintError(18);
            return 0;
        }

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
    //build load vector
    bool sfBuildLoadVector(double *lv) //lv is the load vector
    {
        if (lv == 0)
        {
            sfPrintError(19);
            return 0;
        }

        int rod = 0, p[2] = {0}, IJ = 0;  //rod is the number of rods, dof is the degree of freedom
        double rf[12] = {0}, t[36] = {0}; //rf is the reaction force matrix, t is the transpose matrix, p is a temperary vector for i0j0

        for (int i = 0; i < NOR; i++)
        {
            NRL[NOL + i] = i + 1;
            KOL[NOL + i] = 2;
            DLB[NOL + i] = LCS[0 * NOR + i];
            VOL[NOL + i] = -1 * AREA[i] * g * ROU[i] * sqrt(1 - LCS[3 * NOR + i] * LCS[3 * NOR + i]);
            PLI[NOL + i] = 0;
        }
        for (int i = 0; i < NOR; i++)
        {
            NRL[NOL + NOR + i] = i + 1;
            KOL[NOL + NOR + i] = 4;
            DLB[NOL + NOR + i] = LCS[0 * NOR + i];
            VOL[NOL + NOR + i] = -1 * AREA[i] * g * ROU[i] * LCS[3 * NOR + i];
            PLI[NOL + NOR + i] = 0;
        }

        printf("%f\n", LCS[3]);

        // for (int i = 0; i < NOL + 2 * NOR; i++)
        // {
        //     printf("%f\n", VOL[i]);
        // }

        for (int i = 0; i < (NOL + 2 * NOR); i++)
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
        for (int i = 0; i < 6; i++)
        {
            if (TNNSD[i] <= 0)
            {
                continue;
            }
            for (int j = 0; j < MAXTNN; j++)
            {
                if (NNSD[i * MAXTNN + j] == 0)
                {
                    continue;
                }
                IJ = 6 * (NNSD[i * MAXTNN + j] - NFIN) - 6;
                if (IJ >= 0)
                {

                    if (VSD[i * MAXTNN + j] == 0)
                    {
                        TS[IV[IJ + i] - 1] += 10000000000;
                    }
                    else
                    {
                        TS[IV[IJ + i] - 1] = 10000000000;
                        lv[IJ + i] = 10000000000 * VSD[i * MAXTNN + j];
                    }
                }
            }
        }

        // sfPrintLine2();
        // for (int i = 0; i < 6 * NFRN; i++)
        // {
        //     printf("%20.7f,", lv[i]);
        //     printf("\n");
        // }
        // sfPrintLine2();

        for (int i = 0; i < 6 * NFRN; i++)
            if (fabs(lv[i]) > MAXLV)
                MAXLV = lv[i];

        return 0;
    }
    //calculate reaction force
    bool sfReactionForce(int i, double *rfb, double *rfe) //i is the number of load, rfb and rfe is the reaction force at begining and end of rods
    {
        if (i < 0)
        {
            sfPrintError(20);
            return 0;
        }
        if (rfb == NULL)
        {
            sfPrintError(20);
            return 0;
        }
        if (rfe == NULL)
        {
            sfPrintError(20);
            return 0;
        }

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
    //solve equation of matrix by cholesky
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
    //solve equation of matrix by conjugate gradient
    bool solve_conjugate_gradient(double *A, double *b, double *x, int N)
    {
        double *r, *p, *z;
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
            for (int i = 0; i < N; i++)
            {
                z[i] = 0.0;
                for (int j = 0; j < N; j++)
                {
                    if (i == j)
                    {
                        z[i] += A[IV[i] - 1] * p[j];
                    }
                    else if (j > i)
                    {
                        if ((IV[j] - j + i) > IV[j - 1])
                        {
                            z[i] += A[IV[j] - j + i - 1] * p[j];
                        }
                        else
                        {
                            z[i] += 0;
                        }
                    }
                    else if (i > j)
                    {
                        if ((IV[i] - i + j) > IV[i - 1])
                            z[i] += A[IV[i] - i + j - 1] * p[j];
                        else
                            z[i] += 0;
                    }
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
    //solve equation of matrix by conjugate gradient parallel
    bool solve_conjugate_gradient_par(double *A, double *b, double *x, int N)
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
        else if (N == 0)
        {
            sfPrintError(12);
            return 1;
        }

        double *r, *p, *z;
        double gamma, gamma_new, alpha, beta;

        r = (double *)malloc(N * sizeof(double));
        memset(r, 0, sizeof(double));
        p = (double *)malloc(N * sizeof(double));
        memset(p, 0, sizeof(double));
        z = (double *)malloc(N * sizeof(double));
        memset(z, 0, sizeof(double));

        for (int i = 0; i < NSI; i++)
        {
            A[i] = A[i] / MAXTS;
        }
        for (int i = 0; i < N; i++)
        {
            b[i] = b[i] / MAXLV;
        }

        // x = [0 ... 0]
        // r = b - A * x
        // p = r
        // gamma = r' * r
        gamma = 0.0;
#pragma omp parallel for reduction(+ \
                                   : gamma)
        for (int i = 0; i < N; ++i)
        {
            x[i] = 0.0;
            r[i] = b[i];
            p[i] = r[i];
            gamma += r[i] * r[i];
        }

        for (int n = 0; true; ++n)
        {
// z = A * p
#pragma omp parallel for
            for (int i = 0; i < N; i++)
            {
                z[i] = 0.0;
                for (int j = 0; j < N; j++)
                {
                    if (i == j)
                    {
                        z[i] += A[IV[i] - 1] * p[j];
                    }
                    else if (j > i)
                    {
                        if ((IV[j] - j + i) > IV[j - 1])
                        {
                            z[i] += A[IV[j] - j + i - 1] * p[j];
                        }
                        else
                        {
                            z[i] += 0;
                        }
                    }
                    else if (i > j)
                    {
                        if ((IV[i] - i + j) > IV[i - 1])
                            z[i] += A[IV[i] - i + j - 1] * p[j];
                        else
                            z[i] += 0;
                    }
                }
            }

            // alpha = gamma / (p' * z)
            alpha = 0.0;
#pragma omp parallel for reduction(+ \
                                   : alpha)
            for (int i = 0; i < N; ++i)
            {
                // printf("%d\n", i);
                alpha += p[i] * z[i];
            }
            alpha = gamma / alpha;

            // x = x + alpha * p
            // r = r - alpha * z
            // gamma_new = r' * r
            gamma_new = 0.0;
#pragma omp parallel for reduction(+ \
                                   : gamma_new)
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
#pragma omp parallel for
            for (int i = 0; i < N; ++i)
            {
                p[i] = r[i] + beta * p[i];
            }

            // gamma = gamma_new
            gamma = gamma_new;
        }

        for (int i = 0; i < N; i++)
        {
            x[i] = x[i] * MAXLV / MAXTS;
        }

        free(r);
        free(p);
        free(z);

        return 0;
    }
    //calculate internal force of rods
    bool sfInternalForce(int m, int k, double xp) //m is the number of sections, k is the actual number of rods, xp is the distance between the section and the begining of rods
    {
        if (m < 0)
        {
            sfPrintError(21);
            return 0;
        }
        if (k < 0)
        {
            sfPrintError(21);
            return 0;
        }

        int n = 6 * (k - 1); //n is the matching place of rods
        double tf[6] = {0};  //tf is temperary variable

        IFS[m] = RFE[n]; //calculate internal force cause by reaction force at the end of rods
        IFS[m + 1] = -RFE[n + 1];
        IFS[m + 2] = -RFE[n + 2];
        IFS[m + 3] = RFE[n + 3];
        IFS[m + 4] = -RFE[n + 4] + RFE[n + 2] * (LCS[0 * NOR + k - 1] - xp);
        IFS[m + 5] = RFE[n + 5] + RFE[n + 1] * (LCS[0 * NOR + k - 1] - xp);

        for (int i = 0; i < (NOL + 2 * NOR); i++) //for every rods
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

        double moment = sqrt(IFS[m + 4] * IFS[m + 4] + IFS[m + 5] * IFS[m + 5]);
        double sigma_1 = moment * sqrt(AREA[k - 1] / 3.1415926) / IMZ[k - 1] + IFS[m] / ELASTIC[k - 1] / AREA[k - 1];
        double sigma_2 = -moment * sqrt(AREA[k - 1] / 3.1415926) / IMZ[k - 1] + IFS[m] / ELASTIC[k - 1] / AREA[k - 1];

        if (sigma_1 > SIGMA_1[m / 6] || sigma_2 < SIGMA_2[m / 6])
            DANGER[m / 6] = 1;

        return 0;
    }
    //calculate internal force of cantilever beam
    bool sfCtlInternalForce(int i, double xp, double *tf) //i is the number of load, xp is the distance between the section and the begining of rod, tf is internal force
    {
        if (i < 0)
        {
            sfPrintError(22);
            return 0;
        }
        if (tf == NULL)
        {
            sfPrintError(22);
            return 0;
        }

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
    //calculate internal force of displacement
    bool sfDisplacementForce(int k, double *tref) //k is the actual number of rods, tref is the end force of rods
    {
        if (k < 1)
        {
            sfPrintError(23);
            return 0;
        }
        if (tref == NULL)
        {
            sfPrintError(23);
            return 0;
        }

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
    //print"----------------------------------------"
    bool sfPrintLine()
    {
        printf("--------------------------------------------------------------------------\n");
        return 0;
    }
    //print"****************************************"
    bool sfPrintLine2()
    {
        printf("**************************************************************************\n");
        return 0;
    }
    //print error
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
        case 15:
            printf("Allocating total stiffness matrix failed!\n");
            break;
        case 16:
            printf("There is something wrong in building unit stiffness matrix!\n");
            break;
        case 17:
            printf("There is something wrong in building local stiffness matrix!\n");
            break;
        case 18:
            printf("There is something wrong in building transpose matrix failed!\n");
            break;
        case 19:
            printf("There is something wrong in building load vector!\n");
            break;
        case 20:
            printf("There is something wrong in calculating reaction force!\n");
            break;
        case 21:
            printf("There is something wrong in calculating internal force!\n");
            break;
        case 22:
            printf("There is something wrong in calculating internal force of cantilever!\n");
            break;
        case 23:
            printf("There is something wrong in calculating internal force of displacement!\n");
            break;
        case 24:
            printf("!\n");
            break;
        case 25:
            printf("!\n");
            break;

        default:
            break;
        }
        printf("There is at least one error in your file, please check it and try it one more time.\n");

        return 0;
    }

public:
    SpaceFrame();
    ~SpaceFrame();

    //read data from .csv
    bool sfInput();
    //calculate
    bool sfCalculate();
    //output data
    bool sfOutput();
};

SpaceFrame::SpaceFrame()
{
    EPS = 1e-15;
    MAXTS = 0;
    MAXLV = 0;

    TNN = 0;  //total number of nodes
    NFIN = 0; //number of fixed nodes
    NFRN = 0; //number of free nodes
    NOR = 0;  //number of rods
    NOL = 0;  //number of loads
    NOS = 0;  //number of sections

    XCN = NULL; //X coordinate of nodes
    YCN = NULL; //Y coordinate of nodes
    ZCN = NULL; //Z coordinate of nodes

    BNR = NULL;        //the beginning node number of rods
    ENR = NULL;        //the end node number of rods
    ELASTIC = NULL; //elastic modulus
    SHEAR = NULL;   //shear modulus
    AREA = NULL;    //area
    IMY = NULL;     //inertia moment of Y axis
    IMZ = NULL;     //inertia moment of Z axis
    THETA = NULL;   //theta the deflection angle of main inertia axis

    NRL = NULL;    //the number of rods with load
    PLI = NULL;    //the plane of the load's in
    KOL = NULL;    //the kind of load
    VOL = NULL; //the value of load
    DLB = NULL; //the distance between load and the beginning node

    NRS = NULL;    //the number of rod with section
    DSB = NULL; //the distance between section and the beginning node

    LCS = NULL; //the length, sine and cosine of rods
    DON = NULL; //the displacement of nodes
    IFS = NULL; //the internal force in the section
    RFE = NULL; //the reaction force of the end node
    SIGMA_1 = NULL;
    SIGMA_2 = NULL;
    DANGER = NULL;

    IV = NULL;  //the location of diagonal element
    NSI = 0;     //upper limit
    MAXIBDW = 0; //half bandwidth

    TS = NULL; //total stiffness

    MAXTNN = 0;     //the max total number of nodes with specify displacement
    TNNSD[6] = {0}; //the total number of nodes with specify displacement.
                        //for example 4,0,0,0,1,0 means that there are four nodes with delta x and one node with theta y
    NNSD = NULL;   //the number of nodes with specify displacement
    VSD = NULL; //the value of specify displacement

    ROU = NULL; //the density of rods
    g = 9.81;    //acceleration of gravity
}

SpaceFrame::~SpaceFrame()
{
    free(XCN);
    free(YCN);
    free(ZCN);
    free(BNR);
    free(ENR);
    free(ELASTIC);
    free(SHEAR);
    free(AREA);
    free(IMY);
    free(IMZ);
    free(THETA);
    free(NRL);
    free(PLI);
    free(KOL);
    free(VOL);
    free(DLB);
    free(NRS);
    free(DSB);
    free(DON);
    free(IFS);
    free(RFE);
}

bool SpaceFrame::sfInput()
{
    FILE *fp = NULL;                   //Define the file point
    char *line = 0, *data = 0;         //Define the line string and separated string
    char temporSpace[1000000];         //Apply for temporary storage space
    int rowIndex = 0, columnIndex = 0; //Reset the number of rows to zero, reset the number of columns to zero
    const char DIVIDE[] = ",";         //Set the separater as a ','

    if ((fp = fopen("source&result/sf_2.csv", "r")) == NULL) //Start the process when the file opens successfully
    {
        return 0;
    }

    fseek(fp, 0L, SEEK_SET);                                             //Locate file point to the first line
    while ((line = fgets(temporSpace, sizeof(temporSpace), fp)) != NULL) //The loop continues when the end of the file is not read
    {
        data = strtok(line, DIVIDE); //Split strings with a ',' as a separator
        while (data != NULL)         //Read the data of each row
        {
            if (strcmp(data, "END") == 0) //When the keyword 'END' is read, the reading process will be shut down
            {
                fclose(fp); //Close the file
                fp = NULL;  //Reset the file point
                return 0;
            }

            if (columnIndex++ == 0) //Skip the saving of the first column
            {
                data = strtok(NULL, DIVIDE); //Reset data
                continue;
            }

            switch (rowIndex) //Store variables of each column in different ways
            {
            case 0:
                break;
            case 1:
                if (columnIndex == 2)
                    TNN = atoi(data);
                break;
            case 2:
                if (columnIndex == 2)
                    NFIN = atoi(data);
                NFRN = TNN - NFIN;
                break;
            case 3:
                if (columnIndex == 2)
                    NOR = atoi(data);
                break;
            case 4:
                if (columnIndex == 2)
                    NOL = atoi(data);
                break;
            case 5:
                if (columnIndex == 2)
                {
                    NOS = atoi(data);
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
                    ROU = (double *)malloc(NOR * sizeof(double));
                    memset(ROU, 0.0, NOR * sizeof(double));
                    IMY = (double *)malloc(NOR * sizeof(double));
                    memset(IMY, 0, NOR * sizeof(double));
                    IMZ = (double *)malloc(NOR * sizeof(double));
                    memset(IMY, 0, NOR * sizeof(double));
                    THETA = (double *)malloc(NOR * sizeof(double));
                    memset(THETA, 0, NOR * sizeof(double));
                    NRL = (int *)malloc((NOL + 2 * NOR) * sizeof(int));
                    memset(NRL, 0, (NOL + 2 * NOR) * sizeof(int));
                    PLI = (int *)malloc((NOL + 2 * NOR) * sizeof(int));
                    memset(PLI, 0, (NOL + 2 * NOR) * sizeof(int));
                    KOL = (int *)malloc((NOL + 2 * NOR) * sizeof(int));
                    memset(KOL, 0, (NOL + 2 * NOR) * sizeof(int));
                    VOL = (double *)malloc((NOL + 2 * NOR) * sizeof(double));
                    memset(VOL, 0, (NOL + 2 * NOR) * sizeof(double));
                    DLB = (double *)malloc((NOL + 2 * NOR) * sizeof(double));
                    memset(DLB, 0, (NOL + 2 * NOR) * sizeof(double));
                    NRS = (int *)malloc(NOS * sizeof(int));
                    memset(NRS, 0, NOS * sizeof(int));
                    DSB = (double *)malloc(NOS * sizeof(double));
                    memset(DSB, 0, NOS * sizeof(double));
                    SIGMA_1 = (double *)malloc(NOS * sizeof(double));
                    memset(SIGMA_1, 0, NOS * sizeof(double));
                    SIGMA_2 = (double *)malloc(NOS * sizeof(double));
                    memset(SIGMA_2, 0, NOS * sizeof(double));
                    DANGER = (int *)malloc(NOS * sizeof(double));
                    memset(DANGER, 0, NOS * sizeof(int));
                    DON = (double *)malloc(6 * NFRN * sizeof(double));
                    memset(DON, 0, 6 * NFRN * sizeof(double));
                    IFS = (double *)malloc(6 * NOS * sizeof(double));
                    memset(IFS, 0, 3 * NOS * sizeof(double));
                    RFE = (double *)malloc(6 * NOR * sizeof(double));
                    memset(RFE, 0, 6 * NOR * sizeof(double));
                }
                break;
            case 6:
                if (columnIndex - 2 < TNN)
                    XCN[columnIndex - 2] = atof(data);
                break;
            case 7:
                if (columnIndex - 2 < TNN)
                    YCN[columnIndex - 2] = atof(data);
                break;
            case 8:
                if (columnIndex - 2 < TNN)
                    ZCN[columnIndex - 2] = atof(data);
                break;
            case 9:
                if (columnIndex - 2 < NOR)
                    BNR[columnIndex - 2] = atoi(data);
                break;
            case 10:
                if (columnIndex - 2 < NOR)
                    ENR[columnIndex - 2] = atoi(data);
                break;
            case 11:
                if (columnIndex - 2 < NOR)
                    ELASTIC[columnIndex - 2] = atof(data);
                break;
            case 12:
                if (columnIndex - 2 < NOR)
                    SHEAR[columnIndex - 2] = atof(data);
                break;
            case 13:
                if (columnIndex - 2 < NOR)
                    AREA[columnIndex - 2] = atof(data);
                break;
            case 14:
                if (columnIndex - 2 < NOR)
                {
                    ROU[columnIndex - 2] = atof(data);
                    // printf("%f\n", ROU[columnIndex - 2]);
                }
                break;
            case 15:
                if (columnIndex - 2 < NOR)
                    IMY[columnIndex - 2] = atof(data);
                break;
            case 16:
                if (columnIndex - 2 < NOR)
                    IMZ[columnIndex - 2] = atof(data);
                break;
            case 17:
                if (columnIndex - 2 < NOR)
                    THETA[columnIndex - 2] = atof(data);
                break;
            case 18:
                if (columnIndex - 2 < NOL)
                    NRL[columnIndex - 2] = atoi(data);
                break;
            case 19:
                if (columnIndex - 2 < NOL)
                    PLI[columnIndex - 2] = atoi(data);
                break;
            case 20:
                if (columnIndex - 2 < NOL)
                    KOL[columnIndex - 2] = atoi(data);
                break;
            case 21:
                if (columnIndex - 2 < NOL)
                    VOL[columnIndex - 2] = atof(data);
                break;
            case 22:
                if (columnIndex - 2 < NOL)
                    DLB[columnIndex - 2] = atof(data);
                break;
            case 23:
                if (columnIndex - 2 < NOS)
                    NRS[columnIndex - 2] = atoi(data);
                break;
            case 24:
                if (columnIndex - 2 < NOS)
                    DSB[columnIndex - 2] = atof(data);
                break;
            case 25:
                if (columnIndex - 2 < NOS)
                    SIGMA_1[columnIndex - 2] = atof(data);
                break;
            case 26:
                if (columnIndex - 2 < NOS)
                    SIGMA_2[columnIndex - 2] = atof(data);
                break;
            case 27:
                if (columnIndex - 2 < 6)
                {
                    TNNSD[columnIndex - 2] = atoi(data);
                    // printf("%d\n", TNNSD[columnIndex - 2]);
                    if (TNNSD[columnIndex - 2] > MAXTNN)
                        MAXTNN = TNNSD[columnIndex - 2];
                }
                break;
            case 28:
                if (NNSD == NULL)
                {
                    NNSD = (int *)malloc(6 * MAXTNN * sizeof(int));
                    memset(NNSD, 0, 6 * MAXTNN * sizeof(int));
                }
                if (columnIndex - 2 < 6 * MAXTNN)
                    NNSD[columnIndex - 2] = atoi(data);
                break;
            case 29:
                if (VSD == NULL)
                {
                    VSD = (double *)malloc(6 * MAXTNN * sizeof(double));
                    memset(VSD, 0, 6 * MAXTNN * sizeof(double));
                }
                if (columnIndex - 2 < MAXTNN)
                    VSD[columnIndex - 2] = atof(data);
                break;
            case 30:
                if (columnIndex == 2)
                    EPS = atof(data);
                break;

                //input finished
            }
            data = strtok(NULL, DIVIDE); //Reset data
        }
        rowIndex++;      //RowIndex steps forward once
        columnIndex = 0; //Reset columnIndex
    }
    fclose(fp); //Close the file
    fp = NULL;  //Reset the file point

    return 0;
}

bool SpaceFrame::sfOutput()
{
    printf("\t\t\tCalculation Of Space Rigid Frame\n");
    sfPrintLine();

    printf("\t\tTNN = %d\t\t\tNFIN = %d\n\t\tNFRN = %d\t\tNOR = %d\n", TNN, NFIN, NFRN, NOR);
    printf("\t\tNOL = %d\t\t\tNOS = %d\n", NOL, NOS);
    sfPrintLine();

    printf("NUMBER OF NODES     Coordinate-X    Coordinate-Y    Coordinate-Z\n");
    for (int i = 0; i < TNN; i++)
        printf("%15d%15.7f%15.7f%15.7f\n", i + 1, XCN[i], YCN[i], ZCN[i]);
    sfPrintLine();

    printf("NUMBER OF NODES     LEFT NODES    RIGHT NODES  Elastic modulus  Shear modulus    Area   Inertia moment Y axis  Inertia moment Z axis\n");
    for (int i = 0; i < NOR; i++)
        printf("%15d%15d%15d%15.0f%15.0f%11.4f%16.5f%23.5f\n", i + 1, BNR[i], ENR[i], ELASTIC[i], SHEAR[i], AREA[i], IMY[i], IMZ[i]);
    sfPrintLine();
    printf("NUMBER OF SECTIONS         PLI            DLB\n");
    for (int i = 0; i < NOS; i++)
        printf("%15d%15d%15.7f\n", i + 1, NRS[i], DSB[i]);
    sfPrintLine();

    printf("Calculating......\nThe results are as follows: \n");
    sfPrintLine();

    printf("NUMBER OF NODES   Displacement-X Displacement-Y Displacement-Z   Diversion-X    Diversion-Y    Diversion-Z\n");
    for (int i = NFIN; i < TNN; i++)
        printf("%15d%15.7f%15.7f%15.7f%15.7f%15.7f%15.7f\n", i + 1, DON[6 * (i - NFIN)], DON[6 * (i - NFIN) + 1], DON[6 * (i - NFIN) + 2], DON[6 * (i - NFIN) + 3], DON[6 * (i - NFIN) + 4], DON[6 * (i - NFIN) + 5]);
    sfPrintLine();

    printf("NUMBER OF SECTIONS Axial force-X  Shear force-Y  Shear force-Z    Torque-X   Bending moment-Y  Bending moment-Z\n");
    for (int i = 0; i < NOS; i++)
        printf("%15d%15.7f%15.7f%15.7f%15.7f%15.7f%15.7f\n", i + 1, IFS[6 * i], IFS[6 * i + 1], IFS[6 * i + 2], IFS[6 * i + 3], IFS[6 * i + 4], IFS[6 * i + 5]);

    int cnt = 0;
    printf("DANGEROUS SECTIONS:");
    for (int i = 0; i < NOS; i++)
    {
        // printf("&&%d&&", DANGER[i]);
        if (DANGER[i])
        {
            printf("%15d", i + 1);
            cnt++;
        }
    }
    if (cnt)
    {
        printf("\n");
    }
    else
    {
        printf("NONE\n");
    }

    FILE *fp = NULL;
    fp = fopen("source&result/sfResultClass.csv", "w");
    fprintf(fp, "TITLE\n");
    fprintf(fp, "TNN,%d\nNFIN,%d\nNFRN,%d\nNOR,%d\nNOL,%d\nNOS,%d", TNN, NFIN, NFRN, NOR, NOL, NOS);

    //------------NODES-------------------------------------------------
    fprintf(fp, "\nNODES,");
    for (int i = 0; i < TNN; i++)
        fprintf(fp, "%d,", i + 1);

    fprintf(fp, "\nCON,");
    for (int i = 0; i < TNN; i++)
        fprintf(fp, "(%f %f %f),", XCN[i], YCN[i], ZCN[i]);

    //------------RODS-------------------------------------------------
    fprintf(fp, "\nRODS,");
    for (int i = 0; i < NOR; i++)
        fprintf(fp, "%d,", i + 1);

    fprintf(fp, "\nBNR->ENR,");
    for (int i = 0; i < NOR; i++)
        fprintf(fp, "p%d -> p%d,", BNR[i], ENR[i]);

    fprintf(fp, "\nELASTIC,");
    for (int i = 0; i < NOR; i++)
        fprintf(fp, "%f,", ELASTIC[i]);

    fprintf(fp, "\nSHEAR,");
    for (int i = 0; i < NOR; i++)
        fprintf(fp, "%f,", SHEAR[i]);

    fprintf(fp, "\nAREA,");
    for (int i = 0; i < NOR; i++)
        fprintf(fp, "%f,", AREA[i]);

    fprintf(fp, "\nIMY,");
    for (int i = 0; i < NOR; i++)
        fprintf(fp, "%f,", IMY[i]);

    fprintf(fp, "\nIMZ,");
    for (int i = 0; i < NOR; i++)
        fprintf(fp, "%f,", IMZ[i]);

    fprintf(fp, "\nTHETA,");
    for (int i = 0; i < NOR; i++)
        fprintf(fp, "%f,", THETA[i]);

    //------------LOADS-------------------------------------------------
    fprintf(fp, "\nLOADS,");
    for (int i = 0; i < NOL; i++)
        fprintf(fp, "%d,", i + 1);

    fprintf(fp, "\nPLI,");
    for (int i = 0; i < NOL; i++)
        fprintf(fp, "%d,", PLI[i]);

    fprintf(fp, "\nNRL,");
    for (int i = 0; i < NOL; i++)
        fprintf(fp, "%d,", NRL[i]);

    fprintf(fp, "\nKOL,");
    for (int i = 0; i < NOL; i++)
        fprintf(fp, "%d,", KOL[i]);

    fprintf(fp, "\nVOL,");
    for (int i = 0; i < NOL; i++)
        fprintf(fp, "%f,", VOL[i]);

    fprintf(fp, "\nDLB,");
    for (int i = 0; i < NOL; i++)
        fprintf(fp, "%f,", DLB[i]);

    //-----------SECTIONS-------------------------------------------------
    fprintf(fp, "\nNOS,");
    for (int i = 0; i < NOS; i++)
        fprintf(fp, "%d,", i + 1);

    fprintf(fp, "\nNRS,");
    for (int i = 0; i < NOS; i++)
        fprintf(fp, "%d,", NRS[i]);

    fprintf(fp, "\nDSB,");
    for (int i = 0; i < NOS; i++)
        fprintf(fp, "%f,", DSB[i]);

    //-----------RESULTS OF NODES-----------------------------------------
    fprintf(fp, "\nNORN,");
    for (int i = 0; i < NFRN; i++)
        fprintf(fp, "x%d,y%d,z%d,", i + NFIN, i + NFIN, i + NFIN);

    fprintf(fp, "\nDISPLACEMENT,");
    for (int i = 0; i < NFRN; i++)
        fprintf(fp, "%f,%f,%f,", DON[6 * i], DON[6 * i + 1], DON[6 * i + 2]);

    fprintf(fp, "\nDIVERSION,");
    for (int i = 0; i < NFRN; i++)
        fprintf(fp, "%f,%f,%f,", DON[6 * i + 3], DON[6 * i + 4], DON[6 * i + 5]);

    //-----------RESULTS OF SECTIONS--------------------------------------
    fprintf(fp, "\nNOS,");
    for (int i = 0; i < NOS; i++)
        fprintf(fp, "x%d(AXIAL),y%d(SHEAR),z%d(SHEAR),", i + NOS + 1, i + NOS + 1, i + NOS + 1);

    fprintf(fp, "\nAXIAL&SHEAR FORCE,");
    for (int i = 0; i < NOS; i++)
        fprintf(fp, "%f,%f,%f,", IFS[6 * i], IFS[6 * i + 1], IFS[6 * i + 2]);

    fprintf(fp, "\nTORQUE&BENDING MOMENT,");
    for (int i = 0; i < NOS; i++)
        fprintf(fp, "%f,%f,%f,", IFS[6 * i + 3], IFS[6 * i + 4], IFS[6 * i + 5]);
    fclose(fp);

    return 0;
}

bool SpaceFrame::sfCalculate()
{
    double *lv = 0; //declare total stiffness and load vector

    // int i = 0;
    // printf("TNN: %d\n", TNN);
    // printf("NFIN: %d\n", NFIN);
    // printf("NOR: %d\n", NOR);
    // printf("NOL:% d\n", NOL);
    // printf("NOS: %d\nXCN ", NOS);
    // for (i = 0; i < TNN; i++)
    // {
    //     printf("%f ", XCN[i]);
    // }
    // printf("\nYCN ");
    // for (i = 0; i < TNN; i++)
    // {
    //     printf("%f ", YCN[i]);
    // }
    // printf("\nZCN ");
    // for (i = 0; i < TNN; i++)
    // {
    //     printf("%f ", ZCN[i]);
    // }
    // printf("\nBNR ");
    // for (i = 0; i < NOR; i++)
    // {
    //     printf("%d ", BNR[i]);
    // }
    // printf("\nENR ");
    // for (i = 0; i < NOR; i++)
    // {
    //     printf("%d ", ENR[i]);
    // }
    // printf("\nELASTIC ");
    // for (i = 0; i < NOR; i++)
    // {
    //     printf("%f ", ELASTIC[i]);
    // }
    // printf("\nSHEAR ");
    // for (i = 0; i < NOR; i++)
    // {
    //     printf("%f ", SHEAR[i]);
    // }
    // printf("\nAREA ");
    // for (i = 0; i < NOR; i++)
    // {
    //     printf("%f ", AREA[i]);
    // }
    // printf("\nIMY ");
    // for (i = 0; i < NOR; i++)
    // {
    //     printf("%.11f ", IMY[i]);
    // }
    // printf("\nIMZ ");
    // for (i = 0; i < NOR; i++)
    // {
    //     printf("%.11f ", IMZ[i]);
    // }
    // printf("\nTHETA ");
    // for (i = 0; i < NOR; i++)
    // {
    //     printf("%f ", THETA[i]);
    // }
    // printf("\nNRL ");
    // for (i = 0; i < NOL; i++)
    // {
    //     printf("%d ", NRL[i]);
    // }
    // printf("\nPLI ");
    // for (i = 0; i < NOL; i++)
    // {
    //     printf("%d ", PLI[i]);
    // }
    // printf("\nKOL ");
    // for (i = 0; i < NOL; i++)
    // {
    //     printf("%d ", KOL[i]);
    // }
    // printf("\nVOL ");
    // for (i = 0; i < NOL; i++)
    // {
    //     printf("%f ", VOL[i]);
    // }
    // printf("\nDLB ");
    // for (i = 0; i < NOL; i++)
    // {
    //     printf("%f ", DLB[i]);
    // }
    // printf("\nNRS ");
    // for (i = 0; i < NOS; i++)
    // {
    //     printf("%d ", NRS[i]);
    // }
    // printf("\nDSB ");
    // for (i = 0; i < NOS; i++)
    // {
    //     printf("%f ", DSB[i]);
    // }
    // printf("\nSIGMA_1");
    // for (i = 0; i < NOS; i++)
    // {
    //     printf("%f ", SIGMA_1[i]);
    // }
    // printf("\nSIGMA_2");
    // for (i = 0; i < NOS; i++)
    // {
    //     printf("%f ", SIGMA_2[i]);
    // }
    // printf("\n");

    int dof = 6 * NFRN;

    lv = (double *)malloc(dof * sizeof(double)); //allocate memory for load vector
    memset(lv, 0, dof * sizeof(double));

    if (sfBuildTotalStiff()) //build total stiffness matrix
    {
        sfPrintError(2);
        return 1;
    }
    else
        printf("Building total stiffness matrix succeeded!\n");

    if (sfBuildLoadVector(lv)) //build load stiffness vector
    {
        sfPrintError(3);
        return 1;
    }
    else
        printf("Building load vector succeeded!\n");

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

    if (solve_conjugate_gradient_par(TS, lv, DON, 6 * NFRN)) //solve matrix equation
    {
        sfPrintError(4);
       return 1;
    }
    else
        printf("Solving equation succeeded!\n");

    free(TS);
    free(lv);

    // sfPrintLine2();
    // for (int i = 0; i < 6 * NFRN; i++)
    // {
    //     printf("%15.15f", DON[i]);
    //     printf("\n");
    // }
    // sfPrintLine2();

    for (int i = 0; i < NOS; i++)
        if (sfInternalForce(6 * i, NRS[i], DSB[i])) //calculate the internal force of each rods
        {
            sfPrintError(5);
            return 1;
        }

    // sfPrintLine2();
    // for (int i = 0; i < 6 * NOS; i += 1)
    // {
    //     printf("%15.15f\n", IFS[i]);
    // }
    // sfPrintLine2();

    return 0;
}

int main()
{
    SpaceFrame Frame;
    system("pause");
    Frame.sfInput();
    Frame.sfCalculate();
    Frame.sfOutput();
    system("pause");

    return 0;
}
