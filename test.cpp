#include <Windows.h>
#include <ctype.h>
#include <iostream>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

using namespace std;

bool sfCircularStructure(int n, int m, int l)
{
    FILE *fp = 0;
    fp = fopen("source&result/sf_test.csv", "w");

    fprintf(fp, "Stress Test, degree of freedom is %d,\n", ((m + 1) * (n + 1) * (l + 1) - (m + 1) * (n + 1)) * 6);
    fprintf(fp, "TNN,%d,\n", (m + 1) * (n + 1) * (l + 1));
    fprintf(fp, "NFIN,%d,\n", (m + 1) * (n + 1));
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

    fclose(fp);
    fp = NULL;

    return 0;
}

class SpaceFrame
{
private:
    double EPS;
    double MAXTS;
    double MAXLV;

    int TNN;  // total number of nodes
    int NFIN; // number of fixed nodes
    int NFRN; // number of free nodes
    int NOR;  // number of rods
    int NOL;  // number of loads
    int NOS;  // number of sections

    struct Node // parameters of nodes
    {
        double XCN; // X coordinate of nodes
        double YCN; // Y coordinate of nodes
        double ZCN; // Z coordinate of nodes
    };
    Node *nodes; // parameters of nodes

    struct Rod // parameters of nodes
    {
        int ENR;        // the end node number of rods
        int BNR;        // the beginning node number of rods
        double ELASTIC; // elastic modulus
        double SHEAR;   // shear modulus
        double AREA;    // area
        double IMY;     // inertia moment of Y axis
        double IMZ;     // inertia moment of Z axis
        double THETA;   // theta the deflection angle of main inertia axis
        double LCS[4];  // the length, sine and cosine of rods
        double RFE[6];  // the reaction force of the end node
    };
    Rod *rods; // parameters of nodes

    struct Load // parameters of loads
    {
        int NRL;    // the number of rods with load
        int PLI;    // the plane of the load's in
        int KOL;    // the kind of load
        double VOL; // the value of load
        double DLB; // the distance between load and the beginning node
    };
    Load *loads; // parameters of loads

    struct Section // parameters of sections
    {
        int NRS;       // the number of rod with section
        double DSB;    // the distance between section and the beginning node
        double IFS[6]; // the internal force in the section
    };
    Section *sections; // parameters of sections

    double *TotalStiffness; // total stiffness
    double *LoadVector;     // load vector
    double *Displacement;   // displacement of nodes

    int *IV;     // the location of diagonal element
    int NSI;     // upper limit
    int MAXIBDW; // half bandwidth

    bool ProgressBar = 0;

    // calculate the length sine and cosine of rods
    bool sfLCosSin()
    {
        for (int k = 0; k < NOR; k++)
        {
            int i = rods[k].BNR - 1, j = rods[k].ENR - 1; // index of beginning and end nodes of rods
            rods[k].LCS[1] = nodes[j].XCN - nodes[i].XCN;
            rods[k].LCS[2] = nodes[j].YCN - nodes[i].YCN;
            rods[k].LCS[3] = nodes[j].ZCN - nodes[i].ZCN;
            rods[k].LCS[0] = sqrt(rods[k].LCS[1] * rods[k].LCS[1] + rods[k].LCS[2] * rods[k].LCS[2] + rods[k].LCS[3] * rods[k].LCS[3]);
            if (rods[k].LCS[0] < EPS) // if the length of rod is too small, then return error
            {
                sfPrintError(8);
                return 1;
            }
            rods[k].LCS[1] = rods[k].LCS[1] / rods[k].LCS[0];
            rods[k].LCS[2] = rods[k].LCS[2] / rods[k].LCS[0];
            rods[k].LCS[3] = rods[k].LCS[3] / rods[k].LCS[0];
        }

        return 0;
    }
    // allocate total stiffness matrix, load vector and displacement vector
    bool sfAllocate()
    {
        int it = 0, mm = 0, dof = 6 * NFRN;
        int *peribdw = new int[TNN](); // bandwidth per line in total stiffness matrix
        IV = new int[dof]();
        for (int i = 0; i < NOR; i++) // for each rod
        {
            if (rods[i].BNR > NFIN)
            {
                mm = rods[i].ENR - rods[i].BNR; // bandwidth is end number minus begin number
                if (mm > peribdw[rods[i].ENR - 1])
                    peribdw[rods[i].ENR - 1] = mm; // find the maximum bandwith per line
            }
        }
        for (int i = NFIN; i < TNN; i++) // for each line in total stiffness matrix
        {
            if (peribdw[i] > MAXIBDW) // find maxim
                MAXIBDW = peribdw[i];
            for (int j = 1; j <= 6; j++)
            {
                it = it + 1;
                if (it == 1)
                    IV[it - 1] = 6 * peribdw[i] + j;
                else
                    IV[it - 1] = IV[it - 2] + 6 * peribdw[i] + j;
            }
        }
        MAXIBDW = 6 * MAXIBDW + 5;
        NSI = IV[dof - 1];
        delete[] peribdw;

        TotalStiffness = new double[NSI](); // allocate memory for total stiffness matrix
        LoadVector = new double[dof]();     // allocate memory for load vector
        Displacement = new double[dof]();   // allocate memory for displacement vector

        return 0;
    }
    // build total stiffness matrix
    bool sfBuildTotalStiff() // ts is total stiffness matrix
    {
        double us[36] = {0}; // unit stiffness matrix
        int p[2] = {0};      // p is a temperary vector for i0j0, dof is the degree of freedom of nods

        for (int k = 0; k < NOR; k++)
        {
            p[0] = 6 * (rods[k].BNR - NFIN - 1); // match the displacement with nods
            p[1] = 6 * (rods[k].ENR - NFIN - 1);

            for (int i = 0; i < 2; i++)
            {
                if (p[i] >= 0) // determine free node
                {
                    if (sfBuildUnitStiff(k, i + 1, us)) // build unit stiffness matrix
                    {
                        sfPrintError(7);
                        return 1;
                    }
                    for (int m = 0; m < 6; m++)
                        for (int n = 0; n <= m; n++)
                            TotalStiffness[IV[(p[i] + m)] + (p[i] + n) - (p[i] + m) - 1] += us[m * 6 + n]; // superpose
                }
            }
            if (p[0] >= 0 && p[1] >= 0)
            {
                if (sfBuildUnitStiff(k, 3 + 1, us)) // build unit stiffness matrix
                {
                    sfPrintError(7);
                    return 1;
                }
                for (int m = 0; m < 6; m++)
                    for (int n = 0; n < 6; n++)
                        TotalStiffness[IV[(p[1] + m)] + (p[0] + n) - (p[1] + m) - 1] += us[m * 6 + n]; // superpose
            }
        }

        for (int i = 0; i < NSI; i++)
            if (fabs(TotalStiffness[i]) > MAXTS)
                MAXTS = TotalStiffness[i];

        return 0;
    }
    // build unit stiffness matrix
    bool sfBuildUnitStiff(int k, int flag, double *us) // k is the number of rods, flag is the index of matrix parts, us is the unit stiffness matrix
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

        double rd[36] = {0}, t[36] = {0}, c[36] = {0}, tmp = 0; // rd is local stiffness matrix, t is transpose matrix, c is a temperary matrix
        memset(us, 0, 36 * sizeof(double));

        if (sfBuildLocalStiff(k, flag, rd)) // build local stiffness matrix
        {
            sfPrintError(9);
            return 1;
        }
        if (sfBuildTrans(k, t)) // build transpose matrix
        {
            sfPrintError(10);
            return 1;
        }

        for (int i = 0; i < 6; i++) // transpose matrix times local stiffness matrix, store the result in c
            for (int m = 0; m < 6; m++)
            {
                tmp = t[i * 6 + m];
                for (int j = 0; j < 6; j++)
                    c[i * 6 + j] += tmp * rd[m * 6 + j];
            }

        for (int i = 0; i < 6; i++) // c times the transposition of transpose matrix, store the result in unit stiff
            for (int j = 0; j < 6; j++)
                for (int m = 0; m < 6; m++)
                    us[i * 6 + j] += c[i * 6 + m] * t[j * 6 + m];

        return 0;
    }
    // build local stiffness matrix
    bool sfBuildLocalStiff(int k, int flag, double *rd) // k is the number of rods, flag is the number of matrix
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

        double a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0, l = rods[k].LCS[0];

        a = rods[k].ELASTIC * rods[k].AREA / l;              // EA/1
        b = rods[k].SHEAR * (rods[k].IMY + rods[k].IMZ) / l; // GJ(p)/1
        c = 4 * rods[k].ELASTIC * rods[k].IMY / l;           // 4EJ(y)/1
        d = c / 2 * 3 / l;                                   // 6EJ(z)/l/l
        e = 2 * d / l;                                       // 12EJ(y)/l/l/l
        f = 4 * rods[k].ELASTIC * rods[k].IMZ / l;           // 4EJ(z)/l
        g = f / 2 * 3 / l;                                   // 6EJ(Z)/l/l
        h = 2 * g / l;                                       // 12EJ(z)/l/l/l

        switch (flag)
        {
        case 1: // k11
            rd[0 * 6 + 0] = a;
            rd[1 * 6 + 1] = h;
            rd[1 * 6 + 5] = rd[5 * 6 + 1] = g;
            rd[2 * 6 + 2] = e;
            rd[2 * 6 + 4] = rd[4 * 6 + 2] = -d;
            rd[3 * 6 + 3] = b;
            rd[4 * 6 + 4] = c;
            rd[5 * 6 + 5] = f;
            break;
        case 2: // k22
            rd[0 * 6 + 0] = a;
            rd[1 * 6 + 1] = h;
            rd[1 * 6 + 5] = rd[5 * 6 + 1] = -g;
            rd[2 * 6 + 2] = e;
            rd[2 * 6 + 4] = rd[4 * 6 + 2] = d;
            rd[3 * 6 + 3] = b;
            rd[4 * 6 + 4] = c;
            rd[5 * 6 + 5] = f;
            break;
        case 3: // k12
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
        case 4: // k21
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
    // build transpose matrix
    bool sfBuildTrans(int k, double *t) // k is the number of rods, t is transpose matrix
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

        double coa = 0, cob = 0, coc = 0, sic = 0, sit = 0, cot = 0, m = 0, n = 0; // co means cosine, si means sine, m and n is temperary variable

        memset(t, 0, 36 * sizeof(double));

        coa = rods[k].LCS[1];     // cosine alpha
        cob = rods[k].LCS[2];     // cosine beta
        coc = rods[k].LCS[3];     // cosine gama
        sit = sin(rods[k].THETA); // sine theta
        cot = cos(rods[k].THETA); // cosine theta

        if (fabs(coc - 1) < EPS) // vertical(z axis positive direction) rods' transpose matrix
        {
            t[2 * 6 + 0] = t[5 * 6 + 3] = 1;
            t[0 * 6 + 1] = t[3 * 6 + 4] = t[1 * 6 + 2] = t[4 * 6 + 5] = sit;
            t[1 * 6 + 1] = t[4 * 6 + 4] = cot;
            t[0 * 6 + 2] = t[3 * 6 + 5] = -cot;
        }
        else if (fabs(coc + 1) < EPS) // vertical(z axis negative direction) rods' transpose matrix
        {
            t[2 * 6 + 0] = t[5 * 6 + 3] = -1;
            t[0 * 6 + 1] = t[3 * 6 + 4] = sit;
            t[1 * 6 + 2] = t[4 * 6 + 5] = -sit;
            t[1 * 6 + 1] = t[4 * 6 + 4] = t[0 * 6 + 2] = t[3 * 6 + 5] = cot;
        }
        else
        {
            sic = sqrt(1 - coc * coc); // sine gama
            m = coa * coc;             // cosine alpha times cosine gama
            n = cob * coc;             // cosine beta times cosine gama

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
    // build load vector
    bool sfBuildLoadVector(double *lv) // lv is the load vector
    {
        if (lv == 0)
        {
            sfPrintError(19);
            return 0;
        }

        int rod = 0, p[2] = {0};          // rod is the number of rods, dof is the degree of freedom
        double rf[12] = {0}, t[36] = {0}; // rf is the reaction force matrix, t is the transpose matrix, p is a temperary vector for i0j0

        for (int i = 0; i < NOL; i++)
        {
            rod = loads[i].NRL - 1;             // the number of rods with load
            memset(rf, 0, 12 * sizeof(double)); // zero clearing

            if (sfReactionForce(i, &rf[0 * 6], &rf[1 * 6])) // calculate reaction force
            {
                sfPrintError(11);
                return 1;
            }
            for (int j = 0; j < 6; j++) // add reaction force to RFE
                rods[rod].RFE[j] += rf[1 * 6 + j];
            if (sfBuildTrans(rod, t)) // build transpose matrix
            {
                sfPrintError(10);
                return 1;
            }

            p[0] = 6 * (rods[rod].BNR - NFIN - 1); // match the displacement with nods
            p[1] = 6 * (rods[rod].ENR - NFIN - 1);

            for (int j = 0; j < 2; j++) // add reaction force to load vector
            {
                if (p[j] >= 0) // determine free node
                    for (int m = 0; m < 6; m++)
                        for (int n = 0; n < 6; n++)
                            lv[p[j] + m] -= t[m * 6 + n] * rf[j * 6 + n];
            }
        }

        for (int i = 0; i < 6 * NFRN; i++)
            if (fabs(lv[i]) > MAXLV)
                MAXLV = lv[i];

        return 0;
    }
    // calculate reaction force
    bool sfReactionForce(int i, double *rfb, double *rfe) // i is the number of load, rfb and rfe is the reaction force at begining and end of rods
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

        double ra = 0, rb = 0, a = 0, b = 0, q = loads[i].VOL, xq = loads[i].DLB; // ra, rb, a and b are middle variable
        int rod = loads[i].NRL - 1, pm = loads[i].PLI, t = 0;                     // rod is the number of rods

        if (pm == 0)      // load is in XY plane
            t = -1;       // The bending moment in the support-reaction equation is positive clockwise, convert it to positive to the coordinate axis
        else if (pm == 1) // load is in XZ plane
            t = 1;        // The bending moment in the support-reaction equation is positive clockwise, convert it to positive to the coordinate axis

        ra = loads[i].DLB / rods[rod].LCS[0]; // x(q) / L
        rb = 1 - ra;                          // 1 - x(q) / L
        switch (loads[i].KOL)
        {
        case 1: // vertical concentrating load
            a = rb * rb;
            rfb[pm + 1] = -q * rb * (1 + ra - 2 * ra * ra);
            rfe[pm + 1] = -q - rfb[pm + 1];
            rfb[5 - pm] = t * q * rb * ra * (rods[rod].LCS[0] - xq);
            rfe[5 - pm] = -t * q * ra * rb * xq;
            break;
        case 2: // vertical uniform load
            a = q * xq;
            b = a * xq / 12;
            rfb[pm + 1] = -a * (1 + 0.5 * ra * ra * ra - ra * ra);
            rfe[pm + 1] = -a - rfb[pm + 1];
            rfb[5 - pm] = t * b * (6 - 8 * ra + 3 * ra * ra);
            rfe[5 - pm] = -t * b * (4 * ra - 3 * ra * ra);
            break;
        case 3: // axial concentrating force when PLI == 0, torque when PLI ==1
            rfb[3 * pm] = -q * rb;
            rfe[3 * pm] = -q * ra;
            break;
        case 4: // axial uniform load
            a = q * xq;
            rfe[3 * pm] = -a * ra / 2;
            rfb[3 * pm] = -a - rfe[3 * pm];
            break;
        case 5: // vertical triangle distributed load
            a = q * xq / 2;
            b = -0.4 * ra * ra;
            rfb[pm + 1] = -2 * a * (0.5 - 0.75 * ra * ra + 0.4 * ra * ra * ra);
            rfe[pm + 1] = -a - rfb[pm + 1];
            rfb[5 - pm] = t * a * (2 / 3 + b - ra);
            rfe[5 - pm] = -t * a * (0.5 * ra + b);
            break;
        case 6: // concentrating bending moment
            rfb[2 - pm] = t * 6 * q * rb * ra / rods[rod].LCS[0];
            rfe[2 - pm] = -rfb[2 - pm];
            rfb[pm + 4] = t * q * rb * (-1 + 3 * ra);
            rfe[pm + 4] = t * q * ra * (2 - 3 * ra);
            break;
        case 7: // unifrom temperature rise
            rfb[0] = q * xq * rods[rod].ELASTIC * rods[rod].AREA;
            rfe[0] = -rfb[0];
            break;
        case 8: // different temperature rise
            if (pm == 0)
                a = rods[rod].IMZ;
            else if (pm == 1)
                a = rods[rod].IMY;
            rfb[5 - pm] = t * q * 2 * rods[rod].ELASTIC * a * xq;
            rfe[5 - pm] = -rfb[5 - pm];
            break;
        default:
            break;
        }

        return 0;
    }
    // solve equation of matrix by conjugate gradient parallel
    bool sfConjugateGradientPar(double *A, double *b, double *x, int N)
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

        double *r = NULL, *p = NULL, *z = NULL;
        double gamma = 0, gamma_new = 0, gamma_new_sqrt = 0, alpha = 0, beta = 0;
        int percent = 0, percent_new = 0;

        if (ProgressBar)
        {
            printf("\rSolving equation      [ 0%% ][                                                 ]");
        }

        r = (double *)malloc(N * sizeof(double));
        memset(r, 0, sizeof(double));
        p = (double *)malloc(N * sizeof(double));
        memset(p, 0, sizeof(double));
        z = (double *)malloc(N * sizeof(double));
        memset(z, 0, sizeof(double));

        for (int i = 0; i < NSI; i++)
            A[i] = A[i] / MAXTS;
        for (int i = 0; i < N; i++)
            b[i] = b[i] / MAXLV;

        //  x = [0 ... 0]
        //  r = b - A * x
        //  p = r
        //  gamma = r' * r
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
//  z = A * p
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
                            z[i] += A[IV[j] - j + i - 1] * p[j];
                        else
                            z[i] += 0;
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

            //  alpha = gamma / (p' * z)
            alpha = 0.0;
#pragma omp parallel for reduction(+ \
                                   : alpha)
            for (int i = 0; i < N; ++i)
                alpha += p[i] * z[i];
            alpha = gamma / alpha;

            //  x = x + alpha * p
            //  r = r - alpha * z
            //  gamma_new = r' * r
            gamma_new = 0.0;
#pragma omp parallel for reduction(+ \
                                   : gamma_new)
            for (int i = 0; i < N; ++i)
            {
                x[i] += alpha * p[i];
                r[i] -= alpha * z[i];
                gamma_new += r[i] * r[i];
            }

            gamma_new_sqrt = sqrt(gamma_new);
            if (gamma_new_sqrt < EPS)
                break;

            if (ProgressBar)
            {
                percent_new = (int)((1 - log10(gamma_new_sqrt * 1e15) / 16) * 100);
                if (percent_new > percent)
                {
                    percent = percent_new;
                    printf("\rSolving equation ");
                    for (int i = 0; i <= 4; i++)
                        if (i <= n % 4)
                            printf(".");
                        else
                            printf(" ");
                    printf("[ %d%% ]", percent_new);
                    printf("[");
                    for (int i = 0; i < 49; i++)
                        if (i < percent / 2)
                            printf("=");
                        else
                            printf(" ");
                    printf("]");
                }
                else
                {
                    printf("\rSolving equation ");
                    for (int i = 0; i <= 4; i++)
                        if (i <= n % 4)
                            printf(".");
                        else
                            printf(" ");
                }
            }

            beta = gamma_new / gamma;

//  p = r + (gamma_new / gamma) * p;
#pragma omp parallel for
            for (int i = 0; i < N; ++i)
                p[i] = r[i] + beta * p[i];

            //  gamma = gamma_new
            gamma = gamma_new;
        }

        if (ProgressBar)
        {
            printf("\rSolving equation done [ 100%% ][=================================================]\n");
        }

        for (int i = 0; i < NSI; i++)
            A[i] = A[i] * MAXTS;
        for (int i = 0; i < N; i++)
            b[i] = b[i] * MAXLV;
        for (int i = 0; i < N; i++)
            x[i] = x[i] * MAXLV / MAXTS;

        free(r);
        free(p);
        free(z);

        return 0;
    }
    // calculate internal force of rods
    bool sfInternalForce(int mm, int k, double xp) // m is the number of sections, k is the actual number of rods, xp is the distance between the section and the begining of rods
    {
        if (mm < 0)
        {
            sfPrintError(21);
            return 0;
        }
        if (k < 0)
        {
            sfPrintError(21);
            return 0;
        }

        double tf[6] = {0}; // tf is temperary variable

        sections[mm].IFS[0] = +rods[k - 1].RFE[0]; // calculate internal force cause by reaction force at the end of rods
        sections[mm].IFS[1] = -rods[k - 1].RFE[1];
        sections[mm].IFS[2] = -rods[k - 1].RFE[2];
        sections[mm].IFS[3] = +rods[k - 1].RFE[3];
        sections[mm].IFS[4] = -rods[k - 1].RFE[4] + rods[k - 1].RFE[2] * (rods[k - 1].LCS[0] - xp);
        sections[mm].IFS[5] = +rods[k - 1].RFE[5] + rods[k - 1].RFE[1] * (rods[k - 1].LCS[0] - xp);

        for (int i = 0; i < NOL; i++) // for every rods
            if (loads[i].NRL == k)    // if load is on rod k
            {
                memset(tf, 0, 6 * sizeof(double)); // zero clear tf
                if (sfCtlInternalForce(i, xp, tf)) //  calculate internal force of cantilever beam
                {
                    sfPrintError(13);
                    return 1;
                }
                for (int j = 0; j < 6; j++) // add internal force of cantilever into IFR
                    sections[mm].IFS[j] += tf[j];
            }

        if (sfDisplacementForce(k, tf)) // calculate end force
        {
            sfPrintError(14);
            return 1;
        }

        sections[mm].IFS[0] -= tf[0]; // calculate section force cause by end force
        sections[mm].IFS[1] += tf[1];
        sections[mm].IFS[2] += tf[2];
        sections[mm].IFS[3] -= tf[3];
        sections[mm].IFS[4] += tf[4] + tf[2] * xp;
        sections[mm].IFS[5] -= tf[5] - tf[1] * xp;

        return 0;
    }
    // calculate internal force of cantilever beam
    bool sfCtlInternalForce(int i, double xp, double *tf) // i is the number of load, xp is the distance between the section and the begining of rod, tf is internal force
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

        double xq = loads[i].DLB, t = xq - xp, r = xp / xq, q = loads[i].VOL; // t and r are temperary variables
        int e = loads[i].PLI;
        switch (loads[i].KOL) // calculate section force according to kind of loads
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
                tf[3 * e] = q;
            break;
        case 4:
            if (xp < xq)
                tf[3 * e] = q * t;
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
                tf[e + 4] = (2 * e - 1) * q;
            break;
        case 7: // temperature change don't generate internal force on cantilever beam
            break;
        case 8:
            break;

        default:
            break;
        }

        return 0;
    }
    // calculate internal force of displacement
    bool sfDisplacementForce(int k, double *tref) // k is the actual number of rods, tref is the end force of rods
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

        int p[2] = {0};                                  // p is a temperary vector for i0j0
        double rd[36] = {0}, rdb[36] = {0}, t[36] = {0}; // rd

        memset(tref, 0, 6 * sizeof(double));

        if (sfBuildTrans(k - 1, t)) // calculate transpose matrix
        {
            sfPrintError(10);
            return 1;
        }

        p[0] = 6 * (rods[k - 1].BNR - NFIN - 1); // match the displacement with nods
        p[1] = 6 * (rods[k - 1].ENR - NFIN - 1);

        for (int i = 0; i < 2; i++)
        {
            if (p[i] >= 0) // determine free node
            {
                if (sfBuildLocalStiff(k - 1, 2 * i + 1, rd)) // build unit stiffness matrix
                {
                    sfPrintError(9);
                    return 1;
                }

                memset(rdb, 0, 36 * sizeof(double)); // zero clean rdb

                for (int j = 0; j < 6; j++) // rd times transposition of transpose matrix
                    for (int m = 0; m < 6; m++)
                        for (int n = 0; n < 6; n++)
                            rdb[j * 6 + m] += rd[j * 6 + n] * t[m * 6 + n];
                for (int j = 0; j < 6; j++) // rdb times DON
                    for (int m = 0; m < 6; m++)
                        tref[j] += rdb[j * 6 + m] * Displacement[p[i] + m];
            }
            else // fixed node
                for (int j = 0; j < 3; j++)
                    tref[j] += 0;
        }

        return 0;
    }
    // print"----------------------------------------"
    bool sfPrintLine()
    {
        printf("--------------------------------------------------------------------------\n");
        return 0;
    }
    // print"****************************************"
    bool sfPrintLine2()
    {
        printf("**************************************************************************\n");
        return 0;
    }
    // print error
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
    SpaceFrame(SpaceFrame &);
    ~SpaceFrame();

    // read data from .csv
    bool sfInput();
    // calculate
    bool sfCalculate(bool);
    // output data
    bool sfOutput();
};

SpaceFrame::SpaceFrame()
{
    EPS = 1e-15;
    MAXTS = 0;
    MAXLV = 0;

    TNN = 0;  // total number of nodes
    NFIN = 0; // number of fixed nodes
    NFRN = 0; // number of free nodes
    NOR = 0;  // number of rods
    NOL = 0;  // number of loads
    NOS = 0;  // number of sections

    nodes = NULL;    // parameters of nodes
    rods = NULL;     // parameters of rods
    loads = NULL;    // parameters of loads
    sections = NULL; // parameters of sections

    TotalStiffness = NULL; // total stiffness
    LoadVector = NULL;     // load vector
    Displacement = NULL;   // the displacement of nodes

    IV = NULL;   // the location of diagonal element
    NSI = 0;     // upper limit
    MAXIBDW = 0; // half bandwidth
}

SpaceFrame::SpaceFrame(SpaceFrame &Frame)
{
    EPS = Frame.EPS;
    MAXTS = Frame.MAXTS;
    MAXLV = Frame.MAXLV;

    TNN = Frame.TNN;
    NFIN = Frame.NFIN;
    NFRN = Frame.NFRN;
    NOR = Frame.NOR;
    NOL = Frame.NOL;
    NOS = Frame.NOS;

    nodes = new Node[TNN]();
    if (Frame.nodes != NULL)
        memcpy(nodes, Frame.nodes, TNN * sizeof(Node));

    rods = new Rod[NOR]();
    if (Frame.rods != NULL)
        memcpy(rods, Frame.rods, NOR * sizeof(Rod));

    loads = new Load[NOL]();
    if (Frame.loads != NULL)
        memcpy(loads, Frame.loads, NOL * sizeof(Load));

    sections = new Section[NOS]();
    if (Frame.sections != NULL)
        memcpy(sections, Frame.sections, NOS * sizeof(Section));

    int dof = 6 * NFRN;
    IV = new int[dof]();
    if (Frame.IV != NULL)
        memcpy(IV, Frame.IV, dof * sizeof(int));

    NSI = Frame.NSI;
    MAXIBDW = Frame.MAXIBDW;

    TotalStiffness = new double[NSI]();
    if (Frame.TotalStiffness != NULL)
        memcpy(TotalStiffness, Frame.TotalStiffness, NSI * sizeof(double));

    LoadVector = new double[dof]();
    if (Frame.LoadVector != NULL)
        memcpy(LoadVector, Frame.LoadVector, dof * sizeof(double));

    Displacement = new double[dof]();
    if (Frame.Displacement != NULL)
        memcpy(Displacement, Frame.Displacement, dof * sizeof(double));
}

SpaceFrame::~SpaceFrame()
{
    delete[] nodes;
    nodes = NULL;
    delete[] rods;
    rods = NULL;
    delete[] loads;
    loads = NULL;
    delete[] sections;
    sections = NULL;
    delete[] TotalStiffness;
    TotalStiffness = NULL;
    delete[] LoadVector;
    LoadVector = NULL;
    delete[] Displacement;
    Displacement = NULL;
    delete[] IV;
    IV = NULL;
}

bool SpaceFrame::sfInput()
{
    FILE *fp = NULL;                   // Define the file point
    char *line = 0, *data = 0;         // Define the line string and separated string
    char temporSpace[1000000];         // Apply for temporary storage space
    int rowIndex = 0, columnIndex = 0; // Reset the number of rows to zero, reset the number of columns to zero
    const char DIVIDE[] = ",";         // Set the separater as a ','

    if ((fp = fopen("source&result/sf_test.csv", "r")) == NULL) // Start the process when the file opens successfully
    {
        printf("There is no such file!");
        return 0;
    }

    fseek(fp, 0L, SEEK_SET);                                             // Locate file point to the first line
    while ((line = fgets(temporSpace, sizeof(temporSpace), fp)) != NULL) // The loop continues when the end of the file is not read
    {
        data = strtok(line, DIVIDE); // Split strings with a ',' as a separator
        while (data != NULL)         // Read the data of each row
        {
            if (strcmp(data, "END") == 0) // When the keyword 'END' is read, the reading process will be shut down
            {
                fclose(fp); // Close the file
                fp = NULL;  // Reset the file point
                cout << "Inputing data succeed!\n";

                return 0;
            }

            if (columnIndex++ == 0) // Skip the saving of the first column
            {
                data = strtok(NULL, DIVIDE); // Reset data
                continue;
            }

            switch (rowIndex) // Store variables of each column in different ways
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

                    if (nodes != NULL)
                        this->~SpaceFrame();

                    nodes = new Node[TNN]();
                    rods = new Rod[NOR]();
                    loads = new Load[NOL]();
                    sections = new Section[NOS]();
                }
                break;
            case 6:
                if (columnIndex - 2 < TNN)
                    nodes[columnIndex - 2].XCN = atof(data);
                break;
            case 7:
                if (columnIndex - 2 < TNN)
                    nodes[columnIndex - 2].YCN = atof(data);
                break;
            case 8:
                if (columnIndex - 2 < TNN)
                    nodes[columnIndex - 2].ZCN = atof(data);
                break;
            case 9:
                if (columnIndex - 2 < NOR)
                    rods[columnIndex - 2].BNR = atoi(data);
                break;
            case 10:
                if (columnIndex - 2 < NOR)
                    rods[columnIndex - 2].ENR = atoi(data);
                break;
            case 11:
                if (columnIndex - 2 < NOR)
                    rods[columnIndex - 2].ELASTIC = atof(data);
                break;
            case 12:
                if (columnIndex - 2 < NOR)
                    rods[columnIndex - 2].SHEAR = atof(data);
                break;
            case 13:
                if (columnIndex - 2 < NOR)
                    rods[columnIndex - 2].AREA = atof(data);
                break;
            case 14:
                if (columnIndex - 2 < NOR)
                    rods[columnIndex - 2].IMY = atof(data);
                break;
            case 15:
                if (columnIndex - 2 < NOR)
                    rods[columnIndex - 2].IMZ = atof(data);
                break;
            case 16:
                if (columnIndex - 2 < NOR)
                    rods[columnIndex - 2].THETA = atof(data);
                break;
            case 17:
                if (columnIndex - 2 < NOL)
                    loads[columnIndex - 2].NRL = atoi(data);
                break;
            case 18:
                if (columnIndex - 2 < NOL)
                    loads[columnIndex - 2].PLI = atoi(data);
                break;
            case 19:
                if (columnIndex - 2 < NOL)
                    loads[columnIndex - 2].KOL = atoi(data);
                break;
            case 20:
                if (columnIndex - 2 < NOL)
                    loads[columnIndex - 2].VOL = atof(data);
                break;
            case 21:
                if (columnIndex - 2 < NOL)
                    loads[columnIndex - 2].DLB = atof(data);
                break;
            case 22:
                if (columnIndex - 2 < NOS)
                    sections[columnIndex - 2].NRS = atoi(data);
                break;
            case 23:
                if (columnIndex - 2 < NOS)
                    sections[columnIndex - 2].DSB = atof(data);
                break;
            } // input finished

            data = strtok(NULL, DIVIDE); // Reset data
        }
        rowIndex++;      // RowIndex steps forward once
        columnIndex = 0; // Reset columnIndex
    }
    fclose(fp); // Close the file
    fp = NULL;  // Reset the file point

    return 0;
}

bool SpaceFrame::sfOutput()
{
    if (false) // console
    {
        printf("\t\t\tCalculation Of Space Rigid Frame\n");
        sfPrintLine();

        printf("\t\tTNN = %d\t\t\tNFIN = %d\n\t\tNFRN = %d\t\tNOR = %d\n", TNN, NFIN, NFRN, NOR);
        printf("\t\tNOL = %d\t\t\tNOS = %d\n", NOL, NOS);
        sfPrintLine();

        printf("NUMBER OF NODES     Coordinate-X    Coordinate-Y    Coordinate-Z\n");
        for (int i = 0; i < TNN; i++)
            printf("%15d%15.7f%15.7f%15.7f\n", i + 1, nodes[i].XCN, nodes[i].YCN, nodes[i].ZCN);
        sfPrintLine();

        printf("NUMBER OF NODES     LEFT NODES    RIGHT NODES  Elastic modulus  Shear modulus    Area   Inertia moment Y axis  Inertia moment Z axis\n");
        for (int i = 0; i < NOR; i++)
            printf("%15d%15d%15d%15.0f%15.0f%11.4f%16.5f%23.5f\n", i + 1, rods[i].BNR, rods[i].ENR, rods[i].ELASTIC, rods[i].SHEAR, rods[i].AREA, rods[i].IMY, rods[i].IMZ);
        sfPrintLine();
        printf("NUMBER OF SECTIONS         PLI            DLB\n");
        for (int i = 0; i < NOS; i++)
            printf("%15d%15d%15.7f\n", i + 1, sections[i].NRS, sections[i].DSB);
        sfPrintLine();

        printf("Calculating......\nThe results are as follows: \n");
        sfPrintLine();

        printf("NUMBER OF NODES   Displacement-X Displacement-Y Displacement-Z   Diversion-X    Diversion-Y    Diversion-Z\n");
        for (int i = NFIN; i < TNN; i++)
            printf("%15d%15.7f%15.7f%15.7f%15.7f%15.7f%15.7f\n", i + 1, Displacement[6 * (i - NFIN)], Displacement[6 * (i - NFIN) + 1], Displacement[6 * (i - NFIN) + 2], Displacement[6 * (i - NFIN) + 3], Displacement[6 * (i - NFIN) + 4], Displacement[6 * (i - NFIN) + 5]);
        sfPrintLine();

        printf("NUMBER OF SECTIONS Axial force-X  Shear force-Y  Shear force-Z    Torque-X   Bending moment-Y  Bending moment-Z\n");
        for (int i = 0; i < NOS; i++)
            printf("%15d%15.7f%15.7f%15.7f%15.7f%15.7f%15.7f\n", i + 1, sections[i].IFS[0], sections[i].IFS[1], sections[i].IFS[2], sections[i].IFS[3], sections[i].IFS[4], sections[i].IFS[5]);
    }

    if (true) // file
    {
        FILE *fp = NULL;
        fp = fopen("source&result/sfResultClass.csv", "w");
        fprintf(fp, "TITLE\n");
        fprintf(fp, "TNN,%d\nNFIN,%d\nNFRN,%d\nNOR,%d\nNOL,%d\nNOS,%d", TNN, NFIN, NFRN, NOR, NOL, NOS);

        // ------------NODES-------------------------------------------------
        fprintf(fp, "\nNODES,");
        for (int i = 0; i < TNN; i++)
            fprintf(fp, "%d,", i + 1);

        fprintf(fp, "\nCON,");
        for (int i = 0; i < TNN; i++)
            fprintf(fp, "(%f %f %f),", nodes[i].XCN, nodes[i].YCN, nodes[i].ZCN);

        // ------------RODS-------------------------------------------------
        fprintf(fp, "\nRODS,");
        for (int i = 0; i < NOR; i++)
            fprintf(fp, "%d,", i + 1);

        fprintf(fp, "\nBNR->ENR,");
        for (int i = 0; i < NOR; i++)
            fprintf(fp, "p%d -> p%d,", rods[i].BNR, rods[i].ENR);

        fprintf(fp, "\nELASTIC,");
        for (int i = 0; i < NOR; i++)
            fprintf(fp, "%f,", rods[i].ELASTIC);

        fprintf(fp, "\nSHEAR,");
        for (int i = 0; i < NOR; i++)
            fprintf(fp, "%f,", rods[i].SHEAR);

        fprintf(fp, "\nAREA,");
        for (int i = 0; i < NOR; i++)
            fprintf(fp, "%f,", rods[i].AREA);

        fprintf(fp, "\nIMY,");
        for (int i = 0; i < NOR; i++)
            fprintf(fp, "%f,", rods[i].IMY);

        fprintf(fp, "\nIMZ,");
        for (int i = 0; i < NOR; i++)
            fprintf(fp, "%f,", rods[i].IMZ);

        fprintf(fp, "\nTHETA,");
        for (int i = 0; i < NOR; i++)
            fprintf(fp, "%f,", rods[i].THETA);

        // ------------LOADS-------------------------------------------------
        fprintf(fp, "\nLOADS,");
        for (int i = 0; i < NOL; i++)
            fprintf(fp, "%d,", i + 1);

        fprintf(fp, "\nPLI,");
        for (int i = 0; i < NOL; i++)
            fprintf(fp, "%d,", loads[i].PLI);

        fprintf(fp, "\nNRL,");
        for (int i = 0; i < NOL; i++)
            fprintf(fp, "%d,", loads[i].NRL);

        fprintf(fp, "\nKOL,");
        for (int i = 0; i < NOL; i++)
            fprintf(fp, "%d,", loads[i].KOL);

        fprintf(fp, "\nVOL,");
        for (int i = 0; i < NOL; i++)
            fprintf(fp, "%f,", loads[i].VOL);

        fprintf(fp, "\nDLB,");
        for (int i = 0; i < NOL; i++)
            fprintf(fp, "%f,", loads[i].DLB);

        // -----------SECTIONS-------------------------------------------------
        fprintf(fp, "\nNOS,");
        for (int i = 0; i < NOS; i++)
            fprintf(fp, "%d,", i + 1);

        fprintf(fp, "\nNRS,");
        for (int i = 0; i < NOS; i++)
            fprintf(fp, "%d,", sections[i].NRS);

        fprintf(fp, "\nDSB,");
        for (int i = 0; i < NOS; i++)
            fprintf(fp, "%f,", sections[i].DSB);

        // -----------RESULTS OF NODES-----------------------------------------
        fprintf(fp, "\nNFRN,");
        for (int i = 0; i < NFRN; i++)
            fprintf(fp, "x%d,y%d,z%d,", i + NFIN, i + NFIN, i + NFIN);

        fprintf(fp, "\nDISPLACEMENT,");
        for (int i = 0; i < NFRN; i++)
            fprintf(fp, "%f,%f,%f,", Displacement[6 * i], Displacement[6 * i + 1], Displacement[6 * i + 2]);

        fprintf(fp, "\nDIVERSION,");
        for (int i = 0; i < NFRN; i++)
            fprintf(fp, "%f,%f,%f,", Displacement[6 * i + 3], Displacement[6 * i + 4], Displacement[6 * i + 5]);

        // -----------RESULTS OF SECTIONS--------------------------------------
        fprintf(fp, "\nNOS,");
        for (int i = 0; i < NOS; i++)
            fprintf(fp, "x%d(AXIAL),y%d(SHEAR),z%d(SHEAR),", i + NOS + 1, i + NOS + 1, i + NOS + 1);

        fprintf(fp, "\nAXIAL&SHEAR FORCE,");
        for (int i = 0; i < NOS; i++)
            fprintf(fp, "%f,%f,%f,", sections[i].IFS[0], sections[i].IFS[1], sections[i].IFS[2]);

        fprintf(fp, "\nTORQUE&BENDING MOMENT,");
        for (int i = 0; i < NOS; i++)
            fprintf(fp, "%f,%f,%f,", sections[i].IFS[3], sections[i].IFS[4], sections[i].IFS[5]);
        fclose(fp);
    }

    return 0;
}

bool SpaceFrame::sfCalculate(bool progress_bar = false)
{
    ProgressBar = progress_bar;

    if (sfLCosSin()) // calculate the length, cosine and sine of all rods
    {
        sfPrintError(6);
        return 1;
    }
    else
        printf("Calculating length, cosine and sine succeed!\n");

    if (sfAllocate())
    {
        sfPrintError(15);
        return 0;
    }
    else
        printf("Allocating Variable Bandwith Matrix succeed!\n");

    if (sfBuildTotalStiff()) // build total stiffness matrix
    {
        sfPrintError(2);
        return 1;
    }
    else
        printf("Building total stiffness matrix succeeded!\n");

    if (sfBuildLoadVector(LoadVector)) // build load stiffness vector
    {
        sfPrintError(3);
        return 1;
    }
    else
        printf("Building load vector succeeded!\n");

    if (sfConjugateGradientPar(TotalStiffness, LoadVector, Displacement, 6 * NFRN)) // solve matrix equation
    {
        sfPrintError(4);
        return 1;
    }
    else
        printf("Solving equation succeeded!\n");

    for (int i = 0; i < NOS; i++)
        if (sfInternalForce(i, sections[i].NRS, sections[i].DSB)) // calculate the internal force of each rods
        {
            sfPrintError(5);
            return 1;
        }

    cout << "Outputing data succeed!\n";

    return 0;
}

int main()
{
    FILE *fp = fopen("source&result/testresult.csv", "w");
    SpaceFrame Frame;
    DWORD start, end;
    fprintf(fp, " DOF      , PB off   , PB on    ,\n");

    for (int i = 3; i < 16; i++)
    {
        sfCircularStructure(i, i, i);
        fprintf(fp, "%9d ,", ((i + 1) * (i + 1) * (i + 1) - (i + 1) * (i + 1)) * 6);

        start = GetTickCount();
        Frame.sfInput();
        Frame.sfCalculate();
        Frame.sfOutput();
        end = GetTickCount();
        fprintf(fp, "%9.2f ,", (double)(end - start) / 1000);
        Frame.~SpaceFrame();

        start = GetTickCount();
        Frame.sfInput();
        Frame.sfCalculate(true);
        Frame.sfOutput();
        end = GetTickCount();
        fprintf(fp, "%9.2f ,\n", (double)(end - start) / 1000);
        Frame.~SpaceFrame();
    }

    fclose(fp);
    fp = NULL;
    return 0;
}