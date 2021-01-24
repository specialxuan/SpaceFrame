// #include <iostream>
// #include <math.h>
// using namespace std;

// class SpaceFrame
// {
// private:
//     int TNN;  //total number of nodes
//     int NFIN; //number of fixed nodes
//     int NFRN; //number of free nodes
//     int NOR;  //number of rods
//     int NOL;  //number of loads
//     int NOS;  //number of sections

//     struct Node //parameters of nodes
//     {
//         double XCN; //X coordinate of nodes
//         double YCN; //Y coordinate of nodes
//         double ZCN; //Z coordinate of nodes
//     };
//     Node *nodes; //parameters of nodes

//     struct Rod //parameters of nodes
//     {
//         int ENR;        //the end node number of rods
//         int BNR;        //the beginning node number of rods
//         double ELASTIC; //elastic modulus
//         double SHEAR;   //shear modulus
//         double AREA;    //area
//         double IMY;     //inertia moment of Y axis
//         double IMZ;     //inertia moment of Z axis
//         double THETA;   //theta the deflection angle of main inertia axis
//         double LCS[4];  //the length, sine and cosine of rods
//         double RFE[6];  //the reaction force of the end node
//     };
//     Rod *rods; //parameters of nodes

//     struct Load //parameters of loads
//     {
//         int NRL;    //the number of rods with load
//         int PLI;    //the plane of the load's in
//         int KOL;    //the kind of load
//         double VOL; //the value of load
//         double DLB; //the distance between load and the beginning node
//     };
//     Load *loads; //parameters of loads

//     struct Section //parameters of sections
//     {
//         int NRS;       //the number of rod with section
//         double DSB;    //the distance between section and the beginning node
//         double IFS[6]; //the internal force in the section
//     };
//     Section *sections; //parameters of sections

//     double *TotalStiffness; //matrix of total stiffness
//     double *LoadVector;     //vector of loads
//     double *Displacement;   //displacement of nodes

//     int *IV;     //the location of diagonal element
//     int NSI;     //upper limit
//     int MAXIBDW; //half bandwidth

//     int sfSolve; //0 is conjugate gradient par, 1 is conjugate gradient, 2 is cholesky
//     double EPS = 1e-15;

//     //calculate the length sine and cosine of rods
//     bool sfLCosSin()
//     {
//         for (int k = 0; k < NOR; k++)
//         {
//             int i = rods[k].BNR - 1, j = rods[k].ENR - 1; //index of beginning and end nodes of rods
//             rods[k].LCS[1] = nodes[j].XCN - nodes[i].XCN;
//             rods[k].LCS[2] = nodes[j].YCN - nodes[i].YCN;
//             rods[k].LCS[3] = nodes[j].ZCN - nodes[i].ZCN;
//             rods[k].LCS[0] = sqrt(rods[k].LCS[1] * rods[k].LCS[1] + rods[k].LCS[2] * rods[k].LCS[2] + rods[k].LCS[3] * rods[k].LCS[3]);
//             if (rods[k].LCS[0] < EPS) //if the length of rod is too small, then return error
//             {
//                 sfPrintError(8);
//                 return 1;
//             }
//             rods[k].LCS[1] = rods[k].LCS[1] / rods[k].LCS[0];
//             rods[k].LCS[2] = rods[k].LCS[2] / rods[k].LCS[0];
//             rods[k].LCS[3] = rods[k].LCS[3] / rods[k].LCS[0];
//         }

//         return 0;
//     }
//     //build variable bandwith matrix
//     bool sfVarBandwith()
//     {
//         int it = 0, mm = 0;
//         int *peribdw = new int[TNN](); //bandwidth per line in total stiffness matrix
//         IV = new int[6 * NFRN]();
//         for (int i = 0; i < NOR; i++) //for each rod
//         {
//             if (rods[i].BNR > NFIN)
//             {
//                 mm = rods[i].ENR - rods[i].BNR; //bandwidth is end number minus begin number
//                 if (mm > peribdw[rods[i].ENR - 1])
//                     peribdw[rods[i].ENR - 1] = mm; //find the maximum bandwith per line
//             }
//         }
//         for (int i = NFIN; i < TNN; i++) //for each line in total stiffness matrix
//         {
//             if (peribdw[i] > MAXIBDW) //find maxim 
//                 MAXIBDW = peribdw[i];
//             for (int j = 1; j <= 6; j++)
//             {
//                 it = it + 1;
//                 if (it == 1)
//                     IV[it - 1] = 6 * peribdw[i] + j;
//                 else
//                     IV[it - 1] = IV[it - 2] + 6 * peribdw[i] + j;
//             }
//         }
//         MAXIBDW = 6 * MAXIBDW + 5;
//         NSI = IV[6 * NFRN - 1];
//         delete[] peribdw;
//         TotalStiffness = new double[NSI];

//         return 0;
//     }
//     //total stiffness matrix
//     inline double& sfTotalStiff(int i, int j)
//     {
//         double &tmp = TotalStiffness[IV[i] - i + j - 1];
//         return tmp;
//     }
//     //build total stiffness matrix
//     bool sfBuildTotalStiff()
//     {
//         sfVarBandwith();


//     }
//     //build unit stiffness matrix
//     bool sfBuildUnitStiff();
//     //build local stiffness matrix
//     bool sfBuildLocalStiff();
//     //build transpose matrix
//     bool sfBuildTrans();
//     //build load vector
//     bool sfBuildLoadVector();
//     //calculate reaction force
//     bool sfReactionForce();
//     //solve equation of matrix by cholesky
//     bool sfCholesky();
//     //solve equation of matrix by conjugate gradient
//     bool sfConjugateGradient();
//     //solve equation of matrix by conjugate gradient parallel
//     bool sfConjugateGradientPar();
//     //calculate internal force of rods
//     bool sfInternalForce();
//     //calculate internal force of cantilever beam
//     bool sfCtlInternalForce();
//     //calculate internal force of displacement
//     bool sfDisplacementForce();
//     //print"----------------------------------------"
//     bool sfPrintLine()
//     {
//         cout << "--------------------------------------------------------------------------\n";
//         return 0;
//     }
//     //print"****************************************"
//     bool sfPrintLine2()
//     {
//         cout << "**************************************************************************\n";
//         return 0;
//     }
//     //print error
//     bool sfPrintError(int error);

// public:
//     //initialize SpaceFrame
//     SpaceFrame();
//     //delete SpaceFrame
//     ~SpaceFrame();

//     //input
//     bool sfInput();
//     //calculate
//     bool sfCalculate();
//     //output
//     bool sfOutput();
//     //set solving method
//     bool sfSetSolve(int);
// };

// int main()
// {
//     SpaceFrame Frame;
//     Frame.sfInput();
//     Frame.sfCalculate();
//     Frame.sfOutput();

//     return 0;
// }

// SpaceFrame::SpaceFrame()
// {
//     TNN = 0;
//     NFIN = 0;
//     NFRN = 0;
//     NOR = 0;
//     NOL = 0;
//     NOS = 0;

//     nodes = NULL;
//     rods = NULL;
//     loads = NULL;
//     sections = NULL;

//     TotalStiffness = NULL;
//     LoadVector = NULL;
//     Displacement = NULL;

//     IV = NULL;
//     NSI = 0;
//     MAXIBDW = 0;

//     sfSolve = 0;
// }

// //-------------------------------------
// //--------- public functions ----------
// //-------------------------------------

// SpaceFrame::~SpaceFrame()
// {
//     delete[] nodes;
//     delete[] rods;
//     delete[] loads;
//     delete[] sections;

//     delete[] TotalStiffness;
//     delete[] LoadVector;
//     delete[] Displacement;
// }

// bool SpaceFrame::sfInput()
// {
//     return 0;
// }

// bool SpaceFrame::sfCalculate()
// {
//     return 0;
// }

// bool SpaceFrame::sfOutput()
// {
//     return 0;
// }

// bool SpaceFrame::sfSetSolve(int solve)
// {
//     if (solve == 0 || solve == 1 || solve == 2)
//     {
//         sfSolve = solve;
//     }
//     else
//     {
//         cout << "Seting solving method failed!\n";
//     }

//     return 0;
// }
