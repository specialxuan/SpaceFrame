#include <iostream>

class SpaceFrame
{
private:
    int TNN;  //total number of nodes
    int NFIN; //number of fixed nodes
    int NFRN; //number of free nodes
    int NOR;  //number of rods
    int NOL;  //number of loads
    int NOS;  //number of sections

    struct Node //parameters of nodes
    {
        double XCN; //X coordinate of nodes
        double YCN; //Y coordinate of nodes
        double ZCN; //Z coordinate of nodes
    };
    Node *nodes; //parameters of nodes

    struct Rod //parameters of nodes
    {
        int ENR;        //the end node number of rods
        int BNR;        //the beginning node number of rods
        double ELASTIC; //elastic modulus
        double SHEAR;   //shear modulus
        double AREA;    //area
        double IMY;     //inertia moment of Y axis
        double IMZ;     //inertia moment of Z axis
        double THETA;   //theta the deflection angle of main inertia axis
        double LCS[4];  //the length, sine and cosine of rods
        double RFE[6];  //the reaction force of the end node
    };
    Rod *rods; //parameters of nodes

    struct Load //parameters of loads
    {
        int NRL;    //the number of rods with load
        int PLI;    //the plane of the load's in
        int KOL;    //the kind of load
        double VOL; //the value of load
        double DLB; //the distance between load and the beginning node
    };
    Load *loads; //parameters of loads

    struct Section //parameters of sections
    {
        int NRS;       //the number of rod with section
        double DSB;    //the distance between section and the beginning node
        double IFS[6]; //the internal force in the section
    };
    Section *sections; //parameters of sections

    double *TotalStiffness; //matrix of total stiffness
    double *LoadVector;     //vector of loads
    double *Displacement;   //displacement of nodes

public:
    //initialize SpaceFrame
    SpaceFrame();
    //delete SpaceFrame
    ~SpaceFrame();

    //input
    bool sfInput();
    //output
    bool sfOutput();
};

int main()
{

    return 0;
}

SpaceFrame::SpaceFrame()
{
    TNN = 0;
    NFIN = 0;
    NFRN = 0;
    NOR = 0;
    NOL = 0;
    NOS = 0;

    nodes = NULL;
    rods = NULL;
    loads = NULL;
    sections = NULL;

    TotalStiffness = NULL;
    LoadVector = NULL;
    Displacement = NULL;
}

SpaceFrame::~SpaceFrame()
{
    delete nodes;
    delete rods;
    delete loads;
    delete sections;
    delete TotalStiffness;
    delete LoadVector;
    delete Displacement;
}

bool SpaceFrame::sfInput()
{
    return 0;
}

bool SpaceFrame::sfOutput()
{
    return 0;
}
