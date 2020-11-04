#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>
#define EPS 1e-7

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

bool sfInput()
{
	FILE *fp = NULL;														//Define the file point
	char *line,*data;														//Define the line string and separated string
	char temporSpace[2077296];												//Apply for temporary storage space
//	char *temporSpace = (char *)malloc(10000 * sizeof(char));
//	memset(temporSpace,0,10000 * sizeof(char));
	int rowIndex = 0;														//Reset the number of rows to zero
	int columnIndex = 0;													//Reset the number of columns to zero
	int maxcol=0;                                                           //Define the maxium number of columns
	const char DIVIDE[] = ",";												//Set the separater as a ','
	if ((fp = fopen("sf.csv","at+")) != NULL)								//Start the process when the file opens successfully
	{
		fseek(fp, 0L, SEEK_SET); 											//Locate file point to the first line
		while ((line = fgets(temporSpace, sizeof(temporSpace), fp))!=NULL)	//The loop continues when the end of the file is not read
		{
			data = strtok(line, DIVIDE);									//Split strings with a ',' as a separator
			while (data != NULL)											//Read the data of each row
			{
				if (strcmp(data, "END") == 0)								//When the keyword 'END' is read, the reading process will be shut down
				{
					return 0;
				}
				//----------------分配内存-------------------------------------------------				
				if(rowIndex == 5)											//Request space for multiple variable when loops to the 5th line
				{
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
    				memset(IMZ, 0, NOR * sizeof(double));
    				THETA = (double *)malloc(NOR * sizeof(double));
    				memset(THETA, 0, NOR * sizeof(double));
    				NRL = (int *)malloc(NOL * sizeof(int));
    				memset(NRL, 0, NOL * sizeof(int));
    				PLI = (int *)malloc(NOL * sizeof(int));
    				memset(PLI, 0, NOL * sizeof(int));
    				KOL = (int *)malloc(NOL * sizeof(int));
    				memset(KOL, 0, NOL * sizeof(int));
    				VOL = (double *)malloc(NOL * sizeof(double));
    				memset(NRL, 0, NOL * sizeof(double));
    				DLB = (double *)malloc(NOL * sizeof(double));
    				memset(DLB, 0, NOL * sizeof(double));
    				NRS = (int *)malloc(NOS * sizeof(int));
    				memset(NRS, 0, NOS * sizeof(int));
    				DSB = (double *)malloc(NOS * sizeof(double));
    				memset(NRL, 0, NOS * sizeof(double));
    				DON = (double *)malloc(3 * NFRN * sizeof(double));
    				memset(DON, 0, 3 * NFRN * sizeof(double));
    				IFS = (double *)malloc(3 * NOS * sizeof(double));
    				memset(IFS, 0, 3 * NOS * sizeof(double));
    				RFE = (double *)malloc(6 * NOR * sizeof(double));
    				memset(NRL, 0, 6 * NOR * sizeof(double));
				}
				//-------------------分配结束-------------------------------------------------
				if(columnIndex++ != 0)											//Skip the saving of the first column
				{
					//--------------------------------数据输入-------------------------------------
					switch(rowIndex)											//Store variables of each column in different ways
					{
						case 0 :  break;
						case 1 :  if(columnIndex == 2)TNN = atoi(data);							break;
						case 2 :  if(columnIndex == 2)NFIN = atoi(data); NFRN = TNN - NFIN;		break;
						case 3 :  if(columnIndex == 2)NOR = atoi(data);							break;
						case 4 :  if(columnIndex == 2)NOL = atoi(data);							break;
						case 5 :  if(columnIndex == 2)NOS = atoi(data);							break;
						case 6 :  XCN[columnIndex-2] = atof(data); maxcol = TNN;			break;
						case 7 :  YCN[columnIndex-2] = atof(data);			break;
						case 8 :  ZCN[columnIndex-2] = atof(data);			break;
						case 9 :  BNR[columnIndex-2] = atoi(data);			break;
						case 10:  ENR[columnIndex-2] = atoi(data);			break;
						case 11:  ELASTIC[columnIndex-2] = atof(data);		break;
						case 12:  SHEAR[columnIndex-2] = atof(data);		break;
						case 13:  AREA[columnIndex-2] = atof(data);			break;
						case 14:  IMY[columnIndex-2] = atof(data);			break;
						case 15:  IMZ[columnIndex-2] = atof(data);			break;
						case 16:  THETA[columnIndex-2] = atof(data);		break;
						case 17:  NRL[columnIndex-2] = atoi(data);			break;
						case 18:  PLI[columnIndex-2] = atoi(data);			break;
						case 19:  KOL[columnIndex-2] = atoi(data);			break;
						case 20:  VOL[columnIndex-2] = atof(data);			break;
						case 21:  DLB[columnIndex-2] = atof(data);			break;
						case 22:  NRS[columnIndex-2] = atoi(data);			break;
						case 23:  DSB[columnIndex-2] = atof(data);			break;
					}
					//--------------------------------录入结束-------------------------------------	
				}
//					if (columnIndex == maxcol) 				 					//Only read the first 'maxcol' columns
//						break;
					data = strtok(NULL, DIVIDE);								//Reset data
			}
			rowIndex ++;														//RowIndex steps forward once
			if(columnIndex-1 != maxcol)
				return 1;
			columnIndex = 0;													//Reset columnIndex
		}
		fclose(fp);																//Close the file
		fp = NULL;																//Reset the file point
	}	
}
int main()
{
	int i=0;
	sfInput();
	printf("%d\n",TNN);
	printf("%d\n",NFIN);
	printf("%d\n",NOR);
	printf("%d\n",NOL);
	printf("%d\nXCN ",NOS);
	for(i=0;i<TNN;i++)
	{
		printf("%f ",XCN[i]);
	}
	printf("\nYCN ");
	for(i=0;i<TNN;i++)
	{
		printf("%f ",YCN[i]);
	}
	printf("\nZCN ");
	for(i=0;i<TNN;i++)
	{
		printf("%f ",ZCN[i]);
	}
	printf("\nBNR ");
	for(i=0;i<NOR;i++)
	{
		printf("%d ",BNR[i]);
	}
	printf("\nENR ");
	for(i=0;i<NOR;i++)
	{
		printf("%d ",ENR[i]);
	}
	printf("\nELASTIC ");
	for(i=0;i<NOR;i++)
	{
		printf("%f ",ELASTIC[i]);
	}
	printf("\nSHEAR ");
	for(i=0;i<NOR;i++)
	{
		printf("%f ",SHEAR[i]);
	}
	printf("\nAREA ");
	for(i=0;i<NOR;i++)
	{
		printf("%f ",AREA[i]);
	}
	printf("\nIMY ");
	for(i=0;i<NOR;i++)
	{
		printf("%.11f ",IMY[i]);
	}
	printf("\nIMZ ");
	for(i=0;i<NOR;i++)
	{
		printf("%.11f ",IMZ[i]);
	}
	printf("\nTHETA ");
	for(i=0;i<NOR;i++)
	{
		printf("%f ",THETA[i]);
	}
	printf("\nNRL ");
	for(i=0;i<NOL;i++)
	{
		printf("%d ",NRL[i]);
	}
	printf("\nPLI ");
	for(i=0;i<NOL;i++)
	{
		printf("%d ",PLI[i]);
	}
	printf("\nKOL ");
	for(i=0;i<NOL;i++)
	{
		printf("%d ",KOL[i]);
	}
	printf("\nVOL ");
	for(i=0;i<NOL;i++)
	{
		printf("%f ",VOL[i]);
	}
	printf("\nDLB ");
	for(i=0;i<NOL;i++)
	{
		printf("%f ",DLB[i]);
	}
	printf("\nNRS ");
	for(i=0;i<NOS;i++)
	{
		printf("%d ",NRS[i]);
	}
	printf("\nDSB ");
	for(i=0;i<NOS;i++)
	{
		printf("%f ",DSB[i]);
	}
	printf("\n");
	return 0;
}
	