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
void KGaaa()
{
	printf("************************************************************************\n");
}
bool sfInput()
{
	FILE *fp = NULL;														//Define the file point
	char *line,*data;														//Define the line string and separated string
	char temporSpace[207729];												//Apply for temporary storage space
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
						case 1 :  if(columnIndex == 2)TNN = atoi(data);	maxcol = 2;						break;
						case 2 :  if(columnIndex == 2)NFIN = atoi(data); NFRN = TNN - NFIN;	maxcol = 2;	break;
						case 3 :  if(columnIndex == 2)NOR = atoi(data);	maxcol = 2;						break;
						case 4 :  if(columnIndex == 2)NOL = atoi(data);	maxcol = 2;						break;
						case 5 :  if(columnIndex == 2)NOS = atoi(data);	maxcol = 2;						break;
						case 6 :  XCN[columnIndex-2] = atof(data); maxcol = TNN;			break;
						case 7 :  YCN[columnIndex-2] = atof(data); maxcol = TNN;			break;
						case 8 :  ZCN[columnIndex-2] = atof(data); maxcol = TNN;			break;
						case 9 :  BNR[columnIndex-2] = atoi(data); maxcol = NOR;			break;
						case 10:  ENR[columnIndex-2] = atoi(data); maxcol = NOR;			break;
						case 11:  ELASTIC[columnIndex-2] = atof(data); maxcol = NOR;		break;
						case 12:  SHEAR[columnIndex-2] = atof(data); maxcol = NOR;			break;
						case 13:  AREA[columnIndex-2] = atof(data); maxcol = NOR;			break;
						case 14:  IMY[columnIndex-2] = atof(data);	maxcol = NOR;			break;
						case 15:  IMZ[columnIndex-2] = atof(data);	maxcol = NOR;			break;
						case 16:  THETA[columnIndex-2] = atof(data); maxcol = NOR;			break;
						case 17:  NRL[columnIndex-2] = atoi(data); maxcol = NOL;			break;
						case 18:  PLI[columnIndex-2] = atoi(data); maxcol = NOL;			break;
						case 19:  KOL[columnIndex-2] = atoi(data); maxcol = NOL;			break;
						case 20:  VOL[columnIndex-2] = atof(data); maxcol = NOL;			break;
						case 21:  DLB[columnIndex-2] = atof(data); maxcol = NOL;			break;
						case 22:  NRS[columnIndex-2] = atoi(data); maxcol = NOS;			break;
						case 23:  DSB[columnIndex-2] = atof(data); maxcol = NOS;			break;
					}
					//--------------------------------录入结束-------------------------------------	
				}
//					if (columnIndex == maxcol) 				 					//Only read the first 'maxcol' columns
//						break;
					data = strtok(NULL, DIVIDE);								//Reset data
			}
			
//			if(columnIndex-1 != maxcol || rowIndex >= 2)
//				return 1;
//			printf("columnIndex %d\n",columnIndex);
//			printf("maxcol %d\n",maxcol);
			rowIndex ++;														//RowIndex steps forward once
			columnIndex = 0;													//Reset columnIndex
		}
		fclose(fp);																//Close the file
		fp = NULL;																//Reset the file point
	}	
}
bool output_and_print()
{
	
	int i,j; //循坏控制変量
	printf("\t\t\t空同剛架结构汁算\n");
	KGaaa();
	
	printf("\t\tTNN =%d\t\t\tNFIN=%d\n\t\tNFRN=%d\t\t\tNOR =%d\n", TNN, NFIN, NFRN, NOR) ;
	printf("\t\tNOL =%d\t\t\tNOS =%d\n" , NOL,NOS);
	KGaaa();
	
	printf("	节点号\tX坐标\t\tY坐标\t\tZ坐标\n");
	for(i=1;i<=TNN; i++)
		printf("%d\t %5.7f\t%5.7f\t85.7f\n",i,XCN[i-1],YCN[i-1],ZCN[i-1]);
	KGaaa();
	
	printf("	杆件号	左节点	右节点	弾性模量E	剪切模量G	截面枳A		慣性矩Jy	慣性矩Jz\n");
	for(i=0; i<NOR;i++)
		printf("	%d\t	%d\t	%d\t %5.0f %5.0f %5.4f %.5f %.5f\n",i+1, BNR[i], ENR[i], ELASTIC[i], SHEAR[i], AREA[i],IMY[i],IMZ[i]);
	printf("截面号\t所在杆件号\t距左端距禽\n");
	for(i=1;i<=NOS;i++)
		printf("%d\t\t%d\t\t %5.7f\n",i,NRS[i-1],DSB[i-1]);
	KGaaa();
	
	printf("\n結果輸出如下: \n");
	KGaaa() ;

	printf("节点号位移X\t位移Y\t位移z\t特角x\t装角Y\t特角z\n");
	for (i=NFIN+1, j=0;i<=TNN;i++,j++)
		printf("%d%5.7f%5.7f%5.7f%5.7f %5.7f %5.7f\n",i,DON[6*j],DON[6*j+1], DON[6*j+2],DON[6*j+3],DON[6*j+4],DON[6*j+5]);
	KGaaa() ;

	printf("截面号	軸力x\t	剪力Y\t剪力z\t	抂矩x\t弯矩Y\t	弯矩 Z\n");
	for (i=0;i<NOS;i++)
		printf("	%d	%5.7f	%5.7f	%5.7f	%5.7f	%5.7f	%5.7f\n",i+1, IFS[6*i], IFS[6*i+1], IFS[6*i+2],IFS[6*i+3],IFS[6*i+4], IFS[6*i+5]);
		
	
}
void writeExcel()
{
	char chy[4]={ 'x' ,'a' ,'h','w' } ;
	int data[4]={ 1 , 3 , 6 ,9	};
	int i , j;
	FILE *fp = NULL ;
	fp = fopen("sfRESULT.csv","w") ;
	fprintf(fp,"TITLE\n") ;
	fprintf(fp,"TNN,%d\nNFIN,%d\nNFRN,%d\nNOR,%d\nNOL,%d\nNOS,%d", TNN, NFIN, NFRN, NOR, NOL, NOS);
	
	//------------NODES-------------------------------------------------
	fprintf(fp,"NODES\n,");
	for(i=0;i<TNN;i++)
		fprintf(fp,"%d,",i+1);
	
	fprintf(fp,"\nCON,");
	for(i=0;i<TNN;i++)
		fprintf(fp,"(%f %f %f),",XCN[i],YCN[i],ZCN[i]);
	
	//------------RODES-------------------------------------------------
	fprintf(fp,"\nRODES,");
	for(i=0;i<NOR;i++)
		fprintf(fp,"%d,",i+1);
	
	fprintf(fp,"\nBNR->ENR,");
	for(i=0;i<NOR;i++)
		fprintf(fp,"p%d->p%d,",BNR[i],ENR[i]);
	
	fprintf(fp,"\nELASTIC,");
	for(i=0;i<NOR;i++)
		fprintf(fp,"%f,",ELASTIC[i]);
	
	fprintf(fp,"\nSHEAR,");
	for(i=0;i<NOR;i++)
		fprintf(fp,"%f,",SHEAR[i]);
	
	fprintf(fp,"\nAREA,");
	for(i=0;i<NOR;i++)
		fprintf(fp,"%f,",AREA[i]);
	
	fprintf(fp,"\nIMY,");
	for(i=0;i<NOR;i++)
		fprintf(fp,"%f,",IMY[i]);
	
	fprintf(fp,"\nIMZ,");
	for(i=0;i<NOR;i++)
		fprintf(fp,"%f,",IMZ[i]);
	
	fprintf(fp,"\nTHETA,");
	for(i=0;i<NOR;i++)
		fprintf(fp,"%f,",THETA[i]);
	
	//------------LOADS-------------------------------------------------
	fprintf(fp,"\nLOADS,");
	for(i=0;i<NOL;i++)
		fprintf(fp,"%d,",i+1);
	
	fprintf(fp,"\nPLI,");
	for(i=0;i<NOL;i++)
		fprintf(fp,"%d,",PLI[i]);
	
	fprintf(fp,"\nNRL,");
	for(i=0;i<NOL;i++)
		fprintf(fp,"%d,",NRL[i]);
	
	fprintf(fp,"\nKOL,");
	for(i=0;i<NOL;i++)
		fprintf(fp,"%d,",KOL[i]);
	
	fprintf(fp,"\nVOL,");
	for(i=0;i<NOL;i++)
		fprintf(fp,"%f,",VOL[i]);
	
	fprintf(fp,"\nDLB,");
	for(i=0;i<NOL;i++)
		fprintf(fp,"%f,",DLB[i]);
	
	//-----------SECTIONS-------------------------------------------------
	fprintf(fp,"\nNOS,");
	for(i=0;i<NOS;i++)
		fprintf(fp,"%d,",i+1);
	
	fprintf(fp,"\nNRS,");
	for(i=0;i<NOS;i++)
		fprintf(fp,"%d,",NRS[i]);
	
	
	fprintf(fp,"\nDSB,");
	for(i=0;i<NOS;i++)
		fprintf(fp,"%f,",DSB[i]);
	
	//-----------RESULTS OF NODES-----------------------------------------
	fprintf(fp,"\nNORN,");
	for(i=0;i<NFRN;i++)
		fprintf(fp,"x%d,y%d,z%d,",i+NFIN,i+NFIN,i+NFIN);
		
	fprintf(fp,"\nDISPLACEMENT,");
	for(i=0;i<NFRN;i++)
		fprintf(fp,"%f,%f,%f,",DON[6*i],DON[6*i+1],DON[6*i+2]);
		
	fprintf(fp,"\nDIVERSION,");
	for(i=0;i<NFRN;i++)
		fprintf(fp,"%f,%f,%f,",DON[6*i+3],DON[6*i+4],DON[6*i+5]);
	
	//-----------RESULTS OF SECTIONS--------------------------------------
	fprintf(fp,"\nNOS,");
	for(i=0;i<NOS;i++)
		fprintf(fp,"x%d,y%d,z%d,",i+NOS+1,i+NOS+1,i+NOS+1);
		
	fprintf(fp,"\nSHEAR FORCE,");
	for(i=0;i<NOS;i++)
		fprintf(fp,"%f,%f,%f,",IFS[6*i],IFS[6*i+1],IFS[6*i+2]);
	
	fprintf(fp,"\nBENDING MOMENT,");
	for(i=0;i<NOS;i++)
		fprintf(fp,"%f,%f,%f,",IFS[6*i+3],IFS[6*i+4],IFS[6*i+5]);
	fclose(fp);
}
					
//	writeExcel()  ;	

int main()
{
	sfInput();
	output_and_print();
	writeExcel()  ;
	return 0;
}
