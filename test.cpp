#include <iostream>
#include <mpi.h>
#include <stdint.h>
using namespace std;

int main(int argc,char* argv[]){
    cout<<"hello"<<endl;
    int myid, numprocs;
    int  namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    printf("%d", argc);
    printf("%s||%s", argv[0], argv[1]);
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Get_processor_name(processor_name,&namelen); 
    std::cout<<"Hello World! Process "<<myid<<" of "<<numprocs<<" on "<<processor_name<<std::endl;
    MPI_Finalize();
    return 0;
}