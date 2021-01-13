/***********************************
V2 IMPLEMENTATION USING RANDOM POINTS
************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <cblas.h>
#include <string.h>
#include <mpi.h>

#include "../inc/knnring.h"
#include "../inc/vptring.h"

#define B 5

int main(int argc, char *argv[]){

  int p_rank;
  int p_size;
  int n;
  int d;
  int k;
  double *X;

  //INITIALIZE MPI COMMUNICATIONS
  MPI_Init(&argc,&argv);
  if(argc<4){
    printf("ERROR: GIVE CORRECT NUMBERS OF ARGUMENTS n,d,k\n");
    exit(1);
  }else{
    n = atoi(argv[1]);
    d = atoi(argv[2]);
    k = atoi(argv[3]);
  }
  MPI_Comm_size(MPI_COMM_WORLD, &p_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
  MPI_Status status;                      // STATUS COMMUNICATION
  int tag = 1;                            //COMMUNICATION TAG

  struct timespec ts_start;
  struct timespec ts_end;
  double time;

  if(p_rank == 0){
    //create all the points and distribut them
    X = (double *)malloc(n*d*sizeof(double));

    for(int i=1;i<p_size;i++){
      create_points(X,n*d);
      MPI_Send(X,n*d,MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
    }

    printf("edw ftanei ???\n");

    //compute points for myself
    create_points(X,n*d);

    clock_gettime(CLOCK_MONOTONIC, &ts_start);
    knnresult results = distrV2(X,n,d,k,B);
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    time = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
    printf("TOTAL V2 time: %lf\n",time,p_rank);

    //write results to file
    //char str[200];
    //snprintf(str,sizeof(str),"results/dist_v2_%d_%d_%d.txt",n,d,k);
    //write_to_file(str,results.nidx,results.ndist,n*k,results.k);

    for(int p=1;p<p_size;p++){
      MPI_Recv(results.nidx, n*k,    MPI_INT,p,tag,MPI_COMM_WORLD, &status);
      MPI_Recv(results.ndist,n*k, MPI_DOUBLE,p,tag,MPI_COMM_WORLD, &status);
      write_to_file(str,results.nidx,results.ndist,n*k,results.k);
    }

  }
  else{

    X = (double *)malloc(n*d*sizeof(double));
    MPI_Recv(X,n*d,MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);

    clock_gettime(CLOCK_MONOTONIC, &ts_start);
    knnresult results = distrV2(X,n,d,k,B);
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    time = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
    //printf("TOTAL V2 time: %lf from process %d\n",time,p_rank);

    //send back the results
    MPI_Send(results.nidx, n*k, MPI_INT   ,0,tag,MPI_COMM_WORLD);
    MPI_Send(results.ndist,n*k, MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
  }

  //execution is over
  MPI_Finalize();
  printf("DONE ....=>EXITING\n");
}
