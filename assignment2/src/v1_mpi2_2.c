/***********************************
V1 IMPLEMENTATION USING THE DATASETS
************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <cblas.h>
#include <string.h>
#include <mpi.h>

#include "../inc/knnring.h"
#include "../inc/readlib.h"

int main(int argc,char *argv[]){
  int p_rank;
  int p_size;
  int n;
  int d;
  int k = 10;
  double *X;

  //CHOOSE WHICH DATASET
  char file_path[200] = "../data/ColorHistogram.asc";

  //INITIALIZE MPI COMMUNICATIONS
  //MPI_Init(&argc,&argv);   doesn't need any argument


  MPI_Init(NULL,NULL);

  MPI_Comm_size(MPI_COMM_WORLD, &p_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
  MPI_Status status;                      // STATUS COMMUNICATION
  int tag = 1;                            //COMMUNICATION TAG

  //time variables to measure time
  struct timespec ts_start;
	struct timespec ts_end;
  double time;


  if(p_rank == 0){

    //read all the data
    double *X_total;
    int n_total;
    X_total = read_data(file_path,&n_total,&d);


    //send n,d variables and the points coordinates to other processes.
    n = n_total/p_size;
    X = (double *)malloc(n*d*sizeof(double));
    for(int i=1;i<p_size;i++){
      for(int j=0;j<n*d;j++){
        X[j] = X_total[i*n*d+j];
      }

      //send to every PROCESS
      MPI_Send(&n,1, MPI_INT,    i, tag, MPI_COMM_WORLD);
      MPI_Send(&d,1, MPI_INT,    i, tag, MPI_COMM_WORLD);
      MPI_Send(X,n*d,MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
    }

    //keep the first n points for myself.
    for(int j=0;j<n*d;j++){
      X[j] = X_total[j];
    }


    //simple knn
    //double *Y;
    //Y = (double *)malloc(10*d*sizeof(double));
    //memcpy(Y,X_total,10*d*sizeof(double));
    //clock_gettime(CLOCK_MONOTONIC, &ts_start);
    //knnresult results_simple_knn = kNN(X_total,Y,n_total,10,d,k);
    //clock_gettime(CLOCK_MONOTONIC, &ts_end);
    //time = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
    //printf("TOTAL V0 time: %lf\n",time);
    //char str1[200];
    //snprintf(str1,sizeof(str1),"results/simple_knn_%d_%d_%d.txt",n,d,k);
    //write_to_file(str1,results_simple_knn.nidx,results_simple_knn.ndist,10*k,k);
    //free(Y);
    free(X_total);


    clock_gettime(CLOCK_MONOTONIC, &ts_start);
    knnresult results = distrAllkNN(X,n,d,k);
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    time = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
    printf("TOTAL V1 time: %lf\n",time);

    //WRITE RESULTS TO FILE
    char str[200];
    snprintf(str,sizeof(str),"results/dist_knn_%d_%d_%d.txt",n,d,k);
    write_to_file(str,results.nidx,results.ndist,n*k,results.k);

    for(int p=1;p<p_size;p++){
      MPI_Recv(results.nidx, n*k,    MPI_INT,p,tag,MPI_COMM_WORLD, &status);
      MPI_Recv(results.ndist,n*k, MPI_DOUBLE,p,tag,MPI_COMM_WORLD, &status);
      write_to_file(str,results.nidx,results.ndist,n*k,results.k);
    }
  }
  else{

    MPI_Recv(&n,1,MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
    MPI_Recv(&d,1,MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
    X = (double *)malloc(n*d*sizeof(double));
    MPI_Recv(X,n*d,MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);

    clock_gettime(CLOCK_MONOTONIC, &ts_start);
    knnresult results = distrAllkNN(X,n,d,k);
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    time = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
    //printf("TOTAL V1 time: %lf from process %d\n",time,p_rank);


    MPI_Send(results.nidx, n*k, MPI_INT   ,0,tag,MPI_COMM_WORLD);
    MPI_Send(results.ndist,n*k, MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
  }

  //execution is over
  MPI_Finalize();
  printf("DONE ....=>EXITING\n");
  return 0;
}
