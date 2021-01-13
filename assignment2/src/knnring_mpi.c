#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <cblas.h>
#include <string.h>
#include <mpi.h>

#include "../inc/knnring.h"


knnresult kNN(double *X, double *Y, int n, int m, int d, int k){

	knnresult k_nearest;

	//initalize the result's variable
	k_nearest.nidx  = (int *)malloc(m*k*sizeof(int));
	k_nearest.ndist = (double *)malloc(m*k*sizeof(double));
	k_nearest.m = m;
	k_nearest.k = k;

  //compute the distance Matrix D
	double *D = compute_D(X,Y,n,m,d,k);

  //select k-nearest neighbors for every point in QUERY
	k_select(D,k_nearest.ndist, k_nearest.nidx,n,m,k);
	return k_nearest;
}


double *compute_D(double *X, double *Y, int n, int m, int d, int k){


	//compute distances using
	//D = (X*X)ee_t - 2XY_t + ee_t(Y*Y)_t
	double *tempD   = (double *)malloc(n*m*sizeof(double));
	double *e_1     = (double *)malloc(d*m*sizeof(double));
	double *XX      = (double *)malloc(n*d*sizeof(double));

	for(int i=0;i<d*m;i++)
		e_1[i]=1;

	for(int i=0;i<n*d;i++)
		XX[i] = X[i]*X[i];

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n,m,d,1,XX,d,e_1,m,0,tempD,m);

	free(e_1);
	free(XX);

	double *e_2 = (double *)malloc(n*d*sizeof(double));
	double *YY = (double *)malloc(m*d*sizeof(double));

	for(int i=0;i<n*d;i++)
		e_2[i]=1;

	for(int i=0;i<m*d;i++)
		YY[i] = Y[i]*Y[i];

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n,m,d,-2,X,d,Y,d,1,tempD,m);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n,m,d,1,e_2,d,YY,d,1,tempD,m);

	for(int i=0;i<n*m;i++)
		tempD[i] = sqrt(fabs(tempD[i]));


	free(e_2);
	free(YY);


	return tempD;
}


void k_select(double *D, double *ndist, int *nidx,int n, int m, int k){
	for(int j=0;j<m;j++){
		for(int i=0;i<k;i++){
			ndist[j*k+i] = D[j+i*m];
			nidx[j*k+i]  = i;
		}
		//sort the initial k neighbors
		quick_sort(ndist,nidx,j*k,k*(j+1)-1);
	}

	//search the rest points to find it there are closer neighbors
	for(int j=0;j<m;j++){
		for(int i=k;i<n;i++){
			if(D[j+i*m] < ndist[k*(j+1)-1]){
				ndist[k*(j+1)-1] = D[j+i*m];
				nidx[k*(j+1)-1] = i;
				quick_sort(ndist,nidx,j*k,(j+1)*k-1);
			}
		}
	}
}


void quick_sort(double *k_dist, int *k_neigh, int low, int high){
	if(low<high){

		//actuall sorting an element
		int pi = partition(k_dist,k_neigh,low,high);

		//recursivelly call quick_sort
		quick_sort(k_dist,k_neigh,low,pi-1);
		quick_sort(k_dist,k_neigh,pi+1,high);
	}
}


int partition(double *k_dist, int *k_neigh,int low, int high){

	double temp_dist;
	int temp_point;
	int checker = 0;

	double pivot = k_dist[high];
	//printf("the pivot is %lf \n",pivot);
	int i = (low-1);
	for(int j=low;j<=high-1;j++){
		if(k_dist[j]<=pivot){
			/*swap the two elements*/

			i=i+1;
			temp_dist  = k_dist[i];
			temp_point = k_neigh[i];

			k_dist[i]  = k_dist[j];
			k_neigh[i] = k_neigh[j];

			k_dist[j]  = temp_dist;
			k_neigh[j] = temp_point;
			checker = 1;
		}
	}

	/* swap the pivot */
	i=i+1;
	temp_dist  = k_dist[i];
	temp_point = k_neigh[i];

	k_dist[i]  = k_dist[high];
	k_neigh[i] = k_neigh[high];

	k_dist[high]  = temp_dist;
	k_neigh[high] = temp_point;


	return i;
}

void create_points(double *X,int size){
	double *temp;
  temp = X;
	srand(time(NULL));
	for(int i=0;i<size;i++)
		temp[i] = (double) 1000*rand()/RAND_MAX;
}


void write_to_file(char *str1, int *nidx, double *ndist, int amount, int k){
  FILE *fp;
  //file for indexes
  fp = fopen(str1,"a+");

  int counter = 1;
  for(int i=0;i<amount;i++){
    fprintf(fp,"%d)%d -- %lf\n",counter,nidx[i],ndist[i]);
    counter ++;
		if((i+1)%k == 0){
      fprintf(fp,"k nearest neighbors for next point \n");
      counter = 1;
    }
  }
	fclose(fp);
}

//######################### V1 FUNCTIONS ############################

knnresult distrAllkNN(double *X, int n, int d, int k){

  //results variables
  knnresult knn_process;
  knnresult knn_temp;

  int p_rank;                            //PROCESS RANK
  int p_size;                            //TOTAL NUMBER OF PROCESSES
  int next,prev;                         //NEXT AND PREVIOUS PROCESSES
  double *receive_buff;                  //BUFFERS FOR COMMUNICATION
  double *Y;                             //QUERY POINTS

  //variables for MPI COMMUNICATION.
  MPI_Status status;
  MPI_Request send_request,receive_request;
  int tag;

	//time variables
	struct timespec ts_start;
  struct timespec ts_end;
	double time;


  //find the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &p_size);

  //find the id of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);

  //define the communication tag
  tag = 1;

  //define the next and the previous process.
  //and create the circular communication
  next = p_rank+1;
  if(next == p_size) next = 0;

  prev = p_rank-1;
  if(prev== -1) prev = p_size-1;

  //initialize Y block
  Y = (double *)malloc(n*d*sizeof(double));

  //in the first iteration every process computes distances between its own points
  memcpy(Y,X,n*d*sizeof(double));

  //initialize receive buffer
  receive_buff = (double *)malloc(n*d*sizeof(double));

  //MASK COMMUNICATION TIME
	//Send and receive point for the next iteration
	clock_gettime(CLOCK_MONOTONIC,&ts_start);
	MPI_Isend(Y,n*d,MPI_DOUBLE,next,tag,MPI_COMM_WORLD, &send_request);
	MPI_Irecv(receive_buff,n*d, MPI_DOUBLE, prev, tag, MPI_COMM_WORLD, &receive_request);

  //find k-nearest neighbors in this set of points.
  knn_process = kNN(X,Y,n,n,d,k);

  for(int i=0;i<n*k;i++) knn_process.nidx[i] = knn_process.nidx[i]+n*p_rank;

  MPI_Wait(&send_request,&status);
  //make sure receive_buffer is ready
  MPI_Wait(&receive_request,&status);
	clock_gettime(CLOCK_MONOTONIC,&ts_end);
	time = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
  printf("communication time%lf\n",time);

  for(int i=0;i<p_size-1;i++){

    memcpy(Y,receive_buff,n*d*sizeof(double));

    //MASK COMMUNICATION TIME
		//send the block to the next node
		clock_gettime(CLOCK_MONOTONIC,&ts_start);
		MPI_Isend(Y,n*d,MPI_DOUBLE,next,tag,MPI_COMM_WORLD, &send_request);
		MPI_Irecv(receive_buff,n*d,MPI_DOUBLE,prev,tag,MPI_COMM_WORLD,&receive_request);

		knn_temp = kNN(Y,X,n,n,d,k);

    //compare k-nearest points
    knn_combine(knn_process, knn_temp,n,k,(p_rank-i-1+p_size)%p_size);

    //make sure Y buffer has been send in order to change it
		MPI_Wait(&send_request,&status);
    //make sure receive_buffer is ready
		MPI_Wait(&receive_request,&status);


		//clock_gettime(CLOCK_MONOTONIC,&ts_end);
		//time = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
	  //printf("communication time and knn compute time%lf\n",time);
  }
  free(Y);
  free(knn_temp.nidx);
  free(knn_temp.ndist);

  for(int i=0;i<knn_process.m;i++) quick_sort(knn_process.ndist,knn_process.nidx,i*k,k*(i+1)-1);
  return knn_process;
}


void knn_combine(knnresult knn_total,knnresult knn_temp, int n, int k, int no_block){
	for(int m=0;m<knn_total.m;m++){
		for(int j=0;j<k;j++){
			if(knn_temp.ndist[m*k+j]>knn_total.ndist[k*(m+1)-1])
				break;
			else{
				knn_total.ndist[k*(m+1)-1] = knn_temp.ndist[m*k+j];
				knn_total.nidx[k*(m+1)-1] = knn_temp.nidx[m*k+j] + no_block*n;

				//sort the new element
				quick_sort(knn_total.ndist,knn_total.nidx,m*k,(m+1)*k-1);
			}
		}
	}
}
