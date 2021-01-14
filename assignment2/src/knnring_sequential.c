#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <cblas.h>
#include <string.h>

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


void write_to_file(char *str1, int *nidx, double *ndist, int amount, int k){
  FILE *fp;
  //file for indexes
  fp = fopen(str1,"w+");

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
