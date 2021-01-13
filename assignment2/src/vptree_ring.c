#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "../inc/knnring.h"
#include "../inc/vptring.h"


tree_struct create_tree(double *X, int *nidx,int n, int d, int B){

  //variables mu_points,step,tag,tree_counter are for debugging
  //and need to get rid of them
  tree_struct tree;

  tree.indexes = (int *)malloc(n*sizeof(int));
  tree.mu = (double *)malloc(n*sizeof(double));
  for(int i=0;i<n;i++)tree.mu[i] = 0;
  tree.mu_points = (double *)malloc(n*sizeof(double));
  tree.coords = X;

  tree.tree_counter = (int *)malloc(sizeof(int));
  tree.tree_counter = 0;
  tree.leaf_counter = (int *)malloc(sizeof(int));
  *tree.leaf_counter = n-1;

  tree.tag = (double *)malloc(sizeof(double));
  *tree.tag = 0;
  tree.step = -(double)1/(2*n);

  int height =0;
  create_subtree(tree,X,nidx,n,d,B,height,n);

  //write the tree to a file
  //char str3[200] = "indexandmu.txt";
  //write_to_file_3(str3,tree.indexes,tree.mu,tree.mu_points,n);
  return tree;
}


void create_subtree(tree_struct tree,double *X, int *nidx, int n, int d, int B, int height, int n_total){
  if(n == 0 ){
    return ;
  }
  //inside a leaf
  else if(n == 1){
    //we have overflow
    if(height>=n_total){
      while(tree.mu[*tree.leaf_counter] !=0){
        if(*tree.leaf_counter<0){
          *tree.leaf_counter = n_total-1;
        }
        *tree.leaf_counter = *tree.leaf_counter - 1;
      }

      tree.indexes[*tree.leaf_counter] = nidx[0];
      tree.mu[*tree.leaf_counter] = -MAX_POINTS;
	  tree.mu_points[*tree.leaf_counter] = 1;
	  double fractpart,intpart;
	  //right leaf
	  if(height%2 == 0){
			fractpart = modf(tree.mu[(height-2)/2],&intpart);
			tree.mu[(height-2)/2] = intpart-(double)((double)(*tree.leaf_counter)/(double)(MAX_POINTS));
			}
	  //left leaf
	  else {
			fractpart = modf(tree.mu[(height-1)/2],&intpart);
			tree.mu[(height-1)/2] = (double)(-1*(*tree.leaf_counter))+fractpart;
			}
	  *tree.leaf_counter = *tree.leaf_counter - 1;
    }
    //dont have overflow
    else{
      tree.indexes[height] = nidx[0];
      tree.mu[height] = -MAX_POINTS;
      tree.mu_points[height] = 1;
    }
  }
  else{
    //new vantage point
    tree.indexes[height] = nidx[n-1];


    //calculate_distances
    //and find mu
    int vantage_point = nidx[n-1];
    double *distances = (double *)malloc((n-1)*sizeof(double));
    double *distances_copy = (double *)malloc((n-1)*sizeof(double));
    calculate_distances(distances,vantage_point,n,nidx,X,d);
    memcpy(distances_copy,distances,(n-1)*sizeof(double));

    double mu;
    int *nidx_left;
    int *nidx_right;
    int len_left;
    int len_right;

    if((n-1)%2==0){
      mu = quick_select(distances_copy,0,n-2,(n-1)/2-1);
      mu += quick_select(distances_copy,0,n-2,(n-1)/2);
      mu = mu*0.5;
      nidx_left  = (int *)malloc(((n-1)/2)*sizeof(int));
      nidx_right = (int *)malloc(((n-1)/2)*sizeof(int));
      len_right = (n-1)/2;
      len_left = (n-1)/2;
    }
    //perito plithos.
    else{
      mu = quick_select(distances_copy,0,n-2,(n-2)/2);
      nidx_left  = (int *)malloc((1+(n-2)/2)*sizeof(int));
      nidx_right = (int *)malloc(((n-2)/2)*sizeof(int));
      len_right = (n-2)/2;
      len_left = (n-2)/2+1;
    }
    free(distances_copy);

		if(n<=B){ //we count it as a leaf
			if(n == 2) tree.mu[height] = (double)(-1*(2*height+1)) - (double)((double)1/(double)(10*MAX_POINTS));
			else tree.mu[height]=(double)(-1*(2*height+1))-(double)((double)(2*height+2)/(double)MAX_POINTS);
		}
		else tree.mu[height] = mu;

    tree.mu_points[height] = n;


    int counter_right=0;
    int counter_left=0;
    int counter_mu=0;

    for(int i=0;i<n-1;i++){
      if(distances[i]==mu)counter_mu++;

      if(counter_right == len_right){
        nidx_left[counter_left]=nidx[i];
        counter_left++;
      }
      else if(counter_left == len_left){
        nidx_right[counter_right]=nidx[i];
        counter_right++;
      }
      else if(distances[i]<=mu){
        nidx_left[counter_left] = nidx[i];
        counter_left++;
      }else{
        nidx_right[counter_right] = nidx[i];
        counter_right++;
      }
    }
    free(distances);


    create_subtree(tree,X,nidx_left,len_left,d,B,2*height+1,n_total);
    free(nidx_left);

    create_subtree(tree,X,nidx_right,len_right,d,B,2*height+2,n_total);
    free(nidx_right);
  }
}


void calculate_distances(double *distances, int vantage_point, int n, int *nidx, double *X, int d){
  double temp;

  for(int i=0;i<n-1;i++){
    temp = 0;
    for(int j=0;j<d;j++){
      temp = temp + pow(X[vantage_point*d+j]-X[nidx[i]*d+j],2);
    }
    distances[i] = sqrt(fabs(temp));
  }
  /*char str2[200]= "distances.txt";
  FILE *fp;
  fp = fopen(str2,"a+");
  fprintf(fp,"VANTAGE_POINT:%d -(%lf,%lf)\n",vantage_point,X[vantage_point*d],X[vantage_point*d+1]);
  for(int i =0;i<n-1;i++)fprintf(fp,"%d(%lf,%lf) -(%lf)\n",nidx[i],X[nidx[i]*d],X[d*nidx[i]+1],distances[i]);
  fclose(fp);*/
}


double quick_select(double *arr, int left, int right, int pos){
	if(pos>=0 && pos<=right-left){
		int index = partition_1(arr,left,right);

		if(index-left == pos){
			return arr[index];
		}

		if(index-left>pos)
			return quick_select(arr,left,index-1,pos);

		if(index-left<pos)
			return quick_select(arr,index+1,right,pos-index+left-1);
	}
}


int partition_1(double *arr, int left, int right){
	double x = arr[right];
	int i = left;
	for(int j = left;j<=right-1;j++){
		if(arr[j] <=x){
			swap(&arr[i],&arr[j]);
			i++;
		}
	}
	swap(&arr[i],&arr[right]);
	return i;
}


void swap(double *a, double *b){
	double temp = *a;
	*a = *b;
	*b = temp;
}

//function to write the vantage point.
//for debbuging purposes.
void write_to_file_3(char *str, int *array,double *array2,double *array3, int amount){
  FILE *fp;
  fp = fopen(str,"w+");
  for(int i=0;i<amount;i++){
    fprintf(fp,"%d -- %lf, n = %lf \n",array[i],array2[i],array3[i]);
  }
	fclose(fp);
}


//##########################SEARCH FUNCTIONS##############################
knnresult search_k_nearest(double *X, int *indexes, double *mu, double *Y, int n, int m, int d, int k){
	
	knnresult knn;

  int    *num_neigh = (int *)malloc(sizeof(int));
  int    *nidx = (int *)malloc(m*k*sizeof(int));
  double *ndist = (double *)malloc(m*k*sizeof(double));
  int height;


  //time variables to measure times
  struct timespec ts_start;
  struct timespec ts_end;
  double time;
  double mean_time=0;


  //search for every element
  for (int i=0;i<m;i++){
		*num_neigh = 0;
		height = 0;
    //initialize distances to INFINITY
    for(int j=0;j<k;j++){
      ndist[i*k+j] = INFINITY;
      nidx[i*k+j] = -1;
    }
    clock_gettime(CLOCK_MONOTONIC, &ts_start);
	search_tree_1(X,indexes,mu,Y,nidx,ndist,n,i,d,k,height,num_neigh);
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    time = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
    mean_time += time;
	}


  //printf("TOTAL TIME ONLY FOR SEARCHING THE TREE IS %lf\n",mean_time);
  //mean_time = mean_time/(double)m;
  //printf("MEAN TIME TO SEARCH ONE ELEMENT IS %lf\n",mean_time);

	//write back the results.
	/*for(int i=0;i<m*k;i++){
		if(i%k==0)printf("\n");
		printf("|index:%d,coords(%lf,%lf),dist:%lf|",nidx[i],X[d*nidx[i]],X[d*nidx[i]+1],ndist[i]);
	}*/
	knn.nidx  = nidx;
	knn.ndist = ndist;
	knn.m = m;
	knn.k = k;
	return knn;
}


void search_tree_1(double *X, int *indexes, double *mu, double *Y, int *nidx,double *ndist, int n, int point, int d, int k, int height, int *num_neigh){

  if(mu[height]<0){
    check_leaf(X,indexes,mu,Y,nidx,ndist,n,point,d,k,height,num_neigh);
    return;
  }
  else{
    double distance = 0;
    for(int j=0;j<d;j++)distance += pow(Y[point*d+j]-X[d*indexes[height]+j],2);
	
	distance = sqrt(fabs(distance));

    insert_point_vpt(nidx,ndist,point,k,num_neigh,indexes[height],distance);
    double tau = ndist[(1+point)*k-1];

    if(distance<mu[height]){
      if(distance-tau<=mu[height])search_tree_1(X,indexes,mu,Y,nidx,ndist,n,point,d,k,2*height+1,num_neigh);  //search the left child
      tau = ndist[(1+point)*k-1];
      if(distance+tau>=mu[height])search_tree_1(X,indexes,mu,Y,nidx,ndist,n,point,d,k,2*height+2,num_neigh);  //search the right child
    }
    else{
      if(distance+tau>=mu[height])search_tree_1(X,indexes,mu,Y,nidx,ndist,n,point,d,k,2*height+2,num_neigh);  //search the right child
      tau = ndist[(1+point)*k-1];
      if(distance-tau<=mu[height])search_tree_1(X,indexes,mu,Y,nidx,ndist,n,point,d,k,2*height+1,num_neigh);  //search the left child
    }
  }
}


void check_leaf(double *X,int *indexes,double *mu,double *Y,int *nidx,double *ndist,int n,int point,int d,int k,int height,int *num_neigh){
	double distance;
	//its a leaf
	if(mu[height] == -MAX_POINTS){
		distance = 0;
		for(int j=0;j<d;j++)distance += pow(Y[point*d+j]-X[d*indexes[height]+j],2);

		distance = sqrt(fabs(distance));
		insert_point_vpt(nidx,ndist,point,k,num_neigh,indexes[height],distance);
		return;
	}
	else if(mu[height]<0){
		distance = 0;
		for(int j=0;j<d;j++)distance += pow(Y[point*d+j]-X[d*indexes[height]+j],2);
		
		distance = sqrt(fabs(distance));
		insert_point_vpt(nidx,ndist,point,k,num_neigh,indexes[height],distance);


		double fractpart,intpart;
		fractpart = modf(mu[height],&intpart);
		check_leaf(X,indexes,mu,Y,nidx,ndist,n,point,d,k,abs(intpart),num_neigh);

		if(round(fractpart*(double)(10*MAX_POINTS)) != -1.0 ){
			check_leaf(X,indexes,mu,Y,nidx,ndist,n,point,d,k,abs((int)round(fractpart*(double)MAX_POINTS)),num_neigh);
		}
	}
}


void insert_point_vpt(int *nidx, double *ndist,int point, int k, int *num_neigh, int new_point, double dist){
	if(*num_neigh<k){
		ndist[*num_neigh+point*k]= dist;
		nidx[*num_neigh+point*k] = new_point;
		quick_sort(ndist,nidx,point*k,*num_neigh+point*k);
		(*num_neigh)++;
	}
	else{
		if(dist>ndist[k*(point+1)-1]){
			return;
		}
		ndist[k*(point+1)-1] = dist;
		nidx[k*(point+1)-1]  = new_point;
		quick_sort(ndist,nidx,point*k,k*(point+1)-1);
	}
}


//######################### MPI FUNCTIONS ###################################
knnresult distrV2(double *X, int n,int d, int k, int B){

  //time variables to measure times
  struct timespec ts_start;
  struct timespec ts_end;
  double time;

  //variables to store results
  knnresult knn;
  knnresult knn_temp;

  int    p_rank;
  int    p_size;
  int    next,prev;
  int    *nidx;

  //variables to store the incoming tree
  int    *receive_buff_indexes;
  double *receive_buff_mu;
  double *receive_buff_coords;

  //store my points
  double *Y;

  //MPI COMMUNICATION
  MPI_Status  status[6];
  MPI_Request requests[6];
  int tag;
  tag = 1;

  MPI_Comm_size(MPI_COMM_WORLD, &p_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);

  //create the ring and find which process is next and previous
  next = p_rank+1;
  if(next == p_size) next = 0;

  prev = p_rank-1;
  if(prev == -1) prev = p_size-1;

  //create tree.
  nidx= (int *)malloc(n*sizeof(int));
  for(int i=0;i<n;i++)nidx[i]=i;

  clock_gettime(CLOCK_MONOTONIC, &ts_start);
  tree_struct vpt_tree = create_tree(X,nidx,n,d,B);
  clock_gettime(CLOCK_MONOTONIC, &ts_end);
  time = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
  printf("create tree time: %lf p_id = %d\n",time, p_rank);

  //every process must keep its points to walk the incoming trees
  Y = (double *)malloc(n*d*sizeof(double));
  memcpy(Y,vpt_tree.coords,n*d*sizeof(double));

  //initialize receive buffers for passing vantage point trees.
  receive_buff_indexes = (int *)malloc(n*sizeof(int));
  receive_buff_mu      = (double *)malloc(n*sizeof(double));
  receive_buff_coords  = (double *)malloc(n*d*sizeof(double));


  MPI_Isend(vpt_tree.indexes, n,MPI_INT,   next,tag,MPI_COMM_WORLD, &requests[3]);
  MPI_Isend(vpt_tree.mu     , n,MPI_DOUBLE,next,tag,MPI_COMM_WORLD, &requests[4]);
  MPI_Isend(vpt_tree.coords,n*d,MPI_DOUBLE,next,tag,MPI_COMM_WORLD, &requests[5]);

  MPI_Irecv(receive_buff_indexes, n,MPI_INT   ,prev,tag,MPI_COMM_WORLD,&requests[0]);
  MPI_Irecv(receive_buff_mu     , n,MPI_DOUBLE,prev,tag,MPI_COMM_WORLD,&requests[1]);
  MPI_Irecv(receive_buff_coords,n*d,MPI_DOUBLE,prev,tag,MPI_COMM_WORLD,&requests[2]);

  //MPI_Waitall(6,requests,status);


  clock_gettime(CLOCK_MONOTONIC, &ts_start);

  knn = search_k_nearest(vpt_tree.coords, vpt_tree.indexes, vpt_tree.mu, Y, n, n, d, k);


  //adapt indexes to point the right element.
  for(int i=0;i<n*k;i++)knn.nidx[i] = knn.nidx[i] + p_rank*n;

  //wait all send and receives in order to compute k nearest to the next vantage point tree
  MPI_Waitall(6,requests,status);

  for(int i=0;i<p_size-1;i++){
    //cope receive buffers in order to use them for the next vantage point tree.
    memcpy(vpt_tree.indexes,receive_buff_indexes,n*sizeof(int));
    memcpy(vpt_tree.mu     ,receive_buff_mu     ,n*sizeof(double));
    memcpy(vpt_tree.coords,receive_buff_coords,n*d*sizeof(double));

    //send the current vantage point to the next one
	MPI_Isend(vpt_tree.indexes, n,MPI_INT   ,next,tag,MPI_COMM_WORLD, &requests[3]);
	MPI_Isend(vpt_tree.mu     , n,MPI_DOUBLE,next,tag,MPI_COMM_WORLD, &requests[4]);
	MPI_Isend(vpt_tree.coords,n*d,MPI_DOUBLE,next,tag,MPI_COMM_WORLD, &requests[5]);

    MPI_Irecv(receive_buff_indexes, n,MPI_INT   ,prev,tag,MPI_COMM_WORLD,&requests[0]);
	MPI_Irecv(receive_buff_mu     , n,MPI_DOUBLE,prev,tag,MPI_COMM_WORLD,&requests[1]);
	MPI_Irecv(receive_buff_coords,n*d,MPI_DOUBLE,prev,tag,MPI_COMM_WORLD,&requests[2]);

    knn_temp = search_k_nearest(vpt_tree.coords,vpt_tree.indexes,vpt_tree.mu,Y,n,n,d,k);
    knn_combine(knn,knn_temp,n,k,(p_rank-i-1+p_size)%p_size);

    MPI_Waitall(6,requests,status);
  }

  clock_gettime(CLOCK_MONOTONIC, &ts_end);
  time = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
  //printf("search tree time: %lf p_id=%d\n",time,p_rank);


  free(Y);
  free(receive_buff_indexes);
  free(receive_buff_mu);
  free(receive_buff_coords);
  free(knn_temp.nidx);
  free(knn_temp.ndist);
  return knn;
}
