#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include <math.h>
#include <time.h>

#include "../lib/mmio.h"
#include "../lib/triangles_library.h"

#define NUM_THREADS 15
pthread_mutex_t mut;


typedef struct{
	uint32_t t_id;
	uint32_t l_bound, u_bound;
	uint32_t workload;
	uint32_t *csc_col;
	uint32_t *csc_row;
	uint32_t *c3;
	uint32_t *val;
	pthread_t *id;
}tag;

//functions for threads.
void Hadamard_algorithm_pthread(uint32_t *csc_col, uint32_t *csc_row, uint32_t * const c3, uint32_t * const val, int n);
void *count_triangles_parallel(void *arg);


int main(int argc, char *argv[])
{
  //#################Read the sparse matrix from file.#################
      int ret_code;
      MM_typecode matcode;
      FILE *f;
      int M, N, nz;


      if (argc < 2){
        fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
        exit(1);
      }else{
          if ((f = fopen(argv[1], "r")) == NULL){
            printf("ERROR: Cannot process the file\n" );
            exit(1);
          }
      }

      if (mm_read_banner(f, &matcode) != 0)
      {
          printf("Could not process Matrix Market banner.\n");
          exit(1);
      }


      if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode) )
      {
          printf("Sorry, this application does not support ");
          printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
          exit(1);
      }

      /* find out size of sparse matrix .... */
      if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0){
        printf("Sorry, this application does not support ");
        exit(1);
      }

      /* reseve memory for matrices */
      uint32_t *I, *J;
      uint32_t *val;
      I = (uint32_t *) malloc(2*nz * sizeof(uint32_t));
      J = (uint32_t *) malloc(2*nz * sizeof(uint32_t));
      val = (uint32_t *) malloc(2*nz * sizeof(uint32_t));

      /*read the market matrix file*/
      for (uint32_t i=0; i<2*nz; i=i+2)
      {
          fscanf(f, "%u %u \n", &I[i], &J[i]);
          /* adjust from 1-based to 0-based */
          I[i]--;
          J[i]--;
          //create the symmetric matrix.
          I[i+1]=J[i];
          J[i+1]=I[i];
          val[i]=0;
      }
      if (f !=stdin) fclose(f);


      //#################CSC FORMAT FOR THE SPARSE MATRIX#################
      uint32_t *csc_col= (uint32_t *)malloc((N + 1) * sizeof(uint32_t));
      uint32_t *csc_row= (uint32_t *)malloc(2*nz * sizeof(uint32_t));
      uint32_t *c3 = (uint32_t *)malloc(N*sizeof(uint32_t));
      coo2csc(csc_row, csc_col, I, J, (uint32_t)2*nz, (uint32_t)N,0);

      //initialize c3's matrices
      for(int i = 0;i<N;i++){
        c3[i]=0;
      }

      /*timespec variables to count the total time of execution */
      struct timespec ts_start;
      struct timespec ts_end;

      printf("====================PTHREAD HADAMARD ALGORITHM==================== \n");
      clock_gettime(CLOCK_MONOTONIC, &ts_start);
      if(N<100){
        printf("FEW NODES THUS SEQUENTIAL ALGORITHM IS SELECTED\n");
        Hadamard_algorithm(csc_col,csc_row,c3,val,N);
      }
      else Hadamard_algorithm_pthread(csc_col,csc_row,c3,val,N);
      clock_gettime(CLOCK_MONOTONIC, &ts_end);

      //#######################WRITE RESULTS TO FILE AND EXIT#######################
      char str[200];
      snprintf(str,sizeof(str),"HADAMARD_PTHREADS.txt");
      double time = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
      write_to_file(str,c3,N,time);
      printf("RESULTS HAVE BEEN WRITTEN\n");
      printf("EXITING...\n");
      free(c3);
      free(val);
      free(csc_col);
      free(csc_row);
      free(I);
      free(J);
      return 0;
}

void Hadamard_algorithm_pthread(uint32_t *csc_col, uint32_t *csc_row, uint32_t * const c3, uint32_t * const val, int n){

  //declare thread variables.
  pthread_t *threads_id;
  threads_id = (pthread_t *)malloc(NUM_THREADS*sizeof(pthread_t));
  uint32_t workload;
  workload = floor((float)n/NUM_THREADS);
  tag *tags = (tag *)malloc(NUM_THREADS*sizeof(tag));

  //initialize and start threads.
  for(uint32_t i=0;i<NUM_THREADS;i++){
    tags[i].t_id = i;
    tags[i].workload = workload;
    tags[i].l_bound = i*workload;
    tags[i].u_bound = (i+1)*workload;
    tags[i].csc_col = csc_col;
    tags[i].csc_row = csc_row;
    tags[i].c3 = c3;
    tags[i].val = val;
    tags[i].id = &threads_id[i];
    pthread_create(tags[i].id,NULL,count_triangles_parallel,(void *)&tags[i]);
  }
  //....................DO WORK FROM MAIN..................
  //COMPUTE TRIANGLES FOR THE REST OF THE NODES.
  if(n-workload*NUM_THREADS>0){
    tag main_tag;
    main_tag.t_id = NUM_THREADS;
    main_tag.workload = n-NUM_THREADS*workload;
    main_tag.l_bound = NUM_THREADS*workload;
    main_tag.u_bound = (uint32_t)n;
    main_tag.csc_row = csc_row;
    main_tag.csc_col = csc_col;
    main_tag.c3 = c3;
    main_tag.val = val;
    main_tag.id = NULL;
    count_triangles_parallel((void *)&main_tag);
  }
  for(int i=0;i<NUM_THREADS;i++){
    pthread_join(*tags[i].id,NULL);
  }
}



void *count_triangles_parallel(void *arg){
	tag *t = (tag *)arg;
	uint32_t i,j,k;
	uint32_t temp1,temp2;
	for(temp1=t->l_bound; temp1<t->u_bound; temp1++){
		for(temp2=0; temp2<t->csc_col[temp1+1]- t->csc_col[temp1]; temp2++){
			i = temp1;
			j = t->csc_row[t->csc_col[temp1]+temp2];
			t->val[t->csc_col[temp1]+temp2]=compute_sum(i, j, t->csc_col, t->csc_row);
			t->c3[temp1] = t->c3[temp1]+t->val[t->csc_col[temp1]+temp2];
		}
		t->c3[temp1]=t->c3[temp1]/2;
	}
}
