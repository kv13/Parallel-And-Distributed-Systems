#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <pthread.h>
#include <time.h>

#include "../lib/mmio.h"
#include "../lib/triangles_library.h"

//global mutex to used by threads.
pthread_mutex_t mut;


//V3 ALGORITHM USING CIKL
void V3_algorithm_cilk(uint32_t const * const csc_col,uint32_t const * const csc_row, uint32_t * const c3, int n);


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
      I = (uint32_t *) malloc(2*nz * sizeof(uint32_t));
      J = (uint32_t *) malloc(2*nz * sizeof(uint32_t));

      /*read the market matrix file*/
      for (uint32_t i=0; i<2*nz; i=i+2)
      {
          if(fscanf(f, "%u %u \n", &I[i], &J[i]))
          ;
          /* adjust from 1-based to 0-based */
          I[i]--;
          J[i]--;
          //create the symmetric matrix.
          I[i+1]=J[i];
          J[i+1]=I[i];
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

      printf("====================PARALLEL V3 ALGORITHM USING CILK==================== \n");
      clock_gettime(CLOCK_MONOTONIC, &ts_start);
      if(N<100){
        printf("FEW NODES THUS SEQUENTIAL ALGORITHM IS SELECTED \n");
        V3_algorithm(csc_col,csc_row,c3,N);
      }
      else V3_algorithm_cilk(csc_col,csc_row,c3,N);
      clock_gettime(CLOCK_MONOTONIC, &ts_end);

      //#######################WRITE RESULTS TO FILE AND EXIT#######################
      char str[200];
      snprintf(str,sizeof(str),"V3_CILK.txt");
      double time = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
      write_to_file(str,c3,N,time);
      printf("RESULTS HAVE BEEN WRITTEN\n");
      printf("EXITING...\n");
      free(c3);
      free(csc_col);
      free(csc_row);
      free(I);
      free(J);
      return 0;
}


//V3 ALGORITHM USING CIKL
void V3_algorithm_cilk(uint32_t const * const csc_col,uint32_t const * const csc_row, uint32_t * const c3, int n){

  //check the number of workers.
  //int numWorkers = __cilkrts_get_nworkers();
  //printf("The number of workers is %d\n",numWorkers);

  //parallel implementation of V3 algorithm.
  cilk_for (uint32_t i = 0;i<n;i++){
    uint32_t j;
    uint32_t k;
    for(uint32_t temp1 = 0;temp1<csc_col[i+1]-csc_col[i];temp1++){
      j = csc_row[csc_col[i]+temp1];
      if(j<i+1)continue;
      for(uint32_t temp2=0;temp2<csc_col[j+1]-csc_col[j];temp2++){
        k = csc_row[csc_col[j]+temp2];
        if(k<j+1)continue;
        if(binary_search(k,i,temp1,csc_col,csc_row) !=0 ){

          //we must use lock to avoid data races.
          pthread_mutex_lock(&mut);
          c3[i]++;
          c3[j]++;
          c3[k]++;
          pthread_mutex_unlock(&mut);
        }
      }
    }
  }
}
