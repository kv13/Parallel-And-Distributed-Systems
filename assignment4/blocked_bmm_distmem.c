/* This implementation is based on blocked matrix multiplication             */
/* It provide two interfaces. One for boolean matrix multiplication          */
/* and one for filtered boolean matrix multiplication. Also to speed         */
/* up the processes openmpi for distributed memory parallelization is used   */
/* Also it can be used shared memory parallelization for further optimization*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

#include "inc/mmio.h"
#include "inc/utils.h"
#include "inc/blocking.h"
#include "inc/bmm.h"
#include "inc/parallel_bmm.h"

void interface_1(int argc, char *argv[]);
void interface_2(int argc, char *argv[]);
void interface_3(int argc, char *argv[]);

int main(int argc, char *argv[]){

    // make sure that are exactly the correct number of arguments.
    if(argc != 5 && argc!=4){
        fprintf(stderr,"BMM: A[martix-market-filename] B[martix-market-filename] [value of b] \n");
        fprintf(stderr,"FILTERED BMM: F[martix-market-filename] A[martix-market-filename] B[martix-market-filename] [value of b] \n");
        exit(1);
    }

    if(argc == 4){
        // BOOLEAN MATRIX MULTIPLICATION
        interface_1(argc, argv);
    }
    if(argc == 5){
        // FILTERED BOOLEAN MATRIX MULTIPLICATION
        interface_2(argc,argv);
    }

}

void interface_1(int argc, char *argv[]){

    // declare some variables
    int b, nz, nb, N, M, *coo_I, *coo_J;
    market_matrix *mtx_B;

    // declare time variables
    struct timespec ts_start;
    struct timespec ts_end;
    double time = 0, time_1;

    // INITIALIZE MPI COMMUNICATIONS
    int p_size, p_rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
    
    MPI_Status status; // status communication
    int tag = 1;       // communication tag

    // matrix B is red by all nodes
    mtx_B = (market_matrix *)malloc(sizeof(market_matrix)); 
    if(p_rank == 0){
        read_mm_matrices(mtx_B, argv[2]);
    }

    // broadcast matrix B information
    
    MPI_Bcast(&mtx_B->N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mtx_B->M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mtx_B->nz,1, MPI_INT, 0, MPI_COMM_WORLD);

    if(p_rank !=0){
        // allocate memory for matrix B in rest of nodes
        mtx_B->coo_I = (int *)malloc(sizeof(int)*mtx_B->nz);
        mtx_B->coo_J = (int *)malloc(sizeof(int)*mtx_B->nz); 
    }

    MPI_Bcast(mtx_B->coo_I, mtx_B->nz, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(mtx_B->coo_J, mtx_B->nz, MPI_INT, 0, MPI_COMM_WORLD);
    
    // master reads and distribute matrix A
    if(p_rank == 0){
        
        // read variable b
        b  = atoi(argv[3]);

        market_matrix *mtx_A;
        mtx_A = (market_matrix *)malloc(sizeof(market_matrix));
        read_mm_matrices(mtx_A, argv[1]);

        if(mtx_A->N % p_size != 0){
            printf("Dimensions Problem \n");
            exit(1);
        }
        
        // compute rows per node and total nonzero element per node
        N            = mtx_A->N / p_size;
        M            = mtx_A->M;
        int *temp_nz = (int *)calloc(p_size,sizeof(int));

        for(int i=0;i<mtx_A->nz;i++)
            temp_nz[mtx_A->coo_I[i] / N]++;

        // send matrix A information
        for(int i=1;i<p_size;i++){
            MPI_Send(&N          , 1, MPI_INT, i, tag, MPI_COMM_WORLD);
            MPI_Send(&mtx_A->M   , 1, MPI_INT, i, tag, MPI_COMM_WORLD);
            MPI_Send(&b          , 1, MPI_INT, i, tag, MPI_COMM_WORLD);
            MPI_Send(&temp_nz[i] , 1, MPI_INT, i, tag, MPI_COMM_WORLD);
        }
        
        // master's matrices
        nz = temp_nz[0];
        int master_idx = 0;
        coo_I = (int *)malloc(sizeof(int) * nz);
        coo_J = (int *)malloc(sizeof(int) * nz);

        // slaves matrices 
        int *temp_idx = (int *) calloc(p_size - 1, sizeof(int));
        int **temp_I  = (int **)malloc((p_size - 1)*sizeof(int *));
        int **temp_J  = (int **)malloc((p_size - 1)*sizeof(int *));

        for(int i=0;i<p_size-1;i++){
            temp_I[i] = (int *)malloc(sizeof(int) * temp_nz[i+1]);
            temp_J[i] = (int *)malloc(sizeof(int) * temp_nz[i+1]);
        }

        // distibute points to correct node
        int t_n;
        for(int i=0;i<mtx_A->nz;i++){
            t_n = mtx_A->coo_I[i] / N;
            if(t_n == 0){
                coo_I[master_idx] = mtx_A->coo_I[i];
                coo_J[master_idx] = mtx_A->coo_J[i];
                master_idx ++;
            }
            else{
                // fix indexes to distribute correct the matrices
                temp_I[t_n-1][temp_idx[t_n-1]] = mtx_A->coo_I[i] - t_n*N;
                temp_J[t_n-1][temp_idx[t_n-1]] = mtx_A->coo_J[i];
                temp_idx[t_n-1]++;
            }
        }

        // actually sending process
        for(int i=0; i<p_size-1; i++){
            MPI_Send(temp_I[i],temp_nz[i+1], MPI_INT, i+1, tag, MPI_COMM_WORLD);
            MPI_Send(temp_J[i],temp_nz[i+1], MPI_INT, i+1, tag, MPI_COMM_WORLD);
        }

        // free useless memory space
        for(int i=0; i<p_size-1;i++){
            free(temp_I[i]);
            free(temp_J[i]);
        }
        free(temp_I);
        free(temp_J);
        free(temp_nz);
        free(temp_idx);

        free(mtx_A->coo_I);
        free(mtx_A->coo_J);
        free(mtx_A);
    }
    else{
        
        // receive matrix A information
        MPI_Recv(&N  , 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&M  , 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&b  , 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&nz , 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);

        coo_I = (int *)malloc(sizeof(int) * nz);
        coo_J = (int *)malloc(sizeof(int) * nz);

        MPI_Recv(coo_I, nz, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(coo_J, nz, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
    }

    // transform submatrix A to csr format
    csr_format *csr_A;
    csr_A = (csr_format *)malloc(sizeof(csr_format));
    csr_init2(csr_A,N,M,nz);
    
    coo2csc(csr_A->csr_col, csr_A->csr_row, coo_J, coo_I, csr_A->nz, csr_A->N, 0);
    
    // free useless memory
    free(coo_I);
    free(coo_J);

    nb = csr_A->M / b;

    // ************** Blocking **************
    csr_format  *csr_AA;
    csr_AA = (csr_format *)malloc(sizeof(csr_format));

    blocked_csr *blk_csr_A;
    blk_csr_A = (blocked_csr *)malloc(sizeof(blocked_csr));

    blocked_csr_init(blk_csr_A, csr_A, b);

    clock_gettime(CLOCK_MONOTONIC, &ts_start);
    //****************************************/

    blocking_csr3(blk_csr_A, csr_A);
    blocked_2_csr(blk_csr_A, csr_AA);

    //****************************************/
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    time_1 = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
    printf("node %d: blocking time for matrix A %lf seconds \n",p_rank, time_1/(double)(1000000));
        
    //free useless memory
    csr_free(csr_A);

    // transform matrix B to csc format
    csc_format *csc_B;
    csc_B = (csc_format *)malloc(sizeof(csc_format));
    csc_init(csc_B, mtx_B);

    coo2csc(csc_B->csc_row, csc_B->csc_col, mtx_B->coo_I, mtx_B->coo_J, mtx_B->nz, mtx_B->M, 0);

    // free useless memory space
    free(mtx_B->coo_I);
    free(mtx_B->coo_J);
    free(mtx_B);

    // ************** Blocking **************
    csc_format *csc_BB;
    csc_BB = (csc_format *)malloc(sizeof(csc_format));

    blocked_csc *blk_csc_B;
    blk_csc_B = (blocked_csc *)malloc(sizeof(blocked_csc));
    blocked_csc_init(blk_csc_B, csc_B, b);

    clock_gettime(CLOCK_MONOTONIC, &ts_start);
    //****************************************/

    blocking_csc3(blk_csc_B, csc_B);
    blocked_2_csc(blk_csc_B, csc_BB);

    //****************************************/
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    time_1 = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
    printf("node %d: blocking time for matrix B %lf seconds \n",p_rank, time_1/(double)(1000000));

    //free useless memory
    csc_free(csc_B);

    // ****************** BMM ******************
    coo_format *coo_C;
    coo_C = (coo_format *)malloc(sizeof(coo_format));
    coo_init(coo_C,blk_csr_A,blk_csc_B);

    clock_gettime(CLOCK_MONOTONIC, &ts_start);
    // choose to use additionally and shared memory parallalization
    //****************************************/

    blocked_bmm(coo_C,blk_csr_A,csr_AA,blk_csc_B,csc_BB);
    //blocked_bmm_shared(coo_C,blk_csr_A,csr_AA,blk_csc_B,csc_BB);
    
    //****************************************/
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    time_1 = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
    printf("node %d: sub matrix multiplication %lf seconds \n",p_rank, time_1/(double)(1000000));

    if(p_rank == 0)
        time += time_1;

    // free matrices.
    csc_free(csc_BB);
    blocked_csc_free(blk_csc_B);
    
    csr_free(csr_AA);
    blocked_csr_free(blk_csr_A);
    
    // gather all results
    if(p_rank == 0){

        int total_size = coo_C->cur_nz;
        int sub_size[p_size-1];

        // compute total number of non zero element
        for(int i=1;i<p_size;i++){
            MPI_Recv(&sub_size[i-1], 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
        }

        for(int i=0;i<p_size-1;i++)
            total_size += sub_size[i];
        
        // define C matrix
        coo_format *total_C = (coo_format *)malloc(sizeof(coo_format));

        total_C->max_nz = total_size;
        total_C->cur_nz = 0;
        total_C->coo_I = (int *)malloc(sizeof(int)*total_C->max_nz);
        total_C->coo_J = (int *)malloc(sizeof(int)*total_C->max_nz);

        // append all points from master first
        for(int k=0;k<coo_C->cur_nz;k++){
            total_C->coo_I[total_C->cur_nz] = coo_C->coo_I[k]; 
            total_C->coo_J[total_C->cur_nz] = coo_C->coo_J[k];
            total_C->cur_nz++;
        }

        coo_free(coo_C);

        // receive all points from slaves
        int *temp_I, *temp_J;
        for(int i=0;i<p_size-1;i++){
            temp_I = (int *)malloc(sizeof(int)*sub_size[i]);
            temp_J = (int *)malloc(sizeof(int)*sub_size[i]);
            MPI_Recv(temp_I, sub_size[i], MPI_INT, i+1, tag, MPI_COMM_WORLD, &status);
            MPI_Recv(temp_J, sub_size[i], MPI_INT, i+1, tag, MPI_COMM_WORLD, &status);

            for(int k=0;k<sub_size[i];k++){
                total_C->coo_I[total_C->cur_nz] = temp_I[k]; 
                total_C->coo_J[total_C->cur_nz] = temp_J[k];
                total_C->cur_nz++;
            }
            free(temp_I);
            free(temp_J);
        }
    }
    else{

        // send number nz elements of submatrix C
        MPI_Send(&coo_C->cur_nz, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);

        // fix index 
        int fix = p_rank * N;
        for(int k=0;k<coo_C->cur_nz;k++)
            coo_C->coo_I[k] = coo_C->coo_I[k] + fix;
        
        // send coordinates
        MPI_Send(coo_C->coo_I, coo_C->cur_nz, MPI_INT, 0, tag, MPI_COMM_WORLD);
        MPI_Send(coo_C->coo_J, coo_C->cur_nz, MPI_INT, 0, tag, MPI_COMM_WORLD);

        //free useless memory
        coo_free(coo_C);

    }

    //execution is over
    MPI_Finalize();

    return ;
}


void interface_2(int argc, char *argv[]){
    
    // declare variables for matrix A
    int b, nz, N, M, *coo_I, *coo_J;
    // declare variables for matrix F
    int F_nz, F_N, F_M, *F_coo_I, *F_coo_J;
    // declare variables for matrix B;
    market_matrix *mtx_B;

    // declare time variables
    struct timespec ts_start;
    struct timespec ts_end;
    double time = 0, time_1;

    // INITIALIZE MPI COMMUNICATIONS
    int p_size, p_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);

    MPI_Status status; // status communication
    int tag = 1;       // communication tag

    // matrix B is red by all nodes
    mtx_B = (market_matrix *)malloc(sizeof(market_matrix)); 

    if(p_rank == 0){
        // read variable b
        b  = atoi(argv[4]);
        read_mm_matrices(mtx_B, argv[3]);
    }

    // broadcast matrix B information
    
    MPI_Bcast(&b       , 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mtx_B->N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mtx_B->M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mtx_B->nz,1, MPI_INT, 0, MPI_COMM_WORLD);

    if(p_rank !=0){
        // allocate memory for matrix B in rest of nodes
        mtx_B->coo_I = (int *)malloc(sizeof(int)*mtx_B->nz);
        mtx_B->coo_J = (int *)malloc(sizeof(int)*mtx_B->nz); 
    }

    MPI_Bcast(mtx_B->coo_I, mtx_B->nz, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(mtx_B->coo_J, mtx_B->nz, MPI_INT, 0, MPI_COMM_WORLD);

    // transform matrix B to csc format
    csc_format *csc_B;
    csc_B = (csc_format *)malloc(sizeof(csc_format));
    csc_init(csc_B, mtx_B);

    coo2csc(csc_B->csc_row, csc_B->csc_col, mtx_B->coo_I, mtx_B->coo_J, mtx_B->nz, mtx_B->M, 0);

    // free useless memory space
    free(mtx_B->coo_I);
    free(mtx_B->coo_J);
    free(mtx_B);

    // ************** Blocking **************
    csc_format *csc_BB;
    csc_BB = (csc_format *)malloc(sizeof(csc_format));

    blocked_csc *blk_csc_B;
    blk_csc_B = (blocked_csc *)malloc(sizeof(blocked_csc));
    blocked_csc_init(blk_csc_B, csc_B, b);

    clock_gettime(CLOCK_MONOTONIC, &ts_start);
    //****************************************/

    blocking_csc3(blk_csc_B, csc_B);
    blocked_2_csc(blk_csc_B, csc_BB);

    //****************************************/
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    time_1 = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
    printf("node %d: blocking time for matrix B %lf seconds \n",p_rank, time_1/(double)(1000000));

    if(p_rank == 0)
        time += time_1;

    //free useless memory
    csc_free(csc_B);

    if(p_rank == 0){

        market_matrix *mtx_A;
        mtx_A = (market_matrix *)malloc(sizeof(market_matrix));
        read_mm_matrices(mtx_A, argv[2]);

        if(mtx_A->N % p_size != 0){
            printf("Dimensions Problem \n");
            exit(1);
        }

        
        // compute rows per node and total nonzero element per node
        N            = mtx_A->N / p_size;
        M            = mtx_A->M;
        int *temp_nz = (int *)calloc(p_size,sizeof(int));

        for(int i=0;i<mtx_A->nz;i++)
            temp_nz[mtx_A->coo_I[i] / N]++;

        // send matrix A information
        for(int i=1;i<p_size;i++){
            MPI_Send(&N          , 1, MPI_INT, i, tag, MPI_COMM_WORLD);
            MPI_Send(&mtx_A->M   , 1, MPI_INT, i, tag, MPI_COMM_WORLD);
            MPI_Send(&temp_nz[i] , 1, MPI_INT, i, tag, MPI_COMM_WORLD);
        }
        
        // master's matrices
        nz = temp_nz[0];
        int master_idx = 0;
        coo_I = (int *)malloc(sizeof(int) * nz);
        coo_J = (int *)malloc(sizeof(int) * nz);

        // slaves matrices 
        int *temp_idx = (int *) calloc(p_size - 1, sizeof(int));
        int **temp_I  = (int **)malloc((p_size - 1)*sizeof(int *));
        int **temp_J  = (int **)malloc((p_size - 1)*sizeof(int *));

        for(int i=0;i<p_size-1;i++){
            temp_I[i] = (int *)malloc(sizeof(int) * temp_nz[i+1]);
            temp_J[i] = (int *)malloc(sizeof(int) * temp_nz[i+1]);
        }

        // distibute points to correct node
        int t_n;
        for(int i=0;i<mtx_A->nz;i++){
            t_n = mtx_A->coo_I[i] / N;
            if(t_n == 0){
                coo_I[master_idx] = mtx_A->coo_I[i];
                coo_J[master_idx] = mtx_A->coo_J[i];
                master_idx ++;
            }
            else{
                // fix indexes to distribute correct the matrices
                temp_I[t_n-1][temp_idx[t_n-1]] = mtx_A->coo_I[i] - t_n*N;
                temp_J[t_n-1][temp_idx[t_n-1]] = mtx_A->coo_J[i];
                temp_idx[t_n-1]++;
            }
        }

        // actually sending process
        for(int i=0; i<p_size-1; i++){
            MPI_Send(temp_I[i],temp_nz[i+1], MPI_INT, i+1, tag, MPI_COMM_WORLD);
            MPI_Send(temp_J[i],temp_nz[i+1], MPI_INT, i+1, tag, MPI_COMM_WORLD);
        }

        // free useless memory space
        for(int i=0; i<p_size-1;i++){
            free(temp_I[i]);
            free(temp_J[i]);
        }
        free(temp_I);
        free(temp_J);
        free(temp_nz);
        free(temp_idx);

        free(mtx_A->coo_I);
        free(mtx_A->coo_J);
        free(mtx_A);

        // read and distribute filter
        market_matrix *mtx_F;
        mtx_F = (market_matrix *)malloc(sizeof(market_matrix));
        read_mm_matrices(mtx_F, argv[1]);
        
        // compute rows per node and total nonzero element per node
        F_N     = mtx_F->N / p_size;
        F_M     = mtx_F->M;
        temp_nz = (int *)calloc(p_size,sizeof(int));

        for(int i=0;i<mtx_F->nz;i++)
            temp_nz[mtx_F->coo_I[i] / F_N]++;

        // send filter F information
        for(int i=1;i<p_size;i++){
            MPI_Send(&F_N       , 1, MPI_INT, i, tag, MPI_COMM_WORLD);
            MPI_Send(&F_M       , 1, MPI_INT, i, tag, MPI_COMM_WORLD);
            MPI_Send(&temp_nz[i], 1, MPI_INT, i, tag, MPI_COMM_WORLD);
        }

        // master's F matrices
        F_nz = temp_nz[0];
        master_idx = 0;
        F_coo_I = (int *)malloc(sizeof(int) * F_nz);
        F_coo_J = (int *)malloc(sizeof(int) * F_nz);

        // slaves matrices 
        temp_idx = (int *) calloc(p_size - 1, sizeof(int));
        temp_I   = (int **)malloc((p_size - 1)*sizeof(int *));
        temp_J   = (int **)malloc((p_size - 1)*sizeof(int *));

        for(int i=0;i<p_size-1;i++){
            temp_I[i] = (int *)malloc(sizeof(int) * temp_nz[i+1]);
            temp_J[i] = (int *)malloc(sizeof(int) * temp_nz[i+1]);
        }

        // distibute points to correct node
        for(int i=0;i<mtx_F->nz;i++){
            t_n = mtx_F->coo_I[i] / F_N;
            if(t_n == 0){
                F_coo_I[master_idx] = mtx_F->coo_I[i];
                F_coo_J[master_idx] = mtx_F->coo_J[i];
                master_idx ++;
            }
            else{
                // fix indexes to distribute correct the matrices
                temp_I[t_n-1][temp_idx[t_n-1]] = mtx_F->coo_I[i] - t_n*F_N;
                temp_J[t_n-1][temp_idx[t_n-1]] = mtx_F->coo_J[i];
                temp_idx[t_n-1]++;
            }
        }

        // actually sending process
        for(int i=0; i<p_size-1; i++){
            MPI_Send(temp_I[i],temp_nz[i+1], MPI_INT, i+1, tag, MPI_COMM_WORLD);
            MPI_Send(temp_J[i],temp_nz[i+1], MPI_INT, i+1, tag, MPI_COMM_WORLD);
        }

        // free useless memory space
        for(int i=0; i<p_size-1;i++){
            free(temp_I[i]);
            free(temp_J[i]);
        }
        free(temp_I);
        free(temp_J);
        free(temp_nz);
        free(temp_idx);

        free(mtx_F->coo_I);
        free(mtx_F->coo_J);
        free(mtx_F);
    }
    else{
        // receive matrix A information
        MPI_Recv(&N  , 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&M  , 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&nz , 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);

        coo_I = (int *)malloc(sizeof(int) * nz);
        coo_J = (int *)malloc(sizeof(int) * nz);

        MPI_Recv(coo_I, nz, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(coo_J, nz, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);

        // receive filter F information
        MPI_Recv(&F_N  , 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&F_M  , 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&F_nz , 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);

        F_coo_I = (int *)malloc(sizeof(int) * F_nz);
        F_coo_J = (int *)malloc(sizeof(int) * F_nz);

        MPI_Recv(F_coo_I, F_nz, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(F_coo_J, F_nz, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);

    }

    // transform submatrix A to csr format
    csr_format *csr_A;
    csr_A = (csr_format *)malloc(sizeof(csr_format));
    csr_init2(csr_A,N,M,nz);
    
    coo2csc(csr_A->csr_col, csr_A->csr_row, coo_J, coo_I, csr_A->nz, csr_A->N, 0);
    
    // free useless memory
    free(coo_I);
    free(coo_J);

    // ************** Blocking **************
    csr_format  *csr_AA;
    csr_AA = (csr_format *)malloc(sizeof(csr_format));

    blocked_csr *blk_csr_A;
    blk_csr_A = (blocked_csr *)malloc(sizeof(blocked_csr));

    blocked_csr_init(blk_csr_A, csr_A, b);

    clock_gettime(CLOCK_MONOTONIC, &ts_start);
    //****************************************/

    blocking_csr3(blk_csr_A, csr_A);
    blocked_2_csr(blk_csr_A, csr_AA);

    //****************************************/
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    time_1 = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
    printf("node %d: blocking time for matrix A %lf seconds \n",p_rank, time_1/(double)(1000000));

    if(p_rank == 0)
        time += time_1;
        
    //free useless memory
    csr_free(csr_A);

    // transform filter F to csr format
    csr_format *csr_F;
    csr_F = (csr_format *)malloc(sizeof(csr_format));
    csr_init2(csr_F,F_N,F_M,F_nz);

    coo2csc(csr_F->csr_col, csr_F->csr_row, F_coo_J, F_coo_I, csr_F->nz, csr_F->N, 0);

    // free useless memory
    free(F_coo_I);
    free(F_coo_J);

    // ************** Blocking **************
    csr_format  *csr_FF;
    csr_FF = (csr_format *)malloc(sizeof(csr_format));
    
    blocked_csr *blk_csr_F;
    blk_csr_F = (blocked_csr *)malloc(sizeof(blocked_csr));
    
    blocked_csr_init(blk_csr_F, csr_F, b);

    clock_gettime(CLOCK_MONOTONIC, &ts_start);
    //****************************************/

    blocking_csr3(blk_csr_F, csr_F);
    blocked_2_csr(blk_csr_F, csr_FF);

    //****************************************/
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    time_1 = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
    printf("node %d: blocking time for matrix F %lf seconds \n",p_rank, time_1/(double)(1000000));
    
    if(p_rank == 0)
        time += time_1;
    
    //free useless memory
    csr_free(csr_F);

    // ************** BMM FILTERED **************
    coo_format *coo_C;
    coo_C = (coo_format *)malloc(sizeof(coo_format));
    coo_init_filter(coo_C,blk_csr_F);

    clock_gettime(CLOCK_MONOTONIC, &ts_start);
    //****************************************/
    
    blocked_bmm_filtered(coo_C,blk_csr_F,csr_FF,blk_csr_A,csr_AA,blk_csc_B,csc_BB);
    //blocked_bmm_filtered_shared(coo_C,blk_csr_F,csr_FF,blk_csr_A,csr_AA,blk_csc_B,csc_BB);

    //****************************************/
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    time_1 = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
    printf("node %d: submatrix filtered multiplication %lf seconds \n",p_rank, time_1/(double)(1000000));

    if(p_rank == 0)
        time += time_1;
    
    // free matrices.
    csr_free(csr_FF);
    blocked_csr_free(blk_csr_F);

    csr_free(csr_AA);
    blocked_csr_free(blk_csr_A);

    csc_free(csc_BB);
    blocked_csc_free(blk_csc_B);
    
    // gather all results
    if(p_rank == 0){

        int total_size = coo_C->cur_nz;
        int sub_size[p_size-1];

        // compute total number of non zero element
        for(int i=1;i<p_size;i++){
            MPI_Recv(&sub_size[i-1], 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
        }

        for(int i=0;i<p_size-1;i++)
            total_size += sub_size[i];
        
        // define C matrix
        coo_format *total_C = (coo_format *)malloc(sizeof(coo_format));

        total_C->max_nz = total_size;
        total_C->cur_nz = 0;
        total_C->coo_I = (int *)malloc(sizeof(int)*total_C->max_nz);
        total_C->coo_J = (int *)malloc(sizeof(int)*total_C->max_nz);

        // append all points from master first
        for(int k=0;k<coo_C->cur_nz;k++){
            total_C->coo_I[total_C->cur_nz] = coo_C->coo_I[k]; 
            total_C->coo_J[total_C->cur_nz] = coo_C->coo_J[k];
            total_C->cur_nz++;
        }

        coo_free(coo_C);

        // receive all points from slaves
        int *temp_I, *temp_J;
        for(int i=0;i<p_size-1;i++){
            if(sub_size[i]!=0){
                temp_I = (int *)malloc(sizeof(int)*sub_size[i]);
                temp_J = (int *)malloc(sizeof(int)*sub_size[i]);
                MPI_Recv(temp_I, sub_size[i], MPI_INT, i+1, tag, MPI_COMM_WORLD, &status);
                MPI_Recv(temp_J, sub_size[i], MPI_INT, i+1, tag, MPI_COMM_WORLD, &status);
                
                for(int k=0;k<sub_size[i];k++){
                    total_C->coo_I[total_C->cur_nz] = temp_I[k]; 
                    total_C->coo_J[total_C->cur_nz] = temp_J[k];
                    total_C->cur_nz++;
                }
                free(temp_I);
                free(temp_J);
            }
        }
    }
    else{
        
        // send number nz elements of submatrix C
        MPI_Send(&coo_C->cur_nz, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);

        // fix index 
        int fix = p_rank * N;
        for(int k=0;k<coo_C->cur_nz;k++)
            coo_C->coo_I[k] = coo_C->coo_I[k] + fix;

        // send coordinates
        if(coo_C->cur_nz != 0){
            MPI_Send(coo_C->coo_I, coo_C->cur_nz, MPI_INT, 0, tag, MPI_COMM_WORLD);
            MPI_Send(coo_C->coo_J, coo_C->cur_nz, MPI_INT, 0, tag, MPI_COMM_WORLD);
        }

        //free useless memory
        coo_free(coo_C);
    }

    //execution is over
    MPI_Finalize();
    return ;
}