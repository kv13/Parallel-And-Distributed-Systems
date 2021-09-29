/* This implementation is based on blocked matrix multiplication    */
/* It provide two interfaces. One for boolean matrix multiplication */
/* and one for filtered boolean matrix multiplication               */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#include "inc/mmio.h"
#include "inc/utils.h"
#include "inc/blocking.h"
#include "inc/bmm.h"

void interface_1(int argc, char *argv[]);
void interface_2(int argc, char *argv[]);

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
    
    // read first matrix
    market_matrix *mtx_A;
    mtx_A = (market_matrix *)malloc(sizeof(market_matrix));
    read_mm_matrices(mtx_A, argv[1]);

    // transform it to csr format
    csr_format *csr_A;
    csr_A = (csr_format *)malloc(sizeof(csr_format));
    csr_init(csr_A,mtx_A);
    
    coo2csc(csr_A->csr_col, csr_A->csr_row, mtx_A->coo_J, mtx_A->coo_I, csr_A->nz, csr_A->N, 0);

    // free useless memory space
    free(mtx_A->coo_I);
    free(mtx_A->coo_J);
    free(mtx_A);

    int b  = atoi(argv[3]);
    int nb = csr_A->M / b;
    
    struct timespec ts_start;
    struct timespec ts_end;
    double time;
    
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
    time = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
    printf("blocking time for matrix A %lf seconds \n",time/(double)(1000000));

    //free useless memory
    csr_free(csr_A);

    // read second matrix
    market_matrix *mtx_B;
    mtx_B = (market_matrix *)malloc(sizeof(market_matrix));
    read_mm_matrices(mtx_B, argv[2]);

    // transform to csc format
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
    time = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
    printf("blocking time for matrix B %lf seconds \n",time/(double)(1000000));

    //free useless memory
    csc_free(csc_B);

    // ****************** BMM ******************
    coo_format *coo_C;
    coo_C = (coo_format *)malloc(sizeof(coo_format));
    coo_init(coo_C,blk_csr_A,blk_csc_B);

    clock_gettime(CLOCK_MONOTONIC, &ts_start);
    //****************************************/

    blocked_bmm(coo_C,blk_csr_A,csr_AA,blk_csc_B,csc_BB);

    //****************************************/
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    time = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
    printf("blocked multiplication took %lf seconds \n",time/(double)(1000000));
    
    // free matrices.
    csc_free(csc_BB);
    blocked_csc_free(blk_csc_B);
    
    csr_free(csr_AA);
    blocked_csr_free(blk_csr_A);

    // ************** test results **************
    // CAUTION: USE ONLY FOR SMALL MATRICES

    // read correct matrix C
    //char path_to_C[100] = "data/matrix_C_100k.mtx";
    //market_matrix *mtx_C_correct;
    //mtx_C_correct = (market_matrix *)malloc(sizeof(market_matrix));
    //read_mm_matrices(mtx_C_correct, path_to_C);

    // test results
    //quick_sort(coo_C->coo_I,coo_C->coo_J,0,coo_C->cur_nz-1);
    //test_results(mtx_C_correct, coo_C);
}

void interface_2(int argc, char *argv[]){
    
    // timining variables
    struct timespec ts_start;
    struct timespec ts_end;
    double time;
    
    // read first matrix
    market_matrix *mtx_A;
    mtx_A = (market_matrix *)malloc(sizeof(market_matrix));
    read_mm_matrices(mtx_A, argv[2]);

    // transform it to csr format
    csr_format *csr_A;
    csr_A = (csr_format *)malloc(sizeof(csr_format));
    csr_init(csr_A,mtx_A);
    
    coo2csc(csr_A->csr_col, csr_A->csr_row, mtx_A->coo_J, mtx_A->coo_I, csr_A->nz, csr_A->N, 0);

    // free useless memory space
    free(mtx_A->coo_I);
    free(mtx_A->coo_J);
    free(mtx_A);

    int b  = atoi(argv[4]);
    int nb = csr_A->M / b;

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
    time = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
    printf("blocking time for matrix A %lf seconds \n",time/(double)(1000000));

    //free useless memory
    csr_free(csr_A);

    // read second matrix
    market_matrix *mtx_B;
    mtx_B = (market_matrix *)malloc(sizeof(market_matrix));
    read_mm_matrices(mtx_B, argv[3]);

    // transform to csc format
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
    time = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
    printf("blocking time for matrix B %lf seconds \n",time/(double)(1000000));

    //free useless memory
    csc_free(csc_B);

    // read filter
    market_matrix *mtx_F;
    mtx_F = (market_matrix *)malloc(sizeof(market_matrix));
    read_mm_matrices(mtx_F, argv[1]);

    // transform it to csr format
    csr_format *csr_F;
    csr_F = (csr_format *)malloc(sizeof(csr_format));
    csr_init(csr_F,mtx_F);
    
    coo2csc(csr_F->csr_col, csr_F->csr_row, mtx_F->coo_J, mtx_F->coo_I, csr_F->nz, csr_F->N, 0);

    // free useless memory space
    free(mtx_F->coo_I);
    free(mtx_F->coo_J);
    free(mtx_F);

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
    time = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
    printf("blocking time for matrix F %lf seconds \n",time/(double)(1000000));

    //free useless memory
    csr_free(csr_F);

    // ************** BMM FILTERED **************
    coo_format *coo_C;
    coo_C = (coo_format *)malloc(sizeof(coo_format));
    coo_init_filter(coo_C,blk_csr_F);

    clock_gettime(CLOCK_MONOTONIC, &ts_start);
    //****************************************/

    blocked_bmm_filtered(coo_C,blk_csr_F,csr_FF,blk_csr_A,csr_AA,blk_csc_B,csc_BB);

    //****************************************/
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    time = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
    printf("blocked filtered multiplication took %lf seconds \n",time/(double)(1000000));
    
    // free matrices.
    csr_free(csr_FF);
    blocked_csr_free(blk_csr_F);

    csr_free(csr_AA);
    blocked_csr_free(blk_csr_A);

    csc_free(csc_BB);
    blocked_csc_free(blk_csc_B);

    // ************** test results **************
    // CAUTION: USE ONLY FOR SMALL MATRICES

    // read correct matrix C
    char path_to_C[150] = "data/matrix_C_f_3M.mtx";
    market_matrix *mtx_C_correct;
    mtx_C_correct = (market_matrix *)malloc(sizeof(market_matrix));
    read_mm_matrices(mtx_C_correct, path_to_C);

    // test results
    quick_sort(coo_C->coo_I,coo_C->coo_J,0,coo_C->cur_nz-1);
    test_results(mtx_C_correct, coo_C);

}