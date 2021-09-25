/* This implementation is based on the simplest algorithm */
/* for matrices multiplication with out any improvements. */
/* It's used only as benchmarck and result testing        */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#include "inc/mmio.h"
#include "inc/utils.h"
#include "inc/blocking.h"
#include "inc/bmm.h"

// ********** functions declaration ********** 
void   init_coo(coo_format *C, csr_format *A, csc_format *B);
void simple_bmm(coo_format *C, csr_format *A, csc_format *B);
void insert_new_element(coo_format *C, int row, int col);

int main(int argc, char *argv[]){

    // ********** Read Matrices **********
    // make sure that are exactly the correct number of arguments.
    if(argc != 4){
        fprintf(stderr,"Usage: [martix-market-filename] [martix-market-filename] [martix-market-filename]\n");
        exit(1);
    }
    
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

    // ************** BMM **************
    struct timespec ts_start;
    struct timespec ts_end;
    double time;

    coo_format *coo_C;
    coo_C = (coo_format *)malloc(sizeof(coo_format));
    init_coo(coo_C, csr_A, csc_B);

    clock_gettime(CLOCK_MONOTONIC, &ts_start);

    //****************************************/
    simple_bmm(coo_C, csr_A, csc_B);
    //****************************************/
    
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    time = 1000000*(double)(ts_end.tv_sec-ts_start.tv_sec)+(double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000;
    printf("boolean multiplication took %lf seconds \n",time/(double)(1000000));

    
    // ************** test results **************
    market_matrix *mtx_C_correct;
    mtx_C_correct = (market_matrix *)malloc(sizeof(market_matrix));
    // read correct matrix C
    read_mm_matrices(mtx_C_correct, argv[3]);
    
    test_results(mtx_C_correct, coo_C);
}


void init_coo(coo_format *C, csr_format *A, csc_format *B){
    
    C->cur_nz = 0;
    C->N      = A->N;
    C->M      = B->M;

    // suppose that the C matrix will have max(A,B) non zero elements
    if(A->nz>B->nz)
        C->max_nz = A->nz;
    else
        C->max_nz = B->nz;
    
    C->coo_I = (int *)malloc(sizeof(int)*C->max_nz);
    C->coo_J = (int *)malloc(sizeof(int)*C->max_nz);
}

void simple_bmm(coo_format *C, csr_format *A, csc_format *B){

    int ptr_A;
    int ptr_B;
    int new_flag;

    for(int row=0;row<C->N;row++){
        for(int col=0;col<C->M;col++){
            ptr_A     = A->csr_row[row];
            ptr_B     = B->csc_col[col];
            new_flag = 0;
            while(ptr_A < A->csr_row[row+1] && ptr_B < B->csc_col[col+1]){
                if(A->csr_col[ptr_A] > B->csc_row[ptr_B]){
                    ptr_B++;
                }
                else if(A->csr_col[ptr_A] < B->csc_row[ptr_B]){
                    ptr_A++;
                }
                else{
                    // found an element
                    new_flag = 1;
                    break;
                }
            }
            if(new_flag == 1){
                insert_new_element(C, row, col);
            }
        }
    }
}

void insert_new_element(coo_format *C, int row, int col){
    // make sure there is enough space for all points
    if(C->cur_nz + 1 > C->max_nz){
        C->max_nz = C->max_nz + 1000;
        C->coo_I  = (int *)realloc(C->coo_I, sizeof(int)*C->max_nz);
        C->coo_J  = (int *)realloc(C->coo_J, sizeof(int)*C->max_nz);
    }

    C->coo_I[C->cur_nz] = row;
    C->coo_J[C->cur_nz] = col;
    C->cur_nz ++ ;
}