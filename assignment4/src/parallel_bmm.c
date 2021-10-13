#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>
#include <time.h>

#include "../inc/utils.h"
#include "../inc/blocking.h"
#include "../inc/bmm.h"
#include "../inc/parallel_bmm.h"

void blocked_bmm_shared(coo_format *C, blocked_csr *blk_A, csr_format *AA, blocked_csc *blk_B, csc_format *BB){

    // dynamic memory allocation for submatrix C
    // We dont know how many non zero elements has each block has at first
    int t_n = 10000;
    int c_n = 0;
    
    int blk_ptr_A, blk_ptr_B;
    int blked_row, blked_col;
    coo2_format *sub_C;

    
    for(blked_row=0;blked_row<AA->N;blked_row ++){
        #pragma omp parallel shared(blked_row, C,blk_A, AA, blk_B, BB) private(sub_C,blk_ptr_A,blk_ptr_B,blked_col)
        {
            #pragma omp for schedule(dynamic) nowait
            for(blked_col=0;blked_col<BB->M;blked_col++){

                sub_C = (coo2_format *)malloc(sizeof(coo2_format));
                
                sub_C->max_nz = t_n;
                sub_C->cur_nz = 0;
                sub_C->coo_I  = (int *)malloc(sizeof(int)*sub_C->max_nz);
                sub_C->coo_J  = (int *)malloc(sizeof(int)*sub_C->max_nz);
                
                blk_ptr_A = AA->csr_row[blked_row];
                blk_ptr_B = BB->csc_col[blked_col];
                
                while(blk_ptr_A < AA->csr_row[blked_row + 1] && blk_ptr_B < BB->csc_col[blked_col+1]){
                    if(AA->csr_col[blk_ptr_A] > BB->csc_row[blk_ptr_B]){
                        blk_ptr_B++;
                    }
                    else if(AA->csr_col[blk_ptr_A]<BB->csc_row[blk_ptr_B]){
                        blk_ptr_A++;
                    }
                    else{
                        submatrix_bmm(blk_A, blk_ptr_A, blk_B, blk_ptr_B, sub_C);
                        blk_ptr_A++;
                        blk_ptr_B++;
                    }
                }
                #pragma omp critical(dataupdate)
                {
                    insert_subcoo_2_coo(C, sub_C->coo_I, sub_C->coo_J,sub_C->cur_nz);
                }
                free(sub_C->coo_I);
                free(sub_C->coo_J);
                free(sub_C);
            }
        }
    }
}

void blocked_bmm_filtered_shared(coo_format *C, blocked_csr *blk_F,csr_format *FF,blocked_csr *blk_A, csr_format *AA, blocked_csc *blk_B, csc_format *BB){
    
    int blk_ptr_A, blk_ptr_B;
    int blked_row, blked_col_ptr;
    coo2_format *sub_C;

    for(blked_row=0;blked_row<FF->N;blked_row++){
        #pragma omp parallel shared(blked_row, C, blk_F, FF, blk_A, AA, blk_B, BB) private(sub_C, blk_ptr_A, blk_ptr_B, blked_col_ptr)
        {
            #pragma omp for schedule(dynamic) nowait
            for(blked_col_ptr=FF->csr_row[blked_row];blked_col_ptr<FF->csr_row[blked_row+1];blked_col_ptr++){
                
                sub_C = (coo2_format *)malloc(sizeof(coo2_format));
                sub_C->max_nz = blk_F->blk_nz[blked_col_ptr+1] - blk_F->blk_nz[blked_col_ptr];
                sub_C->cur_nz = 0;
                sub_C->coo_I  = (int *)malloc(sizeof(int) * sub_C->max_nz);
                sub_C->coo_J  = (int *)malloc(sizeof(int) * sub_C->max_nz);

                blk_ptr_A = AA->csr_row[blked_row];
                blk_ptr_B = BB->csc_col[FF->csr_col[blked_col_ptr]];

                while(blk_ptr_A < AA->csr_row[blked_row + 1] && blk_ptr_B < BB->csc_col[FF->csr_col[blked_col_ptr]+1]){
                    if(AA->csr_col[blk_ptr_A] > BB->csc_row[blk_ptr_B]){
                        blk_ptr_B++;
                    }
                    else if(AA->csr_col[blk_ptr_A]<BB->csc_row[blk_ptr_B]){
                        blk_ptr_A++;
                    }
                    else{
                        submatrix_bmmfiltered(blk_F, blked_col_ptr, blk_A, blk_ptr_A, blk_B, blk_ptr_B, sub_C);
                        blk_ptr_A++;
                        blk_ptr_B++;
                    }
                }

                #pragma omp critical(dataupdate_filtered)
                {
                    insert_subcoo_2_coo(C, sub_C->coo_I, sub_C->coo_J,sub_C->cur_nz);
                }

                free(sub_C->coo_I);
                free(sub_C->coo_J);
                free(sub_C);
            }
        }
    }
}
