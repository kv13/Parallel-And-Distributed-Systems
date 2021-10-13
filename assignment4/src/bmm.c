#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#include "../inc/utils.h"
#include "../inc/blocking.h"
#include "../inc/bmm.h"

void coo_init(coo_format *C, blocked_csr *A, blocked_csc *B){
    
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

void coo_init_filter(coo_format *C, blocked_csr *F){
    
    C->cur_nz = 0;
    C->N      = F->N;
    C->M      = F->M;

    // we know at first the maximum value of non zero elements
    C->max_nz = F->nz;
    
    C->coo_I = (int *)malloc(sizeof(int)*C->max_nz);
    C->coo_J = (int *)malloc(sizeof(int)*C->max_nz);
}

void coo_free(coo_format *C){
    free(C->coo_I);
    free(C->coo_J);
    free(C);
}

void blocked_bmm(coo_format *C, blocked_csr *blk_A, csr_format *AA, blocked_csc *blk_B, csc_format *BB){
    
    // dynamic memory allocation for submatrix C
    // We dont know how many non zero elements has each block has at first
    int t_n = 10000;
    int c_n = 0;
    
    int blk_ptr_A;
    int blk_ptr_B;

    for(int blked_row=0;blked_row<AA->N;blked_row ++){
        for(int blked_col=0;blked_col<BB->M;blked_col++){

            coo2_format *sub_C;
            sub_C = (coo2_format *)malloc(sizeof(coo2_format));
            
            sub_C->max_nz = t_n;
            sub_C->cur_nz = 0;
            sub_C->coo_I = (int *)malloc(sizeof(int)*sub_C->max_nz);
            sub_C->coo_J = (int *)malloc(sizeof(int)*sub_C->max_nz);
            
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

            insert_subcoo_2_coo(C, sub_C->coo_I, sub_C->coo_J,sub_C->cur_nz);
            free(sub_C->coo_I);
            free(sub_C->coo_J);
            free(sub_C);
        }
    }
}

void submatrix_bmm(blocked_csr *blk_A, int A_ptr, blocked_csc *blk_B, int B_ptr, coo2_format *sub_C){
    
    int new_flag = 0;
    int ptr_1, ptr_2;
    int cum_1=0, cum_2;
    
    for(int row=0;row<blk_A->b;row++){
        cum_2 = 0;
        for(int col=0;col<blk_B->b;col++){
            new_flag = 0;
            ptr_1    = 0;
            ptr_2    = 0;
            while(ptr_1<blk_A->csr_row[A_ptr*(blk_A->b+1)+row+1] && ptr_2<blk_B->csc_col[B_ptr*(blk_B->b+1)+col+1]){
                if(blk_A->csr_col[blk_A->blk_nz[A_ptr]+cum_1+ptr_1]>blk_B->csc_row[blk_B->blk_nz[B_ptr]+cum_2+ptr_2]){
                    ptr_2++;
                }
                else if(blk_A->csr_col[blk_A->blk_nz[A_ptr]+cum_1+ptr_1]<blk_B->csc_row[blk_B->blk_nz[B_ptr]+cum_2+ptr_2]){
                    ptr_1++;
                }
                else{
                    new_flag = 1;
                    ptr_1++;
                    ptr_2++; 
                    break;
                }
            }
            
            if(new_flag == 1){
                // insert the element
                int total_col = (blk_B->blk_ids[B_ptr] / blk_B->n_b) * blk_B->b + col;
                int total_row = (blk_A->blk_ids[A_ptr] / blk_A->n_b) * blk_A->b + row;
                
                // search if the element has already been inserted
                int ptr = binary_search2(sub_C->coo_I,sub_C->coo_J, 0, sub_C->cur_nz-1, total_row, total_col);
                
                if(ptr == -1){
                    sub_C->coo_I[0] = total_row;
                    sub_C->coo_J[0] = total_col;
                    sub_C->cur_nz ++;
                }
                else if(sub_C->coo_I[ptr] != total_row || sub_C->coo_J[ptr] != total_col){
                    
                    if(sub_C->cur_nz == sub_C->max_nz){
                        // need to reallocate more memory
                        sub_C->max_nz  = sub_C->max_nz + 10000;
                        sub_C->coo_J = (int *)realloc(sub_C->coo_J, sizeof(int)*sub_C->max_nz);
                        sub_C->coo_I = (int *)realloc(sub_C->coo_I, sizeof(int)*sub_C->max_nz);
                    }

                    if(sub_C->coo_I[ptr] < total_row || (sub_C->coo_I[ptr] == total_row && sub_C->coo_J[ptr]<total_col))
                        ptr ++;

                    for(int i= sub_C->cur_nz-1; i>ptr-1; i--){
                        sub_C->coo_I[i+1] = sub_C->coo_I[i];
                        sub_C->coo_J[i+1] = sub_C->coo_J[i]; 
                    }
    
                    sub_C->coo_I[ptr] = total_row;
                    sub_C->coo_J[ptr] = total_col;
                    sub_C->cur_nz ++;
                }
            }
            cum_2 += blk_B->csc_col[B_ptr*(blk_B->b+1)+col+1];
        }
        cum_1 += blk_A->csr_row[A_ptr*(blk_A->b+1)+row+1];
    }
}

void insert_subcoo_2_coo(coo_format *C, int *sub_I, int *sub_J, int c_n){
    // make sure there is enough space for all points
    if(C->cur_nz + c_n > C->max_nz){
        C->max_nz = C->max_nz + c_n + C->max_nz / 2;
        C->coo_I  = (int *)realloc(C->coo_I, sizeof(int)*C->max_nz);
        C->coo_J  = (int *)realloc(C->coo_J, sizeof(int)*C->max_nz);
    }

    for(int i=0;i<c_n;i++){
        C->coo_I[C->cur_nz] = sub_I[i];
        C->coo_J[C->cur_nz] = sub_J[i];
        C->cur_nz ++ ;
    }
}

void test_results(market_matrix *C_correct, coo_format *C){

    int correct = 0;
    if(C_correct->nz != C->cur_nz){
        printf("correct:%d, mine:%d",C_correct->nz, C->cur_nz);
        printf("Test Failed 1\n");
        return;
    }
    for(int i=0;i<C_correct->nz;i++){
        if(C_correct->coo_I[i] != C->coo_I[i] || C_correct->coo_J[i]!=C->coo_J[i]){
            printf("Test Failesd 2\n");
            return;
        }
    }
    printf("Test Passed!!!\n");
}

// ******************* FILTER FUNCTIONS **************************
void blocked_bmm_filtered(coo_format *C, blocked_csr *blk_F,csr_format *FF,blocked_csr *blk_A, csr_format *AA, blocked_csc *blk_B, csc_format *BB){

    int blk_ptr_A;
    int blk_ptr_B;

    for(int blked_row=0;blked_row<FF->N;blked_row++){
        for(int blked_col_ptr=FF->csr_row[blked_row];blked_col_ptr<FF->csr_row[blked_row+1];blked_col_ptr++){

            coo2_format *sub_C;
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
            insert_subcoo_2_coo(C, sub_C->coo_I, sub_C->coo_J,sub_C->cur_nz);
            free(sub_C->coo_I);
            free(sub_C->coo_J);
            free(sub_C);
        }
    }
    return;
}

void submatrix_bmmfiltered(blocked_csr *blk_F, int F_ptr, blocked_csr *blk_A, int A_ptr, blocked_csc *blk_B, int B_ptr, coo2_format *sub_C){
    int new_flag=0;
    int col;
    int ptr_1, ptr_2;
    int *cum_sum_F, *cum_sum_A, *cum_sum_B;
    
    cum_sum_F = (int *)calloc(blk_F->b+1, sizeof(int));
    cum_sum_A = (int *)calloc(blk_F->b+1, sizeof(int));
    cum_sum_B = (int *)calloc(blk_F->b+1, sizeof(int));

    for(int i=0;i<blk_F->b;i++){
        cum_sum_F[i+1] = cum_sum_F[i] + blk_F->csr_row[F_ptr*(blk_F->b+1)+i+1];
        cum_sum_A[i+1] = cum_sum_A[i] + blk_A->csr_row[A_ptr*(blk_A->b+1)+i+1];
        cum_sum_B[i+1] = cum_sum_B[i] + blk_B->csc_col[B_ptr*(blk_B->b+1)+i+1];
    }

    for(int row=0;row<blk_F->b;row++){
        for(int col_ptr=0; col_ptr<blk_F->csr_row[F_ptr*(blk_F->b+1)+row+1]; col_ptr++){
            col     = blk_F->csr_col[blk_F->blk_nz[F_ptr]+cum_sum_F[row]+col_ptr];
            new_flag = 0;
            
            ptr_1   = 0;
            ptr_2   = 0;
            while(ptr_1 < blk_A->csr_row[A_ptr*(blk_A->b+1)+row+1] && ptr_2 < blk_B->csc_col[B_ptr*(blk_B->b+1)+col+1]){
                if(blk_A->csr_col[blk_A->blk_nz[A_ptr]+cum_sum_A[row]+ptr_1] > blk_B->csc_row[blk_B->blk_nz[B_ptr]+cum_sum_B[col]+ptr_2]){
                    ptr_2++;
                }
                else if(blk_A->csr_col[blk_A->blk_nz[A_ptr]+cum_sum_A[row]+ptr_1] < blk_B->csc_row[blk_B->blk_nz[B_ptr]+cum_sum_B[col]+ptr_2]){
                    ptr_1++;
                }
                else{
                    new_flag = 1;
                    break;
                }
            }
            if(new_flag == 1){
                
                // insert the element
                int total_col   = (blk_B->blk_ids[B_ptr] / blk_B->n_b) * blk_B->b + col; 
                int total_row   = (blk_A->blk_ids[A_ptr] / blk_A->n_b) * blk_A->b + row;
                int ptr = binary_search2(sub_C->coo_I,sub_C->coo_J, 0, sub_C->cur_nz-1,total_row, total_col);
                
                if(ptr==-1){
                    sub_C->coo_I[0] = total_row;
                    sub_C->coo_J[0] = total_col;
                    sub_C->cur_nz++;
                }
                else if(sub_C->coo_I[ptr] != total_row || sub_C->coo_J[ptr] != total_col){
                    
                    if(sub_C->coo_I[ptr] < total_row || (sub_C->coo_I[ptr] == total_row && sub_C->coo_J[ptr]<total_col))
                        ptr ++;
                    
                    for(int i= sub_C->cur_nz-1; i>ptr-1; i--){
                        sub_C->coo_I[i+1] = sub_C->coo_I[i];
                        sub_C->coo_J[i+1] = sub_C->coo_J[i]; 
                    }

                    sub_C->coo_I[ptr] = total_row;
                    sub_C->coo_J[ptr] = total_col;
                    sub_C->cur_nz ++;
                }
            }
        }
    }
}