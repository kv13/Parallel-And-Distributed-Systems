#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#include "../inc/utils.h"
#include "../inc/blocking.h"


void blocked_csr_init(blocked_csr *blk_A, csr_format *csr_A, int b){
    
    blk_A->nz      = csr_A->nz;
    blk_A->N       = csr_A->N;
    blk_A->M       = csr_A->M;
    blk_A->b       = b;
    blk_A->n_b     = csr_A->M / b;
}

void blocked_csr_free(blocked_csr *blk_A){
    free(blk_A->blk_ids);
    free(blk_A->blk_nz);
    free(blk_A->csr_col);
    free(blk_A->csr_row);
    free(blk_A);
}

void blocked_csc_init(blocked_csc *blk_B, csc_format *csc_B, int b){
    
    blk_B->nz      = csc_B->nz;
    blk_B->N       = csc_B->N;
    blk_B->M       = csc_B->M;
    blk_B->b       = b;
    blk_B->n_b     = csc_B->N/b;
}

void blocked_csc_free(blocked_csc *blk_B){
    free(blk_B->blk_ids);
    free(blk_B->blk_nz);
    free(blk_B->csc_col);
    free(blk_B->csc_row);
    free(blk_B);
}

void blocked_2_csc(blocked_csc *blk_B, csc_format *csc_BB){
    
    csc_BB->nz = blk_B->nzb;
    csc_BB->N  = blk_B->N/blk_B->b;
    csc_BB->M  = blk_B->M/blk_B->b; 

    csc_BB->csc_col = (uint32_t *)malloc(sizeof(uint32_t)*(csc_BB->M+1));
    csc_BB->csc_row = (uint32_t *)malloc(sizeof(uint32_t)*csc_BB->nz);

    uint32_t *temp_I = (uint32_t *)malloc(sizeof(uint32_t)*csc_BB->nz);
    uint32_t *temp_J = (uint32_t *)malloc(sizeof(uint32_t)*csc_BB->nz);
    
    for(int i=0;i<csc_BB->nz;i++){
        temp_I[i] = blk_B->blk_ids[i]%blk_B->n_b;
        temp_J[i] = blk_B->blk_ids[i]/blk_B->n_b;
    }
    coo2csc(csc_BB->csc_row,csc_BB->csc_col,temp_I,temp_J,csc_BB->nz,csc_BB->M,0);
}

void blocking_csc3(blocked_csc *blk_B, csc_format *csc_B){
    
    block *arr = (block *)malloc(sizeof(block)*blk_B->nz);
    
    int temp_1,temp_2,temp_3;  
    int nzb = 0;

    long int size,t_ptr;

    for(int col=0;col<csc_B->M;col++){
        temp_1 = col/blk_B->b;
        for(int k=csc_B->csc_col[col];k<csc_B->csc_col[col+1];k++){
            temp_2 = csc_B->csc_row[k]/blk_B->b;
            temp_3 = temp_1 *blk_B->n_b + temp_2;
            insert_block( &nzb, arr, temp_3);
        }
    }

    blk_B->nzb     = nzb;
    blk_B->blk_ids = (uint32_t *)malloc(sizeof(uint32_t)*nzb);
    blk_B->blk_nz  = (uint32_t *)calloc(nzb+1, sizeof(uint32_t));

    for(int i=0;i<nzb;i++){
        blk_B->blk_ids[i]  = arr[i].blk_id;
        blk_B->blk_nz[i+1] = blk_B->blk_nz[i] + arr[i].nzs;
    } 
    free(arr);

    size = blk_B->nzb*(blk_B->b+1);

    blk_B->csc_row = (uint32_t *)calloc(blk_B->nz,sizeof(uint32_t));
    blk_B->csc_col = (uint32_t *)calloc(size,sizeof(uint32_t));

    uint32_t *temp_row_idx = (uint32_t *)calloc(blk_B->nzb,sizeof(uint32_t));
    for(int i=0;i<blk_B->nzb;i++)
        temp_row_idx[i] = blk_B->blk_nz[i];

    int row;
    int blk_id,blk_ptr;
    int blk_row_start,blk_col_start;
    int rltv_row,rltv_col;

    for(int col=0;col<blk_B->M;col++){
        temp_1 = col/blk_B->b;
        for(int k=csc_B->csc_col[col];k<csc_B->csc_col[col+1];k++){
            row    = csc_B->csc_row[k];
            temp_2 = row / blk_B->b;

            blk_id  = temp_1*blk_B->n_b + temp_2;
            blk_ptr = binary_search(blk_B->blk_ids,0,blk_B->nzb-1,blk_id);
            if(blk_ptr <0 || blk_ptr >blk_B->nzb-1){
                printf("that should not happened\n");
                exit(1);
            }

            blk_col_start = (blk_id / blk_B->n_b) * blk_B->b; 
            blk_row_start = (blk_id % blk_B->n_b) * blk_B->b;

            rltv_row      = find_relevant_row(row,blk_row_start);
            rltv_col      = find_relevant_col(col,blk_col_start);

            blk_B->csc_row[temp_row_idx[blk_ptr]] = rltv_row;
            temp_row_idx[blk_ptr] = temp_row_idx[blk_ptr] + 1;

            t_ptr = blk_ptr*(blk_B->b+1) + rltv_col + 1;
            blk_B->csc_col[t_ptr] = blk_B->csc_col[t_ptr] + 1; 
        }
    }
}

void blocked_2_csr(blocked_csr *blk_A, csr_format *csr_AA){

    csr_AA->nz = blk_A->nzb;
    csr_AA->N  = blk_A->N/blk_A->b;
    csr_AA->M  = blk_A->M/blk_A->b;

    csr_AA->csr_row = (uint32_t *)malloc(sizeof(uint32_t)*(csr_AA->N+1));
    csr_AA->csr_col = (uint32_t *)malloc(sizeof(uint32_t)*csr_AA->nz);

    uint32_t *temp_I = (uint32_t *)malloc(sizeof(uint32_t)*csr_AA->nz);
    uint32_t *temp_J = (uint32_t *)malloc(sizeof(uint32_t)*csr_AA->nz);

    for(int i=0;i<csr_AA->nz;i++){
        temp_I[i] = blk_A->blk_ids[i]/blk_A->n_b;
        temp_J[i] = blk_A->blk_ids[i]%blk_A->n_b;
    }

    coo2csc(csr_AA->csr_col,csr_AA->csr_row,temp_J,temp_I,csr_AA->nz,csr_AA->N,0);
}

void blocking_csr3(blocked_csr *blk_A, csr_format *csr_A){

    block *arr = (block *)malloc(sizeof(block)*blk_A->nz);
    
    int temp_1,temp_2,temp_3;
    int nzb = 0;
    long int size;
    long int t_ptr;
    
    for(int row =0;row<csr_A->N;row++){
        temp_1 = row/blk_A->b;
        for(int k=csr_A->csr_row[row];k<csr_A->csr_row[row+1];k++){
            temp_2 = csr_A->csr_col[k]/blk_A->b;
            temp_3 = temp_1*blk_A->n_b + temp_2;
            insert_block( &nzb, arr, temp_3);
        }
    }
    blk_A->nzb     = nzb;
    blk_A->blk_ids = (uint32_t *)malloc(sizeof(uint32_t)*nzb);
    blk_A->blk_nz  = (uint32_t *)calloc(nzb+1, sizeof(uint32_t));

    for(int i=0;i<nzb;i++){
        blk_A->blk_ids[i]  = arr[i].blk_id;
        blk_A->blk_nz[i+1] = blk_A->blk_nz[i] + arr[i].nzs;
    }
    free(arr);
    
    size = blk_A->nzb*(blk_A->b+1);
    
    blk_A->csr_col = (uint32_t *)calloc(blk_A->nz,sizeof(uint32_t));
    blk_A->csr_row = (uint32_t *)calloc(size,sizeof(uint32_t));

    uint32_t *temp_col_idx = (uint32_t *)calloc(blk_A->nzb,sizeof(uint32_t));
    
    for(int i=0;i<blk_A->nzb;i++)
        temp_col_idx[i] = blk_A->blk_nz[i];

    int col;
    int blk_id,blk_ptr;
    int blk_row_start,blk_col_start;
    int rltv_row,rltv_col;

    for(int row=0;row<blk_A->N;row++){
        temp_1 = row/blk_A->b;
        for(int k=csr_A->csr_row[row];k<csr_A->csr_row[row+1];k++){
            col     = csr_A->csr_col[k];
            temp_2  = col/blk_A->b;

            blk_id  = temp_1*blk_A->n_b + temp_2;
            blk_ptr = binary_search(blk_A->blk_ids,0,blk_A->nzb-1,blk_id);
            if(blk_ptr <0 || blk_ptr >blk_A->nzb-1){
                printf("that should not happened\n");
                exit(1);
            }

            blk_row_start = (blk_id / blk_A->n_b) * blk_A->b; 
            blk_col_start = (blk_id % blk_A->n_b) * blk_A->b;

            rltv_row      = find_relevant_row(row,blk_row_start);
            rltv_col      = find_relevant_col(col,blk_col_start);

            blk_A->csr_col[temp_col_idx[blk_ptr]] = rltv_col;
            temp_col_idx[blk_ptr] = temp_col_idx[blk_ptr] + 1;
            t_ptr = blk_ptr*(blk_A->b+1) + rltv_row + 1 ;

            blk_A->csr_row[t_ptr] = blk_A->csr_row[t_ptr] + 1;
        }
    }
}

void insert_block(int *nzb, block *nzb_array, int blk_id){
    if(*nzb ==0){
        nzb_array[*nzb].blk_id = blk_id;
        nzb_array[*nzb].nzs    = 1;
        *nzb = *nzb + 1;
    }
    else{
        int exist_flag = 0;
        int pstion_ptr = 0;
        pstion_ptr = search_block(&exist_flag, blk_id, 0, *nzb-1, nzb_array);

        if(exist_flag == 1){
            nzb_array[pstion_ptr].nzs ++;
            return;
        }
        else if(blk_id > nzb_array[pstion_ptr].blk_id){
            for(int i = *nzb-1;i>pstion_ptr;i--){
                nzb_array[i+1] = nzb_array[i];
            }
            nzb_array[pstion_ptr+1].blk_id = blk_id;
            nzb_array[pstion_ptr+1].nzs    = 1;
        }
        else if(blk_id < nzb_array[pstion_ptr].blk_id){
            for(int i=*nzb-1;i>pstion_ptr-1;i--){
                nzb_array[i+1] = nzb_array[i];
            }
            nzb_array[pstion_ptr].blk_id = blk_id;
            nzb_array[pstion_ptr].nzs    = 1;
        }
        *nzb = *nzb + 1;
    }
}

int search_block( int *exist_flag, int blk_id, int l, int r, block *nnzb_array){
    int m=-1;

    while(l<=r){
        
        // find the middle pointer
        m = l+(r-l)/2;

        // check if blk_id exists on nnzb_array
        if(nnzb_array[m].blk_id == blk_id){
            *exist_flag = 1;
            return m;
        }
        else if(nnzb_array[m].blk_id < blk_id){
            l = m + 1;

        }
        else{
            r = m - 1;
        }
    }
    return m;
}

int find_row(int blk_id, int n_b, int b){
    return (blk_id / n_b) * b;
}

int find_col(int blk_id, int n_b, int b){
    return (blk_id % n_b) * b;
}

int find_relevant_row(int current_row,int blk_row_start){
    return current_row - blk_row_start;
}

int find_relevant_col(int current_col, int blk_col_start){
    return current_col - blk_col_start;
}