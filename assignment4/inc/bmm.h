#ifndef BMM_H
#define BMM_H

// ******* struct definitions *******
typedef struct 
{
    int *coo_I, *coo_J;
    int N;      // number of rows
    int M;      // number of columns
    int max_nz; // max number of non zero elements
    int cur_nz; // current number of non zero elements
}coo_format;

// to do at the end =>replace this struct with the above
typedef struct 
{
    int *coo_I, *coo_J;
    int max_nz; // max number of non zero elements
    int cur_nz; // current number of non zero elements
}coo2_format;


// ******* functions declaration *******

// ******* blocked bmm functions
void coo_init(coo_format *C, blocked_csr *A, blocked_csc *B);
void coo_free(coo_format *C);
void blocked_bmm(coo_format *C, blocked_csr *blk_A, csr_format *AA, blocked_csc *blk_B, csc_format *BB);
void submatrix_bmm(blocked_csr *blk_A, int A_ptr, blocked_csc *blk_B, int B_ptr, coo2_format *sub_C);
void insert_subcoo_2_coo(coo_format *C, int *sub_I, int *sub_J, int c_n);

// ******* filtered blocked bmm functions
void coo_init_filter(coo_format *C, blocked_csr *F);
void blocked_bmm_filtered(coo_format *C, blocked_csr *blk_F,csr_format *FF,blocked_csr *blk_A, csr_format *AA, blocked_csc *blk_B, csc_format *BB);
void submatrix_bmmfiltered(blocked_csr *blk_F, int F_ptr, blocked_csr *blk_A, int A_ptr, blocked_csc *blk_B, int B_ptr, coo2_format *sub_C);

void test_results(market_matrix *C_correct, coo_format *C);
#endif