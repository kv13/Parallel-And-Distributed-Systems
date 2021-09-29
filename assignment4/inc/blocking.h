#ifndef BLOCKING_H
#define BLOCKING_H


// ******* struct definitions *******
typedef struct{
    uint32_t *blk_nz;
    uint32_t *blk_ids;
    uint32_t *csc_col;
    uint32_t *csc_row;
    uint32_t nz, nzb;
    uint32_t N, M, n_b, b;
} blocked_csc;

typedef struct{
    uint32_t *blk_nz;
    uint32_t *blk_ids;
    uint32_t *csr_col;
    uint32_t *csr_row;
    uint32_t nz, nzb;
    uint32_t N, M, n_b, b;
} blocked_csr;

typedef struct{
    int blk_id;
    int nzs;
}block;

// ******* functions declaration *******
void blocked_csr_init(blocked_csr *blk_A, csr_format *csr_A, int b);
void blocked_csr_free(blocked_csr *blk_A);

void blocked_csc_init(blocked_csc *blk_B, csc_format *csc_B, int b);
void blocked_csc_free(blocked_csc *blk_B);

void blocking_csr3(blocked_csr *blk_A, csr_format *csr_A);
void blocked_2_csr(blocked_csr *blk_A, csr_format *csr_AA);

void blocking_csc3(blocked_csc *blk_B, csc_format *csc_B);
void blocked_2_csc(blocked_csc *blk_B, csc_format *csc_BB);

void insert_block(int *nzb, block *nzb_array, int blk_id);
int  search_block(int *exist_flag, int blk_id, int l, int r, block *nnzb_array);

int  find_row(int blk_id, int n_b, int b);
int  find_col(int blk_id, int n_b, int b);
int  find_relevant_row(int current_row,int blk_row_start);
int  find_relevant_col(int current_col, int blk_col_start);

#endif