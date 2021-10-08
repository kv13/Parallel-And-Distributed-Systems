#ifndef UTILS_H
#define UTILS_H

// ******* struct definitions *******

// struct to read coo matrices from maret matrix
typedef struct 
{
    int *coo_I, *coo_J;
    int N;  // number of rows
    int M;  // number of columns
    int nz; // number of non zero elements

} market_matrix;

typedef struct 
{
    uint32_t *csr_col;
    uint32_t *csr_row;
    uint32_t N; // number of rows
    uint32_t M; // number of columns
    uint32_t nz;
} csr_format;

typedef struct 
{
    uint32_t *csc_col;
    uint32_t *csc_row;
    uint32_t N; // number of rows
    uint32_t M; // number of columns
    uint32_t nz;
} csc_format;

// ******* functions declaration *******
int   binary_search(uint32_t *arr , int l, int r, int s_obj);
int  binary_search2(uint32_t *arr1, uint32_t *arr2, int l, int r, int s_obj1, int s_obj2);
int  partition(int *arr_1, int *arr_2, int l, int h);
void read_mm_matrices(market_matrix *mtx, char *name);
void coo2csc(uint32_t *const row, uint32_t *const col, uint32_t const *const row_coo, uint32_t const *const col_coo, uint32_t const nnz, uint32_t const n, uint32_t const isOneBased);
void csr_init(csr_format *mtx1, market_matrix *mtx2);
void csr_init2(csr_format *mtx1, int N, int M, int nz);
void csr_free(csr_format *mtx1);
void csc_init(csc_format *mtx1, market_matrix *mtx2);
void csc_free(csc_format *mtx1);
void visualization_csc(uint32_t *csc_row,uint32_t *csc_col, int N, int M);
void visualization_csr(uint32_t *csr_row,uint32_t *csr_col, int N, int M);
void quick_sort(int *arr_1, int *arr_2, int l, int h);
void swap(int *a, int *b);

#endif