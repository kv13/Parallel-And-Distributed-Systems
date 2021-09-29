#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "../inc/mmio.h"
#include "../inc/utils.h"

// initialize csr matrices
void csr_init(csr_format *mtx1, market_matrix *mtx2){
    
    mtx1->csr_col = (uint32_t *)malloc(sizeof(uint32_t)*mtx2->nz);
    mtx1->csr_row = (uint32_t *)malloc(sizeof(uint32_t)*(mtx2->N+1));

    mtx1->N  = mtx2->N;
    mtx1->M  = mtx2->M;
    mtx1->nz = mtx2->nz;
}

void csr_free(csr_format *mtx1){
    free(mtx1->csr_col);
    free(mtx1->csr_row);
    free(mtx1);
}

// initialize csc matrices
void csc_init(csc_format *mtx1, market_matrix *mtx2){
    
    mtx1->csc_row = (uint32_t *)malloc(sizeof(uint32_t)*mtx2->nz);
    mtx1->csc_col = (uint32_t *)malloc(sizeof(uint32_t)*(mtx2->M+1));

    mtx1->N  = mtx2->N;
    mtx1->M  = mtx2->M;
    mtx1->nz = mtx2->nz;
}

void csc_free(csc_format *mtx1){
    free(mtx1->csc_col);
    free(mtx1->csc_row);
    free(mtx1);
}

void quick_sort(int *arr_1, int *arr_2, int l, int h){
    if(l<h){
        int pi = partition(arr_1, arr_2, l, h);

        quick_sort(arr_1, arr_2, l, pi-1);
        quick_sort(arr_1, arr_2, pi+1, h);
    }

}

int partition(int *arr_1, int *arr_2, int l, int h){
    
    int pivot_1 = arr_1[h];
    int pivot_2 = arr_2[h];

    int i = l - 1;

    for(int j=l; j<=h; j++){
        if(arr_1[j]<pivot_1 || (arr_1[j]==pivot_1 && arr_2[j]<pivot_2)){
            i++;
            swap(&arr_1[i],&arr_1[j]);
            swap(&arr_2[i],&arr_2[j]);
        }

    }
    swap(&arr_1[i+1],&arr_1[h]);
    swap(&arr_2[i+1],&arr_2[h]);
    
    return (i+1);
}

void swap(int *a, int *b){
    int t = *a;
    *a = *b;
    *b = t;
}

int binary_search(uint32_t *arr , int l, int r, int s_obj){
   
    int mid = -1;
    
    while(l<=r){
        
        mid = l+(r-l)/2;

        if(arr[mid] == s_obj)
            return mid;

        if(arr[mid] < s_obj)
            l = mid+1;
        
        else
            r = mid-1;
    }

    return mid;
}

int  binary_search2(uint32_t *arr1, uint32_t *arr2, int l, int r, int s_obj1, int s_obj2){
    int mid = -1;

    while(l<=r){
        
        mid = l+(r-l)/2;

        if(arr1[mid] == s_obj1 && arr2[mid] == s_obj2)
            return mid;

        if(arr1[mid] < s_obj1 || (arr1[mid] == s_obj1 && arr2[mid] < s_obj2))
            l = mid+1;
        
        else
            r = mid-1;
    }
    return mid;
}
// visualize csc/csr format in order to make sure that matrix is properly readen
void visualization_csc(uint32_t *csc_row, uint32_t *csc_col, int N, int M){
    
    int8_t **array;
    array = (int8_t **)malloc(N*sizeof(int8_t *));
    for(int i=0;i<N;i++){
        array[i] = (int8_t *)malloc(M*sizeof(int8_t));
        for(int j=0;j<M;j++){
            array[i][j]=0;
        }
    }

    for(int j=0;j<M;j++){
        // find non zero elements for column j
        for(int i=csc_col[j];i<csc_col[j+1];i++){
            array[csc_row[i]][j]=1;
        }
    }

    for(int i=0;i<N;i++){
        for(int j=0;j<M;j++){
            printf("%d,",array[i][j]);
        }
        printf("\n");
    }

}

void visualization_csr(uint32_t *csr_row,uint32_t *csr_col, int N, int M){
    
    int8_t **array;
    array = (int8_t **)malloc(N*sizeof(int8_t *));
    for(int i=0;i<N;i++){
        array[i] = (int8_t *)malloc(M*sizeof(int8_t));
        for(int j=0;j<M;j++){
            array[i][j]=0;
        }
    }

    for(int i=0;i<N;i++){
        for(int j=csr_row[i];j<csr_row[i+1];j++){
            array[i][csr_col[j]]=1;
        }
    }
    
    for(int i=0;i<N;i++){
        for(int j=0;j<M;j++){
            printf("%d,",array[i][j]);
        }
        printf("\n");
    }
}
//function to read matrices from matrix market
void read_mm_matrices(market_matrix *mtx, char *name){

    
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz; 
    int i, *I, *J;  
    
    if ((f = fopen(name, "r")) == NULL) 
        exit(1);

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0){
        printf("Cannot read file\n");
        exit(1);
    }
        


    /* reseve memory for matrices */

    I = (int *) malloc(nz * sizeof(int));
    J = (int *) malloc(nz * sizeof(int));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    for (i=0; i<nz; i++)
    {
        fscanf(f, "%d %d \n", &I[i], &J[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }

    if (f !=stdin) fclose(f);

    /************************/
    /* now write out matrix */
    /************************/

    /*
    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);
    for (i=0; i<nz; i++)
        fprintf(stdout, "%d %d \n", I[i]+1, J[i]+1);
    */
   
   mtx->coo_I = I;
   mtx->coo_J = J;
   mtx->M     = M;
   mtx->N     = N;
   mtx->nz    = nz;
}

//function to transform coo format to csc/csr format
void coo2csc( uint32_t *const row, uint32_t *const col, uint32_t const *const row_coo, uint32_t const *const col_coo, uint32_t const nnz, uint32_t const n, uint32_t const isOneBased){

  // ----- cannot assume that input is already 0!
  for (uint32_t l = 0; l < n+1; l++) col[l] = 0;


  // ----- find the correct column sizes
  for (uint32_t l = 0; l < nnz; l++)
    col[col_coo[l] - isOneBased]++;

  // ----- cumulative sum
  for (uint32_t i = 0, cumsum = 0; i < n; i++) {
    uint32_t temp = col[i];
    col[i] = cumsum;
    cumsum += temp;
  }
  col[n] = nnz;
  // ----- copy the row indices to the correct place
  for (uint32_t l = 0; l < nnz; l++) {
    uint32_t col_l;
    col_l = col_coo[l] - isOneBased;

    uint32_t dst = col[col_l];
    row[dst] = row_coo[l] - isOneBased;

    col[col_l]++;
  }
  // ----- revert the column pointers
  for (uint32_t i = 0, last = 0; i < n; i++) {
    uint32_t temp = col[i];
    col[i] = last;
    last = temp;
  }
}