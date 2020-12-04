#ifndef TRIANGLES_H_INCLUDED
#define TRIANGLES_H_INCLUDED

//gather all functions are necessary for
//counting number of triangles.
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

//coo2csc from github.
void coo2csc(
  uint32_t       * const row,       /*!< CSC row start indices */
  uint32_t       * const col,       /*!< CSC column indices */
  uint32_t const * const row_coo,   /*!< COO row indices */
  uint32_t const * const col_coo,   /*!< COO column indices */
  uint32_t const         nnz,       /*!< Number of nonzero elements */
  uint32_t const         n,         /*!< Number of rows/columns */
  uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
);


//V3 ALGORITHM
void V3_algorithm(uint32_t const * const csc_col,uint32_t const * const csc_row, uint32_t * const c3, int n);


//binary search to find if an edge exists between two nodes
uint32_t binary_search(uint32_t k, uint32_t i, uint32_t t_j, uint32_t const * const csc_col,uint32_t const * const csc_row);


//Hadamard algorithm
void Hadamard_algorithm(uint32_t const * const csc_col, uint32_t const * const csc_row, uint32_t * const c3, uint32_t * const val, int n);


//algorithm to merge to sorted list in O(N+M) complexity
uint32_t compute_sum(uint32_t const i, uint32_t const j, uint32_t const * const csc_col, uint32_t const * const csc_row);


//write to file function
void write_to_file(char *str,uint32_t const * const c3,int n,double u_sec);
#endif
