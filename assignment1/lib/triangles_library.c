#include "triangles_library.h"


//function to transform coo format to csc format
void coo2csc(
  uint32_t       * const row,       /*!< CSC row start indices */
  uint32_t       * const col,       /*!< CSC column indices */
  uint32_t const * const row_coo,   /*!< COO row indices */
  uint32_t const * const col_coo,   /*!< COO column indices */
  uint32_t const         nnz,       /*!< Number of nonzero elements */
  uint32_t const         n,         /*!< Number of rows/columns */
  uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
) {

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


//V3 algorithm for computing triangles.
void V3_algorithm(uint32_t const * const csc_col,uint32_t const * const csc_row, uint32_t * const c3, int n){
  uint32_t i,k,j;
  for(i = 0;i<n;i++){
    for(uint32_t temp1 = 0;temp1<csc_col[i+1]-csc_col[i];temp1++){
      j = csc_row[csc_col[i]+temp1];
      if(j<i+1)continue;
      for(uint32_t temp2=0;temp2<csc_col[j+1]-csc_col[j];temp2++){
        k = csc_row[csc_col[j]+temp2];
        if(k<j+1)continue;
        if(binary_search(k,i,temp1,csc_col,csc_row) !=0 ){
          c3[i]++;
          c3[j]++;
          c3[k]++;
        }
      }
    }
  }
}

//binary search to find if exist edge from node k to node i.
uint32_t binary_search(uint32_t k, uint32_t i, uint32_t t_j, uint32_t const * const csc_col,uint32_t const * const csc_row){
	uint32_t first,last;
	uint32_t middle;
	first = t_j+1;
	last = csc_col[i+1]-csc_col[i]-1;
	middle = first+(last-first)/2;
	while(first<=last){
		if(csc_row[csc_col[i]+middle] < k){
			first = middle + 1;
		}
		else if(csc_row[csc_col[i]+middle] == k){
			return middle;
		}
		else{
			last = middle - 1;
		}
		middle = first+(last-first)/2;
	}
  return 0;
}


//Hadamard algorithm
void Hadamard_algorithm(uint32_t const * const csc_col, uint32_t const * const csc_row, uint32_t * const c3, uint32_t * const val, int n){
  for (uint32_t temp1=0;temp1<n;temp1++){
    for (uint32_t temp2=0;temp2<csc_col[temp1+1]-csc_col[temp1];temp2++){
      uint32_t i = temp1;
      uint32_t j = csc_row[csc_col[temp1]+temp2];
      if(j<i) val[csc_col[temp1]+temp2] = val[csc_col[j]+binary_search(i,j,0,csc_col,csc_row)];
      else val[csc_col[temp1]+temp2]=compute_sum(temp1,j,csc_col,csc_row);
      c3[temp1] = c3[temp1]+val[csc_col[temp1]+temp2];
    }
    c3[temp1]=c3[temp1]/2;
  }
}


//algorithm to merge to sorted list in O(N+M) complexity
uint32_t compute_sum(uint32_t const i, uint32_t const j, uint32_t const * const csc_col, uint32_t const * const csc_row){
	uint32_t sum=0;
	uint32_t size_i = csc_col[i+1]-csc_col[i];
	uint32_t size_j = csc_col[j+1]-csc_col[j];
	uint32_t counter_i=0;
	uint32_t counter_j=0;
	uint32_t counter=0;
	uint32_t end =0;
	for(;;){

		if(csc_row[csc_col[i]+counter_i] == csc_row[csc_col[j]+counter_j]){
			sum++;
			counter_i++;
			if(counter_i ==size_i){
				end =1;
			}
			counter_j++;
			if(counter_j == size_j){
				end =1;
			}
			counter = counter + 2;
			if(counter == size_i + size_j){
				end = 1;
			}
		}
		else if(csc_row[csc_col[i]+counter_i]<csc_row[csc_col[j]+counter_j]){
			counter_i++;
			counter ++;
			if(counter_i == size_i){
				end=1;
			}
		}
		else{
			counter_j++;
			counter++;
			if(counter_j == size_j){
				end=1;
			}
		}
		if(counter == size_i+size_j){
			end =1;
		}
		if(end == 1)break;
	}
	return sum;
}


//write to file function
void write_to_file(char *str, uint32_t const * const c3, int n, double u_sec){
	FILE *fp;
	fp = fopen(str,"w+");
	if(fp == NULL){
		fprintf(stderr,"ERROR: CANNOT OPEN FILES FOR WRITING THE RESULTS...EXITING\n");
		exit(1);
	}
	for(uint32_t i=0;i<n;i++){
		fprintf(fp,"id = %u, c3 = %u \n",i,c3[i]);
	}
	fprintf(fp,"TOTAL TIME in microseconds: %lf u_sec",u_sec);
	fclose(fp);
}
