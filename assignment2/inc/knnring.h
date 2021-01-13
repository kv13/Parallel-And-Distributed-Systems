#ifndef KNNLRING_H
#define KNNLRING_H



//STRUCT DEFINITION
typedef struct{
	int    *nidx;           //!< Indices (0-based) of nearest neighbors [m-by-k]
	double *ndist;          //!< Distance of nearest neighbors          [m-by-k]
	int    m;               //!< Number of query points                 [scalar]
	int    k;               //!< Number of nearest neighbors            [scalar]
}knnresult;


//function to find ,for each point in a query set Y, the k nearest neighbors in the corpus set X
knnresult kNN(double *X,double *Y, int n, int m, int d, int k);

//function to compute distance matrix D
double *compute_D(double *X, double *Y, int n, int m, int d, int k);

//function to find the k nearest neighbors for all the point query Y
void k_select(double *D, double *ndist, int *nidx,int n, int m, int k);

//quick_sort
void quick_sort(double *k_dist, int *k_neigh, int low, int high);

//helper function for quick sort
int partition(double *k_dist, int *k_neigh,int low, int high);

//function to create random points
void create_points(double *X,int size);

//function to write indexes and distances to file
void write_to_file(char *str1, int *nidx, double *ndist, int amount, int k);

//***************************MPI FUNCTIONS ************************************

//mpi interface
knnresult distrAllkNN(double *X, int n, int d, int k);

void knn_combine(knnresult knn_total,knnresult knn_temp, int n, int k, int no_block);

#endif
