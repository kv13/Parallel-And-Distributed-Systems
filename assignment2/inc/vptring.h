#ifndef VPTRING_H
#define VPTRING_H

#define MAX_POINTS 100000

typedef struct{
	int    *leaf_counter;
	int    *tree_counter;
	int    *indexes;
	double *tag;
  double *coords;
  double *mu;
  double *mu_points;
	double step;
}tree_struct;

//function to initalize variables for the tree struct
tree_struct create_tree(double *X, int *nidx,int n, int d, int B);

//function to compute and create the tree
void create_subtree(tree_struct tree,double *X, int *nidx, int n, int d, int B, int height, int n_total);

//function to calculate distances through vantage point and the rest corpus points
void calculate_distances(double *distances, int vantage_point, int n, int *nidx, double *X, int d);

//quick select function
double quick_select(double *arr, int left, int right, int pos);

//partition function
int partition_1(double *arr, int left, int right);

//swap function
void swap(double *a, double *b);

//write to file indexes median and points for every vantage point
void write_to_file_3(char *str, int *array,double *array2,double *array3, int amount);

//############################ SEARCH FUNCTIONS##############################

//insert point of searching
knnresult search_k_nearest(double *X, int *indexes, double *mu, double *Y, int n, int m, int d, int k);

//function to calculate k nearest.
void insert_point_vpt(int *nidx, double *ndist,int point, int k, int *num_neigh, int new_point, double dist);

//search function
void search_tree_1(double *X, int *indexes, double *mu, double *Y, int *nidx,double *ndist, int n, int point, int d, int k, int height, int *num_neigh);

//######################### MPI FUNCTIONS ###################################
knnresult distrV2(double *X, int n,int d, int k, int B);


#endif
