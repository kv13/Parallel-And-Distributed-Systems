#ifndef READLIB_H
#define READLIB_H


//read 1st dataset
double *read_1(char *file_path, int *n, int *d);

//read 2nd dataset
double *read_2(char *file_path, int *n, int *d);

//read 3rd dataset
double *read_3(char *file_path, int *n, int *d);

//read 4th dataset
double *read_4(char *file_path, int *n, int *d);

//main func to read files
double *read_data(char *str, int *n, int *d);

#endif
