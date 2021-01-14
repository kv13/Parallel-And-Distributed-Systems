//%%%%% dumb read datasets functions
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../inc/readlib.h"

double *read_1(char *file_path, int *n, int *d){
  if(strstr(file_path,"ColorH")!=NULL){
    *n=68040;
  	*d=32;
  }else if(strstr(file_path,"ColorM")!=NULL){
    *n = 68040;
    *d = 9;
  }else{
    printf("ONLY ColorHistogram and ColorMoments read\n");
    exit(1);
  }
  double *X;
  X = (double *)malloc((*n) * (*d) * sizeof(double));

  //open the file
  FILE *fp;
  fp = fopen(file_path,"r");
	if(fp == NULL){
		printf("CANNOT OPEN THE FILE...EXITING\n");
		exit(1);
	}


  //initialize variables
  char *line = NULL;
	size_t  len = 0;
	ssize_t read;

	char delim[] = " ";
	char *ptr ;

  //READ THE FILE
  int counter=0;
	while(read = getline(&line, &len, fp)!= -1){
		ptr = strtok(line,delim);

		//skip first element which is the index;
		ptr = strtok(NULL,delim);

		while(ptr!=NULL){

			X[counter] = atof(ptr);
			counter ++;
			ptr = strtok(NULL,delim);

		}
	}
  return X;
}


double *read_2(char *file_path, int *n, int *d){

  //set the number of points and the dimensions
  *n = 130065;
	*d = 50;

  double *X;
  X = (double *)malloc((*n) * (*d) * sizeof(double));

  //open the file
  FILE *fp;
	fp = fopen(file_path,"r");
	if(fp == NULL){
		printf("CANNOT OPEN THE FILE...EXITING\n");
		exit(1);
	}

  //initialize variables
  char *line = NULL;
	size_t  len = 0;
	ssize_t read;

	char delim[] = "  ";
	char *ptr ;

  //READ data

	//skip the first line
	read = getline(&line,&len,fp);
	//printf("%s\n",line);

  int counter = 0;
	while(read = getline(&line,&len,fp)!=-1){
		ptr = strtok(line,delim);
		while(ptr!=NULL){
			X[counter] = atof(ptr);
			counter++;
			ptr = strtok(NULL,delim);
		}
	}

  return X;
}


double *read_3(char *file_path, int *n, int *d){
  //set the number of points and the dimensions
  *n = 106574;
	*d = 518;

  double *X;
  X = (double *)malloc((*n) * (*d) * sizeof(double));

  //open the file
  FILE *fp;
  fp = fopen(file_path,"r");
  if(fp == NULL){
    printf("CANNOT OPEN THE FILE...EXITING\n");
    exit(1);
  }

  //initialize variables
  char    *line = NULL;
	size_t  len = 0;
	ssize_t read;

	char delim[] = ",";
	char *ptr;

  //READ data

  //skip the first 4 line
	read = getline(&line,&len,fp);
	//printf("%s\n",line);

	read = getline(&line,&len,fp);
	//printf("%s\n",line);

	read = getline(&line,&len,fp);
	//printf("%s\n",line);

	read = getline(&line,&len,fp);
	//printf("%s\n",line);

  int counter = 0;
	while(read = getline(&line,&len,fp)!=-1){

		ptr = strtok(line,delim);

		//skip the first element which is the index
		ptr = strtok(NULL,delim);
		while(ptr != NULL){
			X[counter] = atof(ptr);
			counter ++;
			ptr = strtok(NULL,delim);
		}
	}

  return X;
}


double *read_4(char *file_path, int *n, int *d){
  if    (strstr(file_path,"BBC") != NULL)       *n = 17720;
  else if    (strstr(file_path,"CNN.") != NULL) *n = 22545;
  else if    (strstr(file_path,"CNNI") != NULL) *n = 33117;
  else if   (strstr(file_path,"NDTV") != NULL)  *n = 17051;
  else if  (strstr(file_path,"TIMES") != NULL)  *n = 39252;

  *d = 17;
  int d_min = 17;

  double *X;
  X = (double *)malloc((*n) * (*d) * sizeof(double));

  //open file
  FILE *fp;
	fp = fopen(file_path,"r");
	if(fp == NULL){
		printf("CANNOT OPEN THE FILE...EXITING\n");
		exit(1);
	}

  //initialize variables
  char    *line = NULL;
	size_t  len = 0;
	ssize_t read;

	char delim[] = ":";
	char delim_2[]= "  ";
	char delim_3[] = " ";
	char *ptr;
	int counter = 0;
	int d_counter;

  while(read = getline(&line,&len,fp)!=-1){

		ptr = strtok(line,delim);
		d_counter = 0;
		while(ptr != NULL){
			ptr = strtok(NULL,delim_3);
			//printf("%s\n",ptr);
			X[counter] = atof(ptr);
			counter ++;
			d_counter ++;
			if(d_counter >= d_min)break;
			ptr = strtok(NULL,delim);
		}
	}
  return X;
}

double *read_data(char *file_path, int *n, int *d){

  double *X;

  if         (strstr(file_path, "Co") != NULL)  X = read_1(file_path,n,d);
  else if  (strstr(file_path, "Mini") != NULL)  X = read_2(file_path,n,d);
  else if  (strstr(file_path, "feat") != NULL)  X = read_3(file_path,n,d);
  else if    (strstr(file_path,"BBC") != NULL)  X = read_4(file_path,n,d);
  else if    (strstr(file_path,"CNN") != NULL)  X = read_4(file_path,n,d);
  else if   (strstr(file_path,"NDTV") != NULL)  X = read_4(file_path,n,d);
  else if  (strstr(file_path,"TIMES") != NULL)  X = read_4(file_path,n,d);
  else{
    printf("CANNOT RECOGNIZE THE FILE...=>EXITING\n");
    exit(1);
  }

  return X;
}
