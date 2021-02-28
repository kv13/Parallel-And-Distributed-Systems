/**************************************************************
* Implementation of non local means using only cpu,
* to prove the correctness of the algorithm
* and to compare timings.
* The algorithm follows the matlab code was given.
**************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


//declaration of functions

//read-write image functions
float *read_image(char *image_path, int image_size);
void   write_image(char *image_path, float *image, int image_size);


//******non local means functions******
//gaussian filter creation function
float *gaussian_filter(int patch_size, float sigma);
float *non_local_means(float *image, float *g_f, int image_size, int patch_size, float filtS);
float  compute_dist(float *pixel_patch, float *temp_patch, float *g_f, int patch_size, float filtS);
void   find_patch(float *image, float *pixel_patch, int patch_size, int im_size, int row, int col);
float *compute_residue(float *init_image, float *d_image, int image_size);

int main(int argc, char *argv[]){


  //pass all parameters needed for computation.
  if(argc !=6){
    printf("Give correct arguments:1)source file(.csv format), 2) image size, 3)filter sigma, 4)patch size, 5)gaussian filter sigma\n");
    exit(1);
  }

  char *image_path      = argv[1];
  int   image_size      = atoi(argv[2]);
  float f_sigma         = atof(argv[3]);
  int   patch_size      = atoi(argv[4]);
  float g_sigma         = atof(argv[5]);

  //allocate memory for denoised image and for residue
  float *denoised_image ;
  float *res            ;

  //read image
  float *I              = read_image(image_path,image_size);

  //******non local means execution.******
  struct timespec ts_start;
  struct timespec ts_end;

  //compute gaussian kernel
  float *GKERNEL = gaussian_filter(patch_size,g_sigma);

  //compute denoised image
  clock_gettime(CLOCK_MONOTONIC, &ts_start);
  denoised_image = non_local_means(I,GKERNEL,image_size,patch_size,f_sigma);
  clock_gettime(CLOCK_MONOTONIC, &ts_end);

  //*****************************************

  //compute residue
  res = compute_residue(I,denoised_image,image_size);


  //some print messages and store results
  double time = 1000*(double)(ts_end.tv_sec-ts_start.tv_sec) + (double)(ts_end.tv_nsec-ts_start.tv_nsec)/1000000;

  printf("*********************************************************\n");
  printf("sequential execution finished total time is %lfmsec\n",time);
  printf("image size is %d\n", image_size);
  printf("patch size is %d\n", patch_size);
  printf("filtes sigma is %f and gauss kernel sigma is %f\n",f_sigma,g_sigma);
  printf("*********************************************************\n");

  char denoised_path[200];
  snprintf(denoised_path,sizeof(denoised_path),"data/denoised_image_%d_%d.csv",image_size,patch_size);
  write_image(denoised_path,denoised_image,image_size);

  char residue_path[200];
  snprintf(residue_path,sizeof(residue_path),"data/residue_%d_%d.csv",image_size,patch_size);
  write_image(residue_path,res,image_size);

  free(I);
  free(res);
  free(GKERNEL);
  free(denoised_image);
}


//********************************** functions **********************************

float *non_local_means(float *image, float *g_f, int image_size, int patch_size, float filtS){

  float Z;
  float W;
  float estimation;

  float *pixel_patch   = (float *)malloc(patch_size*patch_size*sizeof(float));
  float *temp_patch    = (float *)malloc(patch_size*patch_size*sizeof(float));
  float *d_image       = (float *)malloc(image_size*image_size*sizeof(float));

  for(int i=0;i<image_size;i++){
    for(int j=0;j<image_size;j++){

      find_patch(image, pixel_patch, patch_size, image_size, i, j);
      d_image[i*image_size+j] = -1;
      estimation = 0;
      Z = 0;

      for(int t_i=0; t_i<image_size; t_i++){
        for(int t_j=0; t_j<image_size; t_j++){
          find_patch(image, temp_patch, patch_size, image_size, t_i, t_j);
          W = compute_dist(pixel_patch, temp_patch, g_f, patch_size, filtS);
          Z += W;
          estimation += W*image[t_i*image_size+t_j];
        }
      }

      d_image[i*image_size+j] = estimation/Z;
    }
  }

  return d_image;
}


void   find_patch(float *image, float *pixel_patch, int patch_size, int im_size, int row, int col){

  int boundary = patch_size/2;
  int temp=0;
  int neigh_row;
  int neigh_col;

  //compute the patch. If pointers are out of bound
  //we use mirroring to find it's value
  for(int x=-boundary;x<=boundary;x++){
    for(int y=-boundary;y<=boundary;y++){

      neigh_row = row+x;
      if(neigh_row<0)        neigh_row = abs(neigh_row)-1;
      if(neigh_row>=im_size) neigh_row = im_size-(neigh_row - im_size)-1;

      neigh_col = col+y;
      if(neigh_col<0)        neigh_col = abs(neigh_col)-1;
      if(neigh_col>=im_size) neigh_col = im_size-(neigh_col - im_size)-1;

      pixel_patch[temp] =image[neigh_row*im_size+neigh_col];
      temp++;
    }
  }

}


float  compute_dist(float *pixel_patch, float *temp_patch, float *g_f, int patch_size, float filtS){

  float temp;
  float sum=0;
  float D;

  for(int i=0;i<patch_size;i++){
    for(int j=0;j<patch_size;j++){

      temp = g_f[i*patch_size+j]*g_f[i*patch_size+j]
             *(pixel_patch[i*patch_size+j]-temp_patch[i*patch_size+j])
             *(pixel_patch[i*patch_size+j]-temp_patch[i*patch_size+j]);

      sum += temp;
    }
  }
  D = expf((-sum)/(filtS*filtS));
  return D;
}


float *gaussian_filter(int patch_size, float sigma){

  //allocate memory for the kernel
  float *kernel   = (float *)malloc(patch_size*patch_size*sizeof(float ));

  //compute variance
  float s         = 2.0*sigma*sigma;

  //sum for normalization
  float sum       = 0.0;

  int boundary    = patch_size/2;
  float r;

  for(int x=-boundary; x<=boundary; x++){
    for(int y=-boundary; y<=boundary; y++){
      r = x*x + y*y;
      kernel[(x+boundary)*patch_size+(y+boundary)] = exp(-r/s);
      sum += kernel[(x+boundary)*patch_size+(y+boundary)];
    }
  }

  //find the max and alos divide with the total sum
  float max = -1;
  for(int x=0;x<patch_size;x++){
    for(int y=0;y<patch_size;y++){
      kernel[x*patch_size+y]= kernel[x*patch_size+y]/sum;
      if(kernel[x*patch_size+y]>max) max = kernel[x*patch_size+y];
    }
  }

  //devide also with the max value
  for(int x=0;x<patch_size;x++){
    for(int y=0;y<patch_size;y++){
      kernel[x*patch_size+y]= kernel[x*patch_size+y]/max;
    }
  }


  return kernel;
}


float *compute_residue(float *init_image, float *d_image, int image_size){

  float *r = (float *)malloc(image_size*image_size*sizeof(float));

  for(int i=0;i<image_size;i++){
    for(int j=0;j<image_size;j++){
      r[i*image_size+j] = init_image[i*image_size +j] - d_image[i*image_size+j];
    }
  }

  return r;
}


float *read_image(char *image_path, int image_size){

  //allocate memory for the image
  float *image = (float *)malloc(image_size*image_size*sizeof(float ));

  //open image file
  FILE *file = fopen(image_path,"r");

  //read file
  for(int i=0;i<image_size;i++){
    for(int j=0;j<image_size;j++){
      fscanf(file,"%f,",image+i*image_size+j);
    }
  }

  fclose(file);
  return image;
}


void write_image(char *image_path, float *image, int image_size){

  //open image file
  FILE *file = fopen(image_path,"w+");

  //write file
  for(int i=0;i<image_size;i++){
    for(int j=0;j<image_size;j++){
      fprintf(file,"%f,",image[i*image_size+j]);
    }
    fprintf(file,"\n");
  }

  fclose(file);
}
