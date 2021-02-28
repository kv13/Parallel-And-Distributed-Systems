/**************************************************************
* Implementation of non local means using gpu,
* In this implementation shared memory is used
* to achieve the maximum speed up
**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PSIZE 3

//declaration of functions

//****host functions****
float **read_image(char *im_path,  int im_size);
void    write_image2(char *d_path , float  *im, int im_size);
float  *gaussian_filter(int patch_size, float sigma);
float  *to_rowmajor(float **image, int im_size);
float **add_padding(float **image, int im_size, int p_size, int *pad_size);
void   checkCuda(cudaError_t result);

//****device functions****
__global__ void filter_kernel(float *image, float *denoised_image, float *gauss_kernel, int p_size, int im_size, float f_s );
__device__ void find_patch(float *pixel_patch, int current_pixel,int im_size,int p_size,float *image);
__device__ float compute_weight(float *pixel_patch,float *temp_patch,float *g_kernel, int p_size, float f_sigma);

int main(int argc, char *argv[]){

  //pass all parameters needed for computation.
  //patch size now is defined because must be known in compile time
  if(argc !=5){
    printf("Give correct arguments:1)source file(.csv format), 2) image size, 3)filter sigma, 4)gaussian filter sigma\n");
    exit(1);
  }

  char *image_path      = argv[1];
  int   image_size      = atoi(argv[2]);
  float f_sigma         = atof(argv[3]);
  float g_sigma         = atof(argv[4]);
  int   patch_size      = PSIZE;
  float *den_image      = (float *)malloc(image_size*image_size*sizeof(float));

  //read 2-D image
  float **I = read_image(image_path,image_size);

  //add padding
  int padded_size;
  float **D=add_padding(I,image_size,patch_size,&padded_size);

  //the initial 2-D image is now useless
  for(int i=0;i<image_size;i++)free(I[i]);
  free(I);

  //write image to rowmajor format
  float *padded_I = to_rowmajor(D,padded_size);

  //the padded 2-D image is now useless
  for(int i=0;i<padded_size;i++)free(D[i]);
  free(D);

  //compute gaussian filte
  float *GKERNEL = gaussian_filter(patch_size,g_sigma);

  //******** gpu section ********

  //define block and grid size
  int block_size = 256;
  int grid_size = image_size*image_size/block_size;

  //device variables
  float *dev_padded_image, *dev_denoised_image, *dev_g_kernel;

  //variables to compute time of execution.
  cudaEvent_t start,stop;

  checkCuda(cudaEventCreate(&start));
  checkCuda(cudaEventCreate(&stop));

  checkCuda(cudaMalloc((void **)&dev_g_kernel,patch_size*patch_size*sizeof(float)));
  checkCuda(cudaMemcpy(dev_g_kernel,GKERNEL,patch_size*patch_size*sizeof(float),cudaMemcpyHostToDevice));

  checkCuda(cudaMalloc((void **)&dev_padded_image,padded_size*padded_size*sizeof(float)));
  checkCuda(cudaMemcpy(dev_padded_image,padded_I,padded_size*padded_size*sizeof(float),cudaMemcpyHostToDevice));

  checkCuda(cudaMalloc((void **)&dev_denoised_image,image_size*image_size*sizeof(float)));

  /*************************nlm execution*************************/
  //declare the size for shared memory.
  int shared_mem_size = padded_size*(block_size/image_size + 2*(patch_size/2));
  cudaEventRecord(start);
  filter_kernel<<<grid_size,block_size,shared_mem_size*sizeof(float)>>>(dev_padded_image,dev_denoised_image,dev_g_kernel,patch_size,padded_size,f_sigma);
  cudaEventRecord(stop);
  /***************************************************************/

  checkCuda(cudaMemcpy(den_image,dev_denoised_image,image_size*image_size*sizeof(float),cudaMemcpyDeviceToHost));

  checkCuda(cudaEventSynchronize(stop));
  float milliseconds = 0;
  checkCuda(cudaEventElapsedTime(&milliseconds,start,stop));


  checkCuda(cudaEventDestroy(start));
  checkCuda(cudaEventDestroy(stop));
  checkCuda(cudaFree(dev_g_kernel));
  checkCuda(cudaFree(dev_padded_image));
  checkCuda(cudaFree(dev_denoised_image));

  //some print messages
  printf("*************************************************************************************\n");
  printf("gpu execution (use of shared memory) finished. Total time is %f msec\n",milliseconds);
  printf("image size is %d\n", image_size);
  printf("patch size is %d\n", patch_size);
  printf("filtes sigma is %f and gauss kernel sigma is %f\n",f_sigma,g_sigma);
  printf("**************************************************************************************\n");

  char denoised_path[200];
  snprintf(denoised_path,sizeof(denoised_path),"../data/gpu_sm_denoised_image_%d_%d.csv",image_size,patch_size);
  write_image2(denoised_path,den_image,image_size);

  free(GKERNEL);
  free(den_image);
  free(padded_I);

}



/***********host functions****************/

float **add_padding(float **image, int im_size, int p_size, int *pad_size){

  int padding     = (p_size-1)/2;
  int padded_size = im_size + padding*2;
  *pad_size = padded_size;

  //allocate memomry for padded image
  float **padded_image =(float **)malloc(padded_size*sizeof(float *));
  for(int i=0;i<padded_size;i++)padded_image[i]=(float *)malloc(padded_size*sizeof(float));

  //initialize all values to -1 for debugging puproses...Later this will not be necessary
  for(int i=0;i<padded_size;i++){
    for(int j=0;j<padded_size;j++){
      padded_image[i][j]=-1;
    }
  }

  //first copy the body of the padded array which will remain the same
  //and also compute outter pixels through mirroring.
  for(int row = padding;row<padding+im_size;row++){
    for(int col=0;col<padded_size;col++){
      if(col<padding)padded_image[row][col] = image[row-padding][padding-col-1];
      else if(col>=padding+im_size)padded_image[row][col] = image[row-padding][im_size-1-(col-padding-im_size)];
      else padded_image[row][col] = image[row-padding][col-padding];
    }
  }

  //mirroring first and last rows.
  for(int row=0;row<padding;row++){
    for(int col = 0;col<padded_size;col++)padded_image[row][col] = padded_image[padding+(padding-row-1)][col];
  }
  for(int row=padding+im_size;row<padded_size;row++){
    for(int col=0;col<padded_size;col++)padded_image[row][col] = padded_image[padding+im_size-1-(row-padding-im_size)][col];
  }


  return padded_image;
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


float *to_rowmajor(float **image, int im_size){

  //allocate memory
  float *row_image = (float *)malloc(im_size*im_size*sizeof(float));

  //copy pixels to 1-D array.
  for(int row=0;row<im_size;row++){
    for(int col=0;col<im_size;col++){
      row_image[row*im_size+col] = image[row][col];
    }
  }

  return row_image;
}


void checkCuda(cudaError_t result){
  if(result != cudaSuccess){
    printf("CUDA Runtime Error :%s\n",cudaGetErrorString(result));
    exit(1);
  }
}


float **read_image(char *im_path,int im_size){

  //open image file
  FILE *file = fopen(im_path,"r");

  //allocate memory
  float **image = (float **)malloc(im_size*sizeof(float *));
  for(int i=0;i<im_size;i++)image[i] = (float *)malloc(im_size*sizeof(float));

  //read the image
  for(int i=0;i<im_size;i++){
    for(int j=0;j<im_size;j++){
      fscanf(file,"%f,",image[i]+j);
    }
  }

  fclose(file);
  return image;
}


void   write_image2(char *d_path , float  *im, int im_size){

  //open file
  FILE *file = fopen(d_path,"w");

  //write image to file
  for(int i=0;i<im_size;i++){
    for(int j=0;j<im_size;j++){
      fprintf(file,"%f,",im[i*im_size +j]);
    }
    fprintf(file,"\n");
  }

  fclose(file);
}


/***********device functions****************/

__global__ void filter_kernel(float *image, float *denoised_image, float *gauss_kernel, int p_size, int im_size, float f_s ){

  //find the id for every thread
  int id = blockIdx.x*blockDim.x+threadIdx.x;

  //the shared memory must have all the pixels for every thread in the block
  //and also the padding in order to compute patches from here without the need
  //to look in global mem.
  extern __shared__ float image_section[];

  //save gauss kernel to shared memory, because we access it very often
  __shared__ float g_k[PSIZE*PSIZE];

  //for every thread calculate which pixel will compute. Taken into account the padding.
  int padded_pixel  = im_size*(p_size-1)/2+blockDim.x/(im_size-2*(p_size-1)/2)*blockIdx.x*im_size+
                      im_size*(threadIdx.x/(im_size-2*(p_size-1)/2))+threadIdx.x%(im_size-2*(p_size-1)/2)+(p_size-1)/2;


  //find for every pixel the corresponding patch
  float pixel_patch[PSIZE*PSIZE];
  find_patch(pixel_patch,padded_pixel,im_size,p_size,image);
  float  temp_patch[PSIZE*PSIZE];

  //write gauss kernel on shared mem.
  if(threadIdx.x == 0){
    for(int i=0;i<p_size;i++){
      for(int j=0;j<p_size;j++)
        g_k[i*p_size+j] = gauss_kernel[i*p_size+j];
    }
  }

  denoised_image[id] = -1;
  float Z =0;
  float W;
  float estimation=0;
  int temp_pixel;

  //load first elements on shared memory
  int shared_rows = blockDim.x/(im_size-2*(p_size-1)/2) + 2*(p_size-1)/2;
  int step = blockDim.x/(im_size-2*(p_size-1)/2);

  if(threadIdx.x<shared_rows){
    for(int i=0;i<im_size;i++)
      image_section[threadIdx.x*im_size+i]=image[threadIdx.x*im_size+i];
  }
  __syncthreads();

  for(int row=(p_size-1)/2;row<im_size-(p_size-1)/2;row = row+step){

    for(int s_row=(p_size-1)/2;s_row<shared_rows-(p_size-1)/2;s_row++){
      for(int s_col=(p_size-1)/1;s_col<im_size-(p_size-1)/2;s_col++){
        temp_pixel = s_row*im_size+s_col;
        find_patch(temp_patch,temp_pixel,im_size,p_size,image_section);
        W = compute_weight(pixel_patch,temp_patch,g_k,p_size,f_s);
        Z = Z + W;
        estimation = estimation + W*image_section[temp_pixel];
      }
    }

    //shift last elements
    __syncthreads();
    if(threadIdx.x == 0){
      for(int i=step;i<shared_rows;i++){
        for(int j=0;j<im_size;j++)
          image_section[im_size*(i-step)+j] = image_section[i*im_size+j];
      }
    }

    //bring the new elements in shared memory shifting the window steps time ahead
    __syncthreads();
    if(threadIdx.x>=shared_rows-step && threadIdx.x<shared_rows){
      for(int k=0;k<im_size;k++)
        image_section[threadIdx.x*im_size+k]=image[(row+step+threadIdx.x-(p_size-1)/2)*im_size+k];
    }
    __syncthreads();
  }

  denoised_image[id] = estimation/Z;
}


__device__ void find_patch(float *pixel_patch, int current_pixel,int im_size,int p_size,float *image){
  int boundary = (p_size-1)/2;
  for(int i=-boundary;i<=boundary;i++){
    for(int j=-boundary;j<=boundary;j++){
      pixel_patch[(i+boundary)*p_size+(j+boundary)] = image[current_pixel+i*im_size+j];
    }
  }
}


__device__ float compute_weight(float *pixel_patch,float *temp_patch,float *g_kernel, int p_size, float f_sigma){
  float temp=0;
  for(int i=0;i<p_size*p_size;i++) temp = temp + g_kernel[i]*g_kernel[i]*(pixel_patch[i]-temp_patch[i])*(pixel_patch[i]-temp_patch[i]);
  return exp(-temp/(f_sigma*f_sigma));
}
