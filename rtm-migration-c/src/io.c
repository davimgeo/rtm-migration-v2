#include <stdlib.h>

#include "io.h"

void write1d(const char* PATH, void* arr, size_t type, int size) 
{
  FILE* bin_data = fopen(PATH, "wb"); 
  if (bin_data == NULL) {
    printf("Could not write binary file.\n");
    exit(-1);
  }

  fwrite(arr, type, size, bin_data); 

  fclose(bin_data);   
}

void write2d(
    const char* PATH, void* arr, 
    size_t type, int height, int width
  ) 
{
  FILE* bin_data = fopen(PATH, "wb"); 
  if (bin_data == NULL) {
    printf("Could not write binary file.\n");
    exit(-1);
  }

  fwrite(arr, type, height*width, bin_data); 

  fclose(bin_data);   
}

float* read2d(const char* PATH, int row, int column) 
{
  float* arr = (float*)malloc(row*column*sizeof(float));

  FILE* bin_data = fopen(PATH, "rb"); 
  if (bin_data == NULL) 
  {
      printf("Could not read binary file.\n");
      exit(-1);
  }

  fread(arr, sizeof(float), row * column, bin_data); 

  fclose(bin_data);   

  return arr;
}

float* read2d_fortran(const char* PATH, int row, int column)
{
  float* arr = (float*)malloc(row*column*sizeof(float));

  FILE* bin_data = fopen(PATH, "rb"); 
  if (bin_data == NULL) 
  {
    printf("Could not read binary file at path: %s.\n",
        PATH);
    exit(-1);
  }

  fread(arr, sizeof(float), row * column, bin_data); 

  fclose(bin_data);   

  float* out = (float*)malloc(row*column*sizeof(float));

  for (int y = 0; y < row; y++) {
    for (int x = 0; x < column; x++) {

        out[y * column + x] = arr[x * row + y];
    }
  }

  return out;
}

