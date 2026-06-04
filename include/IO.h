#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <stdlib.h>

#define BUFFER_SIZE 256

static void IO_read_file_separed_by_comma(FILE* fptr, int total_size, float* result)
{
  char buff[BUFFER_SIZE];
  int len = 0, idx = 0;

  char ch;
  while ((ch = fgetc(fptr)) != EOF && idx < total_size) 
  {
    if (ch == '#') 
    {
      while (ch != '\n' && ch != EOF) 
        ch = fgetc(fptr);
      continue;
    }

    if (ch == ' ' || ch == '\n' || ch == ',') continue;

    len = 0;
    while (ch != ',' && ch != '\n' && 
        ch != EOF && ch != ' ') 
    {
        buff[len++] = ch;
        ch = fgetc(fptr);
    }
    buff[len] = '\0';

    result[idx++] = atof(buff);
  }
}

static inline void write1d(const char* PATH, void* arr, size_t type, int size) 
{
  FILE* bin_data = fopen(PATH, "wb"); 
  if (bin_data == NULL) {
    printf("Could not write binary file.\n");
    exit(-1);
  }

  fwrite(arr, type, size, bin_data); 

  fclose(bin_data);   
}

static inline void write2d(
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

static inline float* read2d(const char* PATH, int row, int column) 
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

static inline float* read2d_fortran(const char* PATH, int row, int column)
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

  float* out = (float*)malloc(row*column*sizeof(float));

  for (int y = 0; y < row; y++) {
    for (int x = 0; x < column; x++) {

        out[y * column + x] = arr[x * row + y];
    }
  }

  return out;

}

#endif
