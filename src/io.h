#ifndef IO_H
#define IO_H

#include <stdio.h>

#define BUFFER_SIZE 256

void write1d(const char* PATH, void* arr, size_t type, int size);

void write2d(
  const char* PATH, void* arr, 
  size_t type, int height, int width
);

float* read2d(const char* PATH, int row, int column);

float* read2d_fortran(const char* PATH, int row, int column);

#endif

