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

#endif
