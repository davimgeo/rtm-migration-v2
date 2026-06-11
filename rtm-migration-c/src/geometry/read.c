#include <stdio.h>
#include <stdlib.h>

#include "read.h"

#define DATA_COL 3

static void read_file_separed_by_comma(FILE* fptr, int total_size, float* result)
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

static void read_receivers(geometry_t* geom, const char* PATH)
{
  FILE* fptr = fopen(PATH, "r");
  if (!fptr) err("Could not open %s\n", PATH);

  char* line = NULL;
  size_t line_buffer_len = 0;
  int f_lines = 0;

  while (getline(&line, &line_buffer_len, fptr) != EOF)
  {
    if (line[0] != '#') f_lines++;
  }

  long size = ftell(fptr); rewind(fptr);

  float* result = allocf(size);

  read_file_separed_by_comma(fptr, size, result);

  geom->rec = alloc_struct(f_lines, geom->rec);

  for (int i = 0; i < f_lines; i++)
  {
    geom->rec->x[i] = result[i * DATA_COL + 1];
    geom->rec->z[i] = result[i * DATA_COL + 2];
  }

  free(result);
  free(line);
  fclose(fptr);
}

static void read_sources(geometry_t* geom, const char* PATH)
{
  FILE* fptr = fopen(PATH, "r");
  if (!fptr) err("Could not open %s\n", PATH);

  char* line = NULL;
  size_t line_buffer_len = 0;
  int f_lines = 0;

  while (getline(&line, &line_buffer_len, fptr) != EOF)
  {
    if (line[0] != '#') f_lines++;
  }

  long size = ftell(fptr); rewind(fptr);

  float* result = allocf(size);

  read_file_separed_by_comma(fptr, size, result);

  for (int i = 0; i < f_lines; i++)
  {
    geom->src->x[i] = result[i * DATA_COL + 1];
    geom->src->z[i] = result[i * DATA_COL + 2];
  }

  free(result);
  free(line);
  fclose(fptr);
}

void geometry_load(
  geometry_t* geom, 
  const char* REC_PATH, 
  const char* SRC_PATH
)
{
  read_receivers(geom, REC_PATH);
  read_sources(geom, SRC_PATH);
}
