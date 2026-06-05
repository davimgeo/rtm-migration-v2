#include <cstdlib>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "utils.h"

#include "../include/IO.h"
#include "../include/geometry.hpp"

#define DATA_COL 3

void Geometry::get()
{
  char* mode = upper(c.geometry_mode);
  if (strcmp(mode, "LOAD") == 0) 
  {
    read_receivers();
    read_sources();
  }
  else if (strcmp(mode, "CREATE") == 0) 
  {
    create_receivers();
    create_sources();

    if(c.SAVE) save();
  }
  else 
  {
    perror("Choose a valid mode. (create, load)");
    exit(-1);
  }
};

void Geometry::create_receivers()
{
  rec = (receiver*)malloc(nrec * sizeof(receiver));
  if (rec == NULL) {
    perror("Could not allocate\n");
    exit(EXIT_FAILURE);
  }

  for (size_t i = 0; i < nrec; i++) 
  {
    rec[i].x = i * c.offset_rec;
    rec[i].z = c.rec_depth_create;
  }
}

void Geometry::create_sources()
{
  src = (sources*)malloc(nsrc * sizeof(sources));
  if (src == NULL) {
    perror("Could not allocate\n");
    exit(EXIT_FAILURE);
  }

  for (size_t i = 0; i < nsrc; i++) 
  {
    src[i].x = i * c.offset_src;
    src[i].z = c.src_depth_create;
  }
}

void Geometry::save()
{
  int BUFFER = 64;

  // save receivers
  FILE *fptr = fopen(c.SAVE_REC_PATH, "wa");

  fprintf(fptr, "# recId, recx, recz\n");

  char line[BUFFER];  
  for (int i = 0; i < nrec; i++) 
  {
    snprintf(
      line,
      BUFFER,
      "%d, %.2f, %.2f\n",
      i, rec[i].x, rec[i].z
    );

    fprintf(fptr, line, BUFFER);
  }

  fclose(fptr);

  printf("test\n");
  // save sources
  FILE *fptr2 = fopen(c.SAVE_SRC_PATH, "wa");

  fprintf(fptr2, "# srcId, srcx, srcz\n");

  char line2[BUFFER];  
  for (int i = 0; i < nsrc; i++) 
  {
    snprintf(
      line2,
      BUFFER,
      "%d, %.2f, %.2f\n",
      i, src[i].x, src[i].z
    );

    fprintf(fptr2, line2, BUFFER);
  }

  fclose(fptr2);
}

void Geometry::read_receivers()
{
  FILE* fptr = fopen(c.receivers_path, "r");
  if (!fptr)
  {
    perror(c.receivers_path);
    exit(EXIT_FAILURE);
  }

  char* line = NULL;
  size_t line_buffer_len = 0;
  int f_lines = 0;

  while (getline(&line, &line_buffer_len, fptr) != EOF)
  {
    if (line[0] != '#') f_lines++;
  }

  long size = ftell(fptr); rewind(fptr);

  float* result = (float*)malloc(size * sizeof(float));
  if (result == NULL) {
    perror("Could not allocate\n");
    exit(EXIT_FAILURE);
  }

  IO_read_file_separed_by_comma(fptr, size, result);

  rec = (receiver*)malloc(f_lines * sizeof(receiver));
  if (rec == NULL) {
    perror("Could not allocate\n");
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < f_lines; i++)
  {
    rec[i].x = result[i * DATA_COL + 1];
    rec[i].z = result[i * DATA_COL + 2];
  }

  free(result);
  free(line);
  fclose(fptr);
}

void Geometry::read_sources()
{
  FILE* fptr = fopen(c.sources_path, "r");
  if (!fptr)
  {
    perror(c.sources_path);
    exit(EXIT_FAILURE);
  }

  char* line = NULL;
  size_t line_buffer_len = 0;
  int f_lines = 0;

  while (getline(&line, &line_buffer_len, fptr) != EOF)
  {
    if (line[0] != '#') f_lines++;
  }

  long size = ftell(fptr); rewind(fptr);

  float* result = (float*)malloc(size * sizeof(float));
  if (result == NULL) {
    perror("Could not allocate\n");
    exit(EXIT_FAILURE);
  }

  IO_read_file_separed_by_comma(fptr, size, result);

  src = (sources*)malloc(f_lines * sizeof(sources));
  if (src == NULL) {
    perror("Could not allocate\n");
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < f_lines; i++)
  {
    src[i].x = result[i * DATA_COL + 1];
    src[i].z = result[i * DATA_COL + 2];
  }

  free(result);
  free(line);
  fclose(fptr);
}
