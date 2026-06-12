#include "../internal.h"
#include "geometry.h"

#define BUFFER_SIZE 256
#define DATA_COL 3

void Geometry_InitLoad(geometry_t* g)
{
  g = alloc_struct(1.0f, g);

  g->nsrc = 0;
  g->nrec = 0;
}

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
  geom->nrec = 0;

  while (getline(&line, &line_buffer_len, fptr) != EOF)
  {
    if (line[0] != '#') geom->nrec++;
  }

  long size = ftell(fptr); rewind(fptr);

  float* result = allocf(size);

  read_file_separed_by_comma(fptr, size, result);

  geom->rec.x = allocf(geom->nrec);
  geom->rec.z = allocf(geom->nrec);

  for (int i = 0; i < geom->nrec; i++)
  {
    geom->rec.x[i] = result[i * DATA_COL + 1];
    geom->rec.z[i] = result[i * DATA_COL + 2];
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
  geom->nsrc = 0;

  while (getline(&line, &line_buffer_len, fptr) != EOF)
  {
    if (line[0] != '#') geom->nsrc++;
  }

  long size = ftell(fptr); rewind(fptr);

  float* result = allocf(size);

  read_file_separed_by_comma(fptr, size, result);

  geom->src.x = allocf(geom->nsrc);
  geom->src.z = allocf(geom->nsrc);

  for (int i = 0; i < geom->nsrc; i++)
  {
    geom->src.x[i] = result[i * DATA_COL + 1];
    geom->src.z[i] = result[i * DATA_COL + 2];
  }

  free(result);
  free(line);
  fclose(fptr);
}

void Geometry_Load(
  geometry_t* geom, 
  const char* REC_PATH, 
  const char* SRC_PATH
)
{
  read_receivers(geom, REC_PATH);
  read_sources(geom, SRC_PATH);
}
