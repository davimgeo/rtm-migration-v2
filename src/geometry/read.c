#include "../internal.h"

#include "geometry.h"

#define BUFFER_SIZE 2048
#define DATA_COL 3

geometry_t* Geometry_InitLoad(geometry_t* g)
{
  g = alloc_struct(1.0f, g);

  g->nsrc = 0;
  g->nrec = 0;

  g->rec.x = allocf(BUFFER_SIZE);
  g->rec.z = allocf(BUFFER_SIZE);

  g->src.x = allocf(BUFFER_SIZE);
  g->src.z = allocf(BUFFER_SIZE);

  return g;
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


  for (int i = 0; i < geom->nrec; i++)
  {
    geom->rec.x[i] = result[i * DATA_COL + 1] / 10.0f;
    geom->rec.z[i] = result[i * DATA_COL + 2] / 10.0f;
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

  for (int i = 0; i < geom->nsrc; i++)
  {
    geom->src.x[i] = result[i * DATA_COL + 1] / 10.0f;
    geom->src.z[i] = result[i * DATA_COL + 2] / 10.0f;
  }

  free(result);
  free(line);
  fclose(fptr);
}

void Geometry_Load(
  geometry_t* geom, 
  const char* REC_PATH, 
  const char* SRC_PATH,
  unsigned int flags
)
{
  if(flags & GEOMETRY_ONLY_RECEIVERS)
  {
    read_receivers(geom, REC_PATH);
  }
  else
  {
    read_receivers(geom, REC_PATH);
    read_sources(geom, SRC_PATH);
  }
}
