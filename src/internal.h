#ifndef INTERNAL_H
#define INTERNAL_H

// make a flag inside ifndef to check
// if user has GNUPLOT
// flag decided inside makefile
// put macros here

#include "io.h"
#include "utils.h"
// do a check if user has GNUPLOT
#include "plot.h"

#define debug(...) \
    debug_impl(__FILE__, __func__, __LINE__, __VA_ARGS__)

#define err(...) \
    err_impl(__FILE__, __func__, __LINE__, __VA_ARGS__)

#define alloc_struct(count, ptr) \
    ((__typeof__(ptr)) allocs_impl((count), sizeof(*(ptr))))

#endif
