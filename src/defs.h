#ifndef DEFS_H
#define DEFS_H

#include "omp.h"

namespace algo{
  enum algo{MINRES = 1};
}

namespace method{
  enum method{SEQUENTAL = 1, OMP = 2, OPENCL = 3, CILKPLUS = 4, CUDA = 5};
}

/**
 * Computes number of active threads in current moment
 */
int omp_thread_count() {
  int n = 0;
#pragma omp parallel reduction(+:n)
  n += 1;
  return n;
}

#endif