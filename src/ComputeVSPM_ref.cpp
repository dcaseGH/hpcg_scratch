
//@HEADER
// ***************************************************
//
// HPCG: High Performance Conjugate Gradient Benchmark
//
// Contact:
// Michael A. Heroux ( maherou@sandia.gov)
// Jack Dongarra     (dongarra@eecs.utk.edu)
// Piotr Luszczek    (luszczek@eecs.utk.edu)
//
// ***************************************************
//@HEADER

/*!
 @file ComputeVSPM_ref.cpp

 HPCG routine

 NOTE - is it possible to do this fast? We don't need values for matrix, so perhaps not an issue??
 */
#include <iostream>

#include "ComputeVSPM_ref.hpp"

#ifndef HPCG_NO_MPI
#include "ExchangeHalo.hpp"
#endif

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif
#include <cassert>

/*!
  Routine to compute matrix vector product y = xA where:
  Precondition: First call exchange_externals to get off-processor values of x

  This is the reference VSPM implementation.  It CANNOT be modified for the
  purposes of this benchmark.

  @param[in]  A the known system matrix
  @param[in]  x the known vector
  @param[out] y the On exit contains the result: xA.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeVSPM
*/
int ComputeVSPM_ref( Vector & x, const SparseMatrix & A, Vector & y) {

  assert(y.localLength>=A.localNumberOfColumns); // Test vector lengths
  assert(x.localLength>=A.localNumberOfRows);

#ifndef HPCG_NO_MPI
    ExchangeHalo(A,x);
    //will this work for non square?? or transpose?? dhc?
#endif
  const double * const xv = x.values;
  double * const yv = y.values;
  const local_int_t nrow = A.localNumberOfRows;
  const local_int_t ncol = A.localNumberOfColumns;
#ifndef HPCG_NO_OPENMP
  #pragma omp parallel for
#endif

  for (local_int_t i=0; i< nrow; i++)  {
//  for (local_int_t i=0; i< ncol; i++)  {
//    double sum = 0.0;
    yv[i] = 0.0;
    const double * const cur_vals = A.matrixValues[i];
    const local_int_t * const cur_inds = A.mtxIndL[i];
    const int cur_nnz = A.nonzerosInRow[i];

    for (int j=0; j< cur_nnz; j++){
        yv[j] = xv[cur_inds[j]] * cur_vals[j];
//      sum += cur_vals[j]*xv[cur_inds[j]];
//    yv[i] = sum;
    }
  }
  return 0;
}
