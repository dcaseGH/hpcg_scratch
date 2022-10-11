
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
 @file ComputeSPMV_ref.cpp

 HPCG routine
 */
#include <iostream>

#include "Compute9PtHorizontalStencil_ref.hpp"

#ifndef HPCG_NO_MPI
#include "ExchangeHalo.hpp"
#endif

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif
#include <cassert>

/*!
  Routine to compute matrix vector product y = Ax where:
  Precondition: First call exchange_externals to get off-processor values of x

  This is the reference SPMV implementation.  It CANNOT be modified for the
  purposes of this benchmark.

  @param[in]  A the known system matrix
  @param[in]  x the known vector
  @param[out] y the On exit contains the result: Ax.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeSPMV
*/



int Compute9PtHorizontalStencil_ref( const SparseMatrix & A, Vector & x, Vector & y) {

  assert(x.localLength>=A.localNumberOfColumns); // Test vector lengths
  assert(y.localLength>=A.localNumberOfRows);

#ifndef HPCG_NO_MPI
    ExchangeHalo(A,x);
#endif
  std::cout << "Doing stencil " << x.localLength << " " <<x.localLength << " " << y.localLength << std::endl;
  std::cout << "Local implementation at first (no mpi) " << std::endl;

  const double * const xv = x.values;
  double * const yv = y.values;

  local_int_t ix = 0;
  local_int_t iy = 0;
  local_int_t nx = A.geom->nx;
  local_int_t ny = A.geom->ny;

  for (ix=1; ix < nx-1; ix++){
      for( iy=1; iy< ny-1; iy++){
            yv[ix+iy*nx] = 26. * xv[ix+iy*nx]
                                 - xv[ix-1+iy*nx]
                                 - xv[ix+1+iy*nx]
                                 - xv[ix+(iy-1)*nx]
                                 - xv[ix+(iy+1)*nx]
                                 - xv[ix-1+(iy-1)*nx]
                                 - xv[ix-1+(iy+1)*nx]
                                 - xv[ix+1+(iy-1)*nx]
                                 - xv[ix+1+(iy+1)*nx] ;
      }
   }
   ix = 0;
   for( iy=1; iy< ny-1; iy++){
          yv[ix+iy*nx] = 26. * xv[ix+iy*nx]
                             - xv[ix+1+iy*nx]
                                 - xv[ix+(iy-1)*nx]
                                 - xv[ix+(iy+1)*nx]
                                 - xv[ix+1+(iy-1)*nx]
                                 - xv[ix+1+(iy+1)*nx] ;
   }
   ix = nx;
   for( iy=1; iy< ny-1; iy++){
       yv[ix+iy*nx] = 26. * xv[ix+iy*nx]
                                 - xv[ix-1+iy*nx]
                                 - xv[ix+(iy-1)*nx]
                                 - xv[ix+(iy+1)*nx]
                                 - xv[ix-1+(iy-1)*nx]
                                 - xv[ix-1+(iy+1)*nx] ;
    }
    iy = 0;
    for (ix=1; ix < nx-1; ix++){
            yv[ix+iy*nx] = 26. * xv[ix+iy*nx]
                                 - xv[ix-1+iy*nx]
                                 - xv[ix+1+iy*nx]
                                 - xv[ix+(iy+1)*nx]
                                 - xv[ix-1+(iy+1)*nx]
                                 - xv[ix+1+(iy+1)*nx] ;
    }
    iy = ny;
    for (ix=1; ix < nx-1; ix++){
        yv[ix+iy*nx] = 26. * xv[ix+iy*nx]
                                 - xv[ix-1+iy*nx]
                                 - xv[ix+1+iy*nx]
                                 - xv[ix+(iy-1)*nx]
                                 - xv[ix-1+(iy-1)*nx]
                                 - xv[ix+1+(iy-1)*nx] ;
    }
    ix, iy = 0,0;
    yv[ix+iy*nx] = 26. * xv[ix+iy*nx]
                         - xv[ix+1+iy*nx]
                         - xv[ix+(iy+1)*nx]
                         - xv[ix+1+(iy+1)*nx] ;
    ix, iy = nx-1,ny-1;
    yv[ix+iy*nx] = 26. * xv[ix+iy*nx]
                         - xv[ix-1+iy*nx]
                         - xv[ix+(iy-1)*nx]
                         - xv[ix-1+(iy-1)*nx] ;
    ix, iy = nx-1,0;
    yv[ix+iy*nx] = 26. * xv[ix+iy*nx]
                         - xv[ix-1+iy*nx]
                         - xv[ix+(iy+1)*nx]
                         - xv[ix-1+(iy+1)*nx] ;
    ix, iy = 0,ny-1;
    yv[ix+iy*nx] = 26. * xv[ix+iy*nx]
                         - xv[ix+1+iy*nx]
                         - xv[ix+(iy-1)*nx]
                         - xv[ix+1+(iy-1)*nx];

  return 0;
}
