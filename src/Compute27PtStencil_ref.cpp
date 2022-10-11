
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

#include "Compute27PtStencil_ref.hpp"

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


//add openmp later
int Compute27PtStencil_ref( const SparseMatrix & A, Vector & x, Vector & y) {

  assert(x.localLength>=A.localNumberOfColumns); // Test vector lengths
  assert(y.localLength>=A.localNumberOfRows);

#ifndef HPCG_NO_MPI
    ExchangeHalo(A,x);
#endif
  std::cout << "Doing 27pt stencil " << x.localLength << " " <<x.localLength << " " << y.localLength << std::endl;

  const double * const xv = x.values;
  double * const yv = y.values;

  local_int_t ix = 0;
  local_int_t iy = 0;
  local_int_t iz = 0;
  local_int_t nx = A.geom->nx;
  local_int_t ny = A.geom->ny;
  local_int_t nz = A.geom->nz;

  std::cout << "Local implementation at first (no mpi) "<< nx << " " << ny << " " << nz << " " << std::endl;

//bulk part
  for( iz=1; iz< nz-1; iz++){
    for( iy=1; iy< ny-1; iy++){
        for (ix=1; ix < nx-1; ix++){
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
- xv[ix+-1 + (iy+-1)*nx + (iz+-1)*nx*ny]
- xv[ix+-1 + (iy+-1)*nx + (iz+0)*nx*ny]
- xv[ix+-1 + (iy+-1)*nx + (iz+1)*nx*ny]
- xv[ix+-1 + (iy+0)*nx + (iz+-1)*nx*ny]
- xv[ix+-1 + (iy+0)*nx + (iz+0)*nx*ny]
- xv[ix+-1 + (iy+0)*nx + (iz+1)*nx*ny]
- xv[ix+-1 + (iy+1)*nx + (iz+-1)*nx*ny]
- xv[ix+-1 + (iy+1)*nx + (iz+0)*nx*ny]
- xv[ix+-1 + (iy+1)*nx + (iz+1)*nx*ny]
- xv[ix+0 + (iy+-1)*nx + (iz+-1)*nx*ny]
- xv[ix+0 + (iy+-1)*nx + (iz+0)*nx*ny]
- xv[ix+0 + (iy+-1)*nx + (iz+1)*nx*ny]
- xv[ix+0 + (iy+0)*nx + (iz+-1)*nx*ny]
- xv[ix+0 + (iy+0)*nx + (iz+1)*nx*ny]
- xv[ix+0 + (iy+1)*nx + (iz+-1)*nx*ny]
- xv[ix+0 + (iy+1)*nx + (iz+0)*nx*ny]
- xv[ix+0 + (iy+1)*nx + (iz+1)*nx*ny]
- xv[ix+1 + (iy+-1)*nx + (iz+-1)*nx*ny]
- xv[ix+1 + (iy+-1)*nx + (iz+0)*nx*ny]
- xv[ix+1 + (iy+-1)*nx + (iz+1)*nx*ny]
- xv[ix+1 + (iy+0)*nx + (iz+-1)*nx*ny]
- xv[ix+1 + (iy+0)*nx + (iz+0)*nx*ny]
- xv[ix+1 + (iy+0)*nx + (iz+1)*nx*ny]
- xv[ix+1 + (iy+1)*nx + (iz+-1)*nx*ny]
- xv[ix+1 + (iy+1)*nx + (iz+0)*nx*ny]
- xv[ix+1 + (iy+1)*nx + (iz+1)*nx*ny]
;}
}
}

// sides
iz = 0;
for (ix=1; ix<nx-1;ix++){
    for (iy=1; iy<ny-1;iy++)
{
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+1)*ny*nx]
;}
}
iz = nz-1;
for (ix=1; ix<nx-1;ix++){
    for (iy=1; iy<ny-1;iy++)
{
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+0)*ny*nx]
;}
}
iy = 0;
for (ix=1; ix<nx-1;ix++){
    for (iz=1; iz<nz-1;iz++)
{
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+1)*ny*nx]
;}
}
iy = ny-1;
for (ix=1; ix<nx-1;ix++){
    for (iz=1; iz<nz-1;iz++)
{
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+1)*ny*nx]
;}
}
ix = 0;
for (iy=1; iy<ny-1;iy++){
    for (iz=1; iz<nz-1;iz++)
{
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+1)*ny*nx]
;}
}
ix = nx-1;
for (iy=1; iy<ny-1;iy++){
    for (iz=1; iz<nz-1;iz++)
{
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+1)*ny*nx]
;}
}
// edges

ix = 0;
iy = 0;
for (iz=1; iz<nz-1;iz++)
{
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+1)*ny*nx]
;}
ix = 0;
iy = ny-1;
for (iz=1; iz<nz-1;iz++)
{
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+1)*ny*nx]
;}
ix = nx-1;
iy = 0;
for (iz=1; iz<nz-1;iz++)
{
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+1)*ny*nx]
;}
ix = nx-1;
iy = ny-1;
for (iz=1; iz<nz-1;iz++)
{
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+1)*ny*nx]
;}
ix = 0;
iz = 0;
for (iy=1; iy<ny-1;iy++)
{
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+1)*ny*nx]
;}
ix = 0;
iz = nz-1;
for (iy=1; iy<ny-1;iy++)
{
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+0)*ny*nx]
;}
ix = nx-1;
iz = 0;
for (iy=1; iy<ny-1;iy++)
{
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+1)*ny*nx]
;}
ix = nx-1;
iz = nz-1;
for (iy=1; iy<ny-1;iy++)
{
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+0)*ny*nx]
;}
iy = 0;
iz = 0;
for (ix=1; ix<nx-1;ix++)
{
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+1)*ny*nx]
;}
iy = 0;
iz = nz-1;
for (ix=1; ix<nx-1;ix++)
{
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+0)*ny*nx]
;}
iy = ny-1;
iz = 0;
for (ix=1; ix<nx-1;ix++)
{
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+1)*ny*nx]
;}
iy = ny-1;
iz = nz-1;
for (ix=1; ix<nx-1;ix++)
{
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+0)*ny*nx]
;}

// Corners

ix = 0;
iy = 0;
iz = 0;
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+1)*ny*nx]
;
ix = 0;
iy = 0;
iz = nz-1;
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+1)*nx+(iz+0)*ny*nx]
;
ix = 0;
iy = ny-1;
iz = 0;
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+1)*ny*nx]
;
ix = 0;
iy = ny-1;
iz = nz-1;
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+1)+(iy+0)*nx+(iz+0)*ny*nx]
;
ix = nx-1;
iy = 0;
iz = 0;
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+1)*ny*nx]
;
ix = nx-1;
iy = 0;
iz = nz-1;
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+1)*nx+(iz+0)*ny*nx]
;
ix = nx-1;
iy = ny-1;
iz = 0;
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+1)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+1)*ny*nx]
;
ix = nx-1;
iy = ny-1;
iz = nz-1;
yv[ix+iy*nx+iz*ny*nx] = 26. * xv[ix+iy*nx+iz*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+-1)*ny*nx]
-xv[(ix+-1)+(iy+0)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+-1)*ny*nx]
-xv[(ix+0)+(iy+-1)*nx+(iz+0)*ny*nx]
-xv[(ix+0)+(iy+0)*nx+(iz+-1)*ny*nx]
;

  return 0;
}

