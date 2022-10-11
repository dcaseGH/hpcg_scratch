
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

#ifndef COMPUTE27PtStencil_REF_HPP
#define COMPUTE27PtStencil_REF_HPP
#include "Vector.hpp"
#include "SparseMatrix.hpp"

int Compute27PtStencil_ref( const SparseMatrix & A, Vector  & x, Vector & y);

#endif  // COMPUTESPMV_REF_HPP
