
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

#ifndef COMPUTE9PTHS_REF_HPP
#define COMPUTE9PTHS_REF_HPP
#include "Vector.hpp"
#include "SparseMatrix.hpp"

int Compute9PtHorizontalStencil_ref( const SparseMatrix & A, Vector  & x, Vector & y);

#endif  // COMPUTESPMV_REF_HPP
