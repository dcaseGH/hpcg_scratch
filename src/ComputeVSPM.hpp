
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

#ifndef COMPUTEVSPM_HPP
#define COMPUTEVSPM_HPP
#include "Vector.hpp"
#include "SparseMatrix.hpp"

int ComputeVSPM( Vector & x, const SparseMatrix & A, Vector & y);

#endif  // COMPUTEVSPM_HPP
