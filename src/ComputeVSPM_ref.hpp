
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

#ifndef COMPUTEVSPM_REF_HPP
#define COMPUTEVSPM_REF_HPP
#include "Vector.hpp"
#include "SparseMatrix.hpp"

int ComputeVSPM_ref( Vector  & x, const SparseMatrix & A, Vector & y);

#endif  // COMPUTEVSPM_REF_HPP
