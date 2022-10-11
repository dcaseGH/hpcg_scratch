
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
 @file CGData.hpp

 HPCG data structure
 */

#ifndef CGDATA_HPP
#define CGDATA_HPP

#include "SparseMatrix.hpp"
#include "Vector.hpp"

struct CGData_STRUCT {
  Vector r; //!< pointer to residual vector
  Vector r2; //!< pointer to residual vector
  Vector z; //!< pointer to preconditioned residual vector
  Vector z2; //!< pointer to preconditioned residual vector
  Vector p; //!< pointer to direction vector
  Vector Ap; //!< pointer to Krylov vector
  Vector x2; //!< pointer to spare x for biCG
  Vector p2A; //!< pointer to transpose vector for biCG
  Vector p2; //!< pointer to transpose vector for biCG
};
typedef struct CGData_STRUCT CGData;

/*!
 Constructor for the data structure of CG vectors.

 @param[in]  A    the data structure that describes the problem matrix and its structure
 @param[out] data the data structure for CG vectors that will be allocated to get it ready for use in CG iterations
 */
inline void InitializeSparseCGData(SparseMatrix & A, CGData & data) {
  local_int_t nrow = A.localNumberOfRows;
  local_int_t ncol = A.localNumberOfColumns;
  InitializeVector(data.r, nrow);
  InitializeVector(data.r2, ncol);
  InitializeVector(data.z, ncol);
  InitializeVector(data.z2, nrow);
  InitializeVector(data.p, ncol);
  InitializeVector(data.Ap, nrow);
  InitializeVector(data.x2, nrow); //check this stuff
  InitializeVector(data.p2, nrow); //check this stuff
  InitializeVector(data.p2A, ncol);
  return;
}

/*!
 Destructor for the CG vectors data.

 @param[inout] data the CG vectors data structure whose storage is deallocated
 */
inline void DeleteCGData(CGData & data) {

  DeleteVector (data.r);
  DeleteVector (data.r2);
  DeleteVector (data.z);
  DeleteVector (data.z2);
  DeleteVector (data.p);
  DeleteVector (data.Ap);
  DeleteVector (data.x2);
  DeleteVector (data.p2);
  DeleteVector (data.p2A);
  return;
}

#endif // CGDATA_HPP

