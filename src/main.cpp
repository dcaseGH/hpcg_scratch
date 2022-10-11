#include <mpi.h>
#include <fstream>
#include <iostream>
#include <cstdlib>
#ifdef HPCG_DETAILED_DEBUG
using std::cin;
#endif
using std::endl;

#include <vector>

#include "hpcg.hpp"

#include "CheckAspectRatio.hpp"
#include "GenerateGeometry.hpp"
#include "GenerateProblem.hpp"
#include "GenerateCoarseProblem.hpp"
#include "SetupHalo.hpp"
#include "CheckProblem.hpp"
#include "ExchangeHalo.hpp"
#include "OptimizeProblem.hpp"
#include "WriteProblem.hpp"
#include "ReportResults.hpp"
#include "mytimer.hpp"
#include "ComputeSPMV_ref.hpp"
#include "ComputeMG_ref.hpp"
#include "ComputeResidual.hpp"
#include "CG.hpp"
#include "CG_ref.hpp"
#include "CG_stencil.hpp"
#include "BiCG.hpp"
#include "BiCG_ref.hpp"
#include "Geometry.hpp"
#include "SparseMatrix.hpp"
#include "Vector.hpp"
#include "CGData.hpp"
#include "TestCG.hpp"
#include "TestSymmetry.hpp"
#include "TestNorms.hpp"
#include "string"
#include "ComputeSPMV.hpp"
#include "Compute27PtStencil.hpp"

//#include "optional"
//#include "hpcg.hpp"

/*
  Solve some small CG test case using matrices etc in hpcg form

*/

void set_vector(Vector &x)
{
   for (int i=0; i< x.localLength; i++)
   {
       x.values[i] = (double) i;
   }
}

//void print_vector(const Vector &x, optional<std::string>& title = std::nullopt)
void print_vector(const Vector &x, const std::string & title )
{
//  const double * const xv = x.values;
  //const double * const xv = x.values;
    std::cout<< "Vec " << title <<" " <<x.localLength << " " ;
//   for (auto q : &(x.values))
   for (int i=0; i< x.localLength; i++)
   {
       std::cout << x.values[i] << " ";
   }
   std::cout << std::endl;
}

void print_matrix(const SparseMatrix &A)
{
// Print the whole matrix (non zeros)

const local_int_t nrow = A.localNumberOfRows;

for (local_int_t i=0; i< nrow; i++)  {
  const double * const cur_vals = A.matrixValues[i];
//  const local_int_t * const cur_inds = A.mtxIndL[i];
  const int cur_nnz = A.nonzerosInRow[i];
  std::cout << "Mat ";
  for (int j=0; j< cur_nnz; j++)
  {
      std::cout << cur_vals[j] << " ";
  }
  std::cout << std::endl;
}
}

void copy_4x4_to_sparse_matrix(const double mat[4][4],  SparseMatrix &A)
{
// copy 4x4
const local_int_t nrow = A.localNumberOfRows;
int col_indx, row_indx;

col_indx = 0;
for (local_int_t i=0; i< nrow; i++)  {
  double * const cur_vals = A.matrixValues[i];
//  const local_int_t * const cur_inds = A.mtxIndL[i];
  const int cur_nnz = A.nonzerosInRow[i];
  std::cout << "Mat ";
  row_indx = 0;
  for (int j=0; j< cur_nnz; j++)
  {
      std::cout << cur_vals[j] << " ";
      cur_vals[j] = mat[col_indx][row_indx];
      row_indx++;
  }
  col_indx++;
  std::cout << std::endl;
}
}


int main(int argc, char * argv[]) {

  MPI_Init(&argc, &argv);

  HPCG_Params params;

  HPCG_Init(&argc, &argv, params);
// just enter some numbers manually for now (can use data file and command line if want later)
//  MPI_Comm_rank( MPI_COMM_WORLD, &params.comm_rank );
//  MPI_Comm_size( MPI_COMM_WORLD, &params.comm_size );


  int size = params.comm_size, rank = params.comm_rank; // Number of MPI processes, My process ID

  local_int_t nx,ny,nz;
  nx = (local_int_t)params.nx;
  ny = (local_int_t)params.ny;
  nz = (local_int_t)params.nz;
  int ierr = 0;  // Used to check return codes on function calls
  std::cout << "Starting the test run " << params.pz << " " << nz <<  std::endl;

//  ierr = CheckAspectRatio(0.125, nx, ny, nz, "local problem", rank==0);
//  if (ierr)
//    return ierr;

  /////////////////////////
  // Problem setup Phase //
  /////////////////////////

  // Construct the geometry and linear system
  Geometry * geom = new Geometry;
  GenerateGeometry(size, rank, params.numThreads, params.pz, params.zl, params.zu, nx, ny, nz, params.npx, params.npy, params.npz, geom);

  SparseMatrix A;
  InitializeSparseMatrix(A, geom);

  Vector b, x, xexact;//, xtemp;
  GenerateProblem(A, &b, &x, &xexact);
  SetupHalo(A);

  CGData data;
  InitializeSparseCGData(A, data);
/*
  //print matrix subroutine??
  print_vector(b,  "b");
  print_vector(x, "x");
  print_vector(xexact, "xexact");
  print_matrix(A);

  print_vector(data.Ap, "data AP"); */
  const int max_iter = 10;
  const double tolerance = 0.0;
  std::vector< double > times(9,0.0);
  double normr = 0.0;
  double normr0 = 0.0;
  int niters=0;

  ZeroVector(xexact); // Zero out x
  set_vector(xexact);
  print_vector(xexact, "xtemp");
  std::cout << " 27 pt stencil " << std::endl;
  Compute27PtStencil(A, xexact, x);
  print_vector(x, "x");

  ZeroVector(xexact); // Zero out x
  set_vector(xexact);
  std::cout << " SPMV " << std::endl;
  ComputeSPMV(A, xexact, x);
  print_vector(x, "x");

  ZeroVector(x); // Zero out x
  // Adding geom as arg for stencil
  ierr = CG_stencil( A, data, b, x, max_iter, tolerance, niters, normr,normr0, &times[0], true);//, geom);
//  ierr = CG( A, data, b, x, max_iter, tolerance, niters, normr,normr0, &times[0], true);
  if (ierr) HPCG_fout << "Error in call to CG: " << ierr << ".\n" << endl;
  if (rank==0) HPCG_fout << "Call [" << 0 << "] Scaled Residual [" << normr/normr0 << "]" << endl;

  print_vector(x, "x");

/*
  ZeroVector(x); // Zero out x
  ierr = BiCG( A, data, b, x, max_iter, tolerance, niters, normr,normr0, &times[0], false);
  if (ierr) HPCG_fout << "Error in call to CG: " << ierr << ".\n" << endl;
  if (rank==0) HPCG_fout << "Call [" << 0 << "] Scaled Residual [" << normr/normr0 << "]" << endl;

  print_vector(x, "x");

  ZeroVector(x); // Zero out x
  double mymat[4][4] = {26.08653784, -0.92259618, -0.95794096, -0.93858691,
       -0.92512332, 26.00558206, -0.91854604, -0.9929825 ,
       -0.93826258, -0.92952712, 26.04142306, -0.97612583,
       -0.99067017, -0.92771146, -0.93300909, 26.01593613};
  copy_4x4_to_sparse_matrix(mymat, A);
  print_matrix(A);
  ierr = BiCG( A, data, b, x, max_iter, tolerance, niters, normr,normr0, &times[0], false);
  if (ierr) HPCG_fout << "Error in call to CG: " << ierr << ".\n" << endl;
  if (rank==0) HPCG_fout << "Call [" << 0 << "] Scaled Residual [" << normr/normr0 << "]" << endl;

  print_vector(x, "x");
*/

  //HPCG_Finalize();
  MPI_Finalize();

  std::cout << "Completed the test run " << params.pz << " " << nz <<  std::endl;

  return 0 ;
}
