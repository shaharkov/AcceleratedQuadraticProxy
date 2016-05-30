//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//% Code implementing the paper "Accelerated Quadratic Proxy for Geometric Optimization", SIGGRAPH 2016.
//% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
//%             Please contact the author to report any bugs.
//% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
//%            Meirav Galun (http://www.wisdom.weizmann.ac.il/~/meirav/)
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "mex.h"
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "mexHelpers.cpp"

using namespace Eigen;

void projBlockRotation(VectorXd &pA, int dim)
{
	int block_size = dim*dim;
	int num_blocks = pA.size() / block_size;
	JacobiSVD<MatrixXd> svdA(dim, dim, (ComputeFullU | ComputeFullV));
	bool flipped;
	MatrixXd U, V;
	Map<MatrixXd> currA(pA.data(), dim, dim);

	// project
	for (int ii = 0; ii < num_blocks; ii++)
	{
		// get current block
		new (&currA) Map<MatrixXd>(pA.data() + ii*block_size, dim, dim);
		// sign of determinant 
		flipped = (currA.determinant() < 0);
		// svd
		svdA.compute(currA);
		// compute frames
		U = svdA.matrixU();
		V = svdA.matrixV();
		// ssvd
		if (flipped)
		{
			U.col(dim - 1) = -U.col(dim - 1);
		}
		// project block
		currA = U*V.transpose();
	}
}

void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])

{
	// assign input
	int A_rows = mxGetM(prhs[0]); // # rows of A
	int A_cols = mxGetN(prhs[0]); // # cols of A
	double *dim;
	const Map<VectorXd> A(mxGetPr(prhs[0]), A_rows, A_cols);
	dim = mxGetPr(prhs[1]);
	
	// init output
	VectorXd pA(A_rows);
	pA = A;

	// project
	projBlockRotation(pA, *dim);
	
	// assign outputs
	mapDenseMatrixToMex(pA, &(plhs[0]));

}