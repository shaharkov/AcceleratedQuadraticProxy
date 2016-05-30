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
#include <svd3x3.h>

using namespace Eigen;
using namespace igl;

void projBlockRotation2x2(VectorXd &pA, int dim)
{
	int block_size = dim*dim;
	int num_blocks = pA.size() / block_size;
	Map<MatrixXd> currA(pA.data(), dim, dim);
	Vector2d b;

	// project
	for (int ii = 0; ii < num_blocks; ii++)
	{
		// get current block
		new (&currA) Map<MatrixXd>(pA.data() + ii*block_size, dim, dim);
		// closest similarity
		b << 0.5*(currA(0, 0) + currA(1, 1)), 0.5*(currA(0, 1) - currA(1, 0)); // first row of B
		// closest rotation
		b = b / b.norm();
		currA << b(0), b(1), -b(1), b(0);
	}
}

void projBlockRotation3x3(VectorXd &pA, int dim)
{
	int block_size = dim*dim;
	int num_blocks = pA.size() / block_size;

	Matrix3f currAf, R;
	Matrix3f U, V;
	Vector3f s;
	Map<Matrix3d> currA(pA.data());

	// project
	for (int ii = 0; ii < num_blocks; ii++)
	{
		// get current block
		new (&currA) Map<MatrixXd>(pA.data() + ii*block_size, dim, dim);
		//
		currAf = currA.cast<float>(); // double -> single
		svd3x3(currAf, U, s, V); // svd
		R = U * V.transpose(); // polar
		currA = R.cast<double>();
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
	if (*dim == 2)
		projBlockRotation2x2(pA, *dim);	
	else if (*dim == 3)
		projBlockRotation3x3(pA, *dim);
	else
		mexErrMsgIdAndTxt("MATLAB:wrong_dimension", "dim must be either 2 or 3");
	
	// assign outputs
	mapDenseMatrixToMex(pA, &(plhs[0]));

}