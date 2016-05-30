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

using namespace Eigen;

void mapSparseMatrixToMex(const SparseMatrix<double>& mat, mxArray **out)
{
	*out = mxCreateSparse(mat.rows(), mat.cols(), mat.nonZeros(), mxREAL);

	mwIndex *Ir, *Jc;
	double *Pr;
	Ir = mxGetIr(*out);
	Jc = mxGetJc(*out);
	Pr = mxGetPr(*out);

	for (int k = 0; k < mat.nonZeros(); k++)
	{
		Pr[k] = (mat.valuePtr())[k];
		Ir[k] = (mat.innerIndexPtr())[k];
	}
	for (int k = 0; k <= mat.cols(); k++)
	{
		Jc[k] = (mat.outerIndexPtr())[k];
	}
}

void mapDenseMatrixToMex(const MatrixXd& mat, mxArray **out)
{
	*out = mxCreateDoubleMatrix(mat.rows(), mat.cols(), mxREAL);
	Map<MatrixXd> temp(mxGetPr(*out), mat.rows(), mat.cols());
	temp = mat; // copy
}