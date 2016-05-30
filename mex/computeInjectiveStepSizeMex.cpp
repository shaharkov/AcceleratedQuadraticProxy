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
#include <complex>
#include "mexHelpers.cpp"

using namespace Eigen;
using namespace std;

double getSmallestPositiveRealQuadRoot(double a, double b, double c, double tol)
{
	// return negative value if no positive real root is found
	double t;

	if (abs(a) <= tol)
		t = -c / b;
	else
	{
		double desc = b*b - 4 * a*c;
		if (desc > 0)
		{
			t = (-b - sqrt(desc)) / (2 * a);
			if (t < 0)
				t = (-b + sqrt(desc)) / (2 * a);
		}
		else // desv<0 ==> imag
			t = -1;
	}
	return t;
}


void computeInjectiveStepSize_2d(const MatrixXd& F, const MatrixXd& x, const MatrixXd& p, double tol, double& t_max)
{
	int n_tri = F.rows();
	double x1, x2, x3, y1, y2, y3;
	double p1, p2, p3, q1, q2, q3;
	double a, b, c, t;

	t_max = -1;
	for (int ii = 0; ii < n_tri; ii++)
	{
		x1 = x(F(ii, 0), 0);
		x2 = x(F(ii, 1), 0);
		x3 = x(F(ii, 2), 0);

		y1 = x(F(ii, 0), 1);
		y2 = x(F(ii, 1), 1);
		y3 = x(F(ii, 2), 1);

		p1 = p(F(ii, 0), 0);
		p2 = p(F(ii, 1), 0);
		p3 = p(F(ii, 2), 0);

		q1 = p(F(ii, 0), 1);
		q2 = p(F(ii, 1), 1);
		q3 = p(F(ii, 2), 1);

		a = p1*q2 - p2*q1 - p1*q3 + p3*q1 + p2*q3 - p3*q2;
		b = p1*y2 - p2*y1 - q1*x2 + q2*x1 - p1*y3 + p3*y1 + q1*x3 - q3*x1 + p2*y3 - p3*y2 - q2*x3 + q3*x2;
		c = x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2;

		t = getSmallestPositiveRealQuadRoot(a, b, c, tol);
		if (t >= 0)
		if ((t_max < 0) | (t_max>t))
			t_max = t;
	}
}


double getSmallestPositiveRealCubicRoot(double a, double b, double c, double d, double tol)
{
	// return negative value if no positive real root is found
	double t = -1;

	if (abs(a) <= tol)
		t = getSmallestPositiveRealQuadRoot(b, c, d, tol);
	else
	{
		complex<double> i(0, 1);
		complex<double> delta0(b*b - 3 * a*c, 0);
		complex<double> delta1(2 * b*b*b - 9 * a*b*c + 27 * a*a*d, 0);
		complex<double> C = pow((delta1 + sqrt(delta1*delta1 - 4.0 * delta0*delta0*delta0)) / 2.0, 1.0 / 3.0);

		complex<double> u2 = (-1.0 + sqrt(3.0)*i) / 2.0;
		complex<double> u3 = (-1.0 - sqrt(3.0)*i) / 2.0;

		complex<double> t1 = (b + C + delta0 / C) / (-3.0*a);
		complex<double> t2 = (b + u2*C + delta0 / (u2*C)) / (-3.0*a);
		complex<double> t3 = (b + u3*C + delta0 / (u3*C)) / (-3.0*a);

		if ((abs(imag(t1))<tol) && (real(t1)>0))
			t = real(t1);
		if ((abs(imag(t2))<tol) && (real(t2)>0) && ((real(t2) < t) || (t < 0)))
			t = real(t2);
		if ((abs(imag(t3))<tol) && (real(t3)>0) && ((real(t3) < t) || (t < 0)))
			t = real(t3);
	}
	return t;
}


void computeInjectiveStepSize_3d(const MatrixXd& F, const MatrixXd& x, const MatrixXd& p, double tol, double& t_max)
{
	int n_tri = F.rows();
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
	double p1, p2, p3, p4, q1, q2, q3, q4, r1, r2, r3, r4;
	double a, b, c, d, t;

	t_max = -1;
	for (int ii = 0; ii < n_tri; ii++)
	{
		x1 = x(F(ii, 0), 0);
		x2 = x(F(ii, 1), 0);
		x3 = x(F(ii, 2), 0);
		x4 = x(F(ii, 3), 0);

		y1 = x(F(ii, 0), 1);
		y2 = x(F(ii, 1), 1);
		y3 = x(F(ii, 2), 1);
		y4 = x(F(ii, 3), 1);

		z1 = x(F(ii, 0), 2);
		z2 = x(F(ii, 1), 2);
		z3 = x(F(ii, 2), 2);
		z4 = x(F(ii, 3), 2);

		p1 = p(F(ii, 0), 0);
		p2 = p(F(ii, 1), 0);
		p3 = p(F(ii, 2), 0);
		p4 = p(F(ii, 3), 0);

		q1 = p(F(ii, 0), 1);
		q2 = p(F(ii, 1), 1);
		q3 = p(F(ii, 2), 1);
		q4 = p(F(ii, 3), 1);

		r1 = p(F(ii, 0), 2);
		r2 = p(F(ii, 1), 2);
		r3 = p(F(ii, 2), 2);
		r4 = p(F(ii, 3), 2);

		a = -p1*q2*r3 + p1*r2*q3 + q1*p2*r3 - q1*r2*p3 - r1*p2*q3 + r1*q2*p3 + p1*q2*r4 - p1*r2*q4 - q1*p2*r4 + q1*r2*p4 + r1*p2*q4 - r1*q2*p4 - p1*q3*r4 + p1*r3*q4 + q1*p3*r4 - q1*r3*p4 - r1*p3*q4 + r1*q3*p4 + p2*q3*r4 - p2*r3*q4 - q2*p3*r4 + q2*r3*p4 + r2*p3*q4 - r2*q3*p4;
		b = -x1*q2*r3 + x1*r2*q3 + y1*p2*r3 - y1*r2*p3 - z1*p2*q3 + z1*q2*p3 + x2*q1*r3 - x2*r1*q3 - y2*p1*r3 + y2*r1*p3 + z2*p1*q3 - z2*q1*p3 - x3*q1*r2 + x3*r1*q2 + y3*p1*r2 - y3*r1*p2 - z3*p1*q2 + z3*q1*p2 + x1*q2*r4 - x1*r2*q4 - y1*p2*r4 + y1*r2*p4 + z1*p2*q4 - z1*q2*p4 - x2*q1*r4 + x2*r1*q4 + y2*p1*r4 - y2*r1*p4 - z2*p1*q4 + z2*q1*p4 + x4*q1*r2 - x4*r1*q2 - y4*p1*r2 + y4*r1*p2 + z4*p1*q2 - z4*q1*p2 - x1*q3*r4 + x1*r3*q4 + y1*p3*r4 - y1*r3*p4 - z1*p3*q4 + z1*q3*p4 + x3*q1*r4 - x3*r1*q4 - y3*p1*r4 + y3*r1*p4 + z3*p1*q4 - z3*q1*p4 - x4*q1*r3 + x4*r1*q3 + y4*p1*r3 - y4*r1*p3 - z4*p1*q3 + z4*q1*p3 + x2*q3*r4 - x2*r3*q4 - y2*p3*r4 + y2*r3*p4 + z2*p3*q4 - z2*q3*p4 - x3*q2*r4 + x3*r2*q4 + y3*p2*r4 - y3*r2*p4 - z3*p2*q4 + z3*q2*p4 + x4*q2*r3 - x4*r2*q3 - y4*p2*r3 + y4*r2*p3 + z4*p2*q3 - z4*q2*p3;
		c = -x1*y2*r3 + x1*z2*q3 + x1*y3*r2 - x1*z3*q2 + y1*x2*r3 - y1*z2*p3 - y1*x3*r2 + y1*z3*p2 - z1*x2*q3 + z1*y2*p3 + z1*x3*q2 - z1*y3*p2 - x2*y3*r1 + x2*z3*q1 + y2*x3*r1 - y2*z3*p1 - z2*x3*q1 + z2*y3*p1 + x1*y2*r4 - x1*z2*q4 - x1*y4*r2 + x1*z4*q2 - y1*x2*r4 + y1*z2*p4 + y1*x4*r2 - y1*z4*p2 + z1*x2*q4 - z1*y2*p4 - z1*x4*q2 + z1*y4*p2 + x2*y4*r1 - x2*z4*q1 - y2*x4*r1 + y2*z4*p1 + z2*x4*q1 - z2*y4*p1 - x1*y3*r4 + x1*z3*q4 + x1*y4*r3 - x1*z4*q3 + y1*x3*r4 - y1*z3*p4 - y1*x4*r3 + y1*z4*p3 - z1*x3*q4 + z1*y3*p4 + z1*x4*q3 - z1*y4*p3 - x3*y4*r1 + x3*z4*q1 + y3*x4*r1 - y3*z4*p1 - z3*x4*q1 + z3*y4*p1 + x2*y3*r4 - x2*z3*q4 - x2*y4*r3 + x2*z4*q3 - y2*x3*r4 + y2*z3*p4 + y2*x4*r3 - y2*z4*p3 + z2*x3*q4 - z2*y3*p4 - z2*x4*q3 + z2*y4*p3 + x3*y4*r2 - x3*z4*q2 - y3*x4*r2 + y3*z4*p2 + z3*x4*q2 - z3*y4*p2;
		d = x1*z2*y3 - x1*y2*z3 + y1*x2*z3 - y1*z2*x3 - z1*x2*y3 + z1*y2*x3 + x1*y2*z4 - x1*z2*y4 - y1*x2*z4 + y1*z2*x4 + z1*x2*y4 - z1*y2*x4 - x1*y3*z4 + x1*z3*y4 + y1*x3*z4 - y1*z3*x4 - z1*x3*y4 + z1*y3*x4 + x2*y3*z4 - x2*z3*y4 - y2*x3*z4 + y2*z3*x4 + z2*x3*y4 - z2*y3*x4;


		t = getSmallestPositiveRealCubicRoot(a, b, c, d, tol);
		if (t >= 0)
		if ((t_max < 0) | (t_max>t))
			t_max = t;
	}
}

void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])
{
	// assign input
	int n_tri = mxGetM(prhs[0]); // # rows of F
	int d_simplex = mxGetN(prhs[0]); // # cols of F
	int dim = d_simplex - 1;	
	int n_vars = mxGetM(prhs[1]); // # rows of x/p
	int n_vert = n_vars / dim;
	const Map<MatrixXd, Aligned> Fmatlab(mxGetPr(prhs[0]), n_tri, d_simplex);
	const Map<MatrixXd, Aligned> x(mxGetPr(prhs[1]), n_vert, dim);
	const Map<MatrixXd, Aligned> p(mxGetPr(prhs[2]), n_vert, dim);
	double *tol;
	tol = mxGetPr(prhs[3]);
	
	// update index numbers to 0-base
	MatrixXd F (Fmatlab);
	F = F.array() - 1;	

	// compute
	double t_max;

	// compute
	if (dim == 2)
		computeInjectiveStepSize_2d(F, x, p, *tol, t_max);
	else if (dim == 3)
		computeInjectiveStepSize_3d(F, x, p, *tol, t_max);
	else
		mexErrMsgIdAndTxt("MATLAB:wrong_dimension", "dim must be either 2 or 3");

	
	// assign outputs
	if (t_max < 0)
		t_max = mxGetInf();
	plhs[0] = mxCreateDoubleScalar(t_max); // functional value
	
}