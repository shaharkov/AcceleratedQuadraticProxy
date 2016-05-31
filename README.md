Accelerated Quadratic Proxy for Geometric Optimization
====

A Matlab implementation of the paper ["Accelerated Quadratic Proxy for Geometric Optimization"](http://www.wisdom.weizmann.ac.il/~shaharko/AcceleratedQuadraticProxy.html).

----
The class `OptimSolverAcclQuadProx.m` implements the accelerated quadratic proxy optimization algorithm described in the paper for the minimization of geometric energies.

This package includes a few examples:
- `example_deformation_2d_gecko_IsoDist.m` demonstrates the minimization of the Isometric Distortion energy for shape deformation (Figure 1 in the paper).
- `example_parameterization_cow_IsoDist.m` demonstrates the minimization of the Isometric Distortion energy for parameterization (the cow of Figure 12 in the paper).
- `example_deformation_2d_bar_ARAP.m` and `example_deformation_3d_boy_ARAP.m` demonstrate the minimization of the As-Rigid-As-Possible energy for shape deformation (see https://goo.gl/skatVH and https://goo.gl/iYXJaP for video clips).

**Compatibility and dependencies:** The code was tested with Matlab (2015a). The code depends on a few MEX functions. Windows (x64) binaries are included; they are compiled with Intel C++ Composer XE 2016 with Microsoft Visual Studio 2013; compilation requires [Eigen](http://eigen.tuxfamily.org/); fast implementation of As-Rigid-As-Possible also relies on [libigl](https://github.com/libigl/libigl). The source code is provided under the `mex/` folder; run `compileAllMex.m` to compile all mex files (only tested under windows).

**Disclaimer:**
The code is provided as-is for academic use only and without any guarantees. Please contact the authors to report any bugs. 
Written by [Shahar Kovalsky](http://www.wisdom.weizmann.ac.il/~shaharko/) and [Meirav Galun](http://www.wisdom.weizmann.ac.il/~/meirav/).