%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Accelerated Quadratic Proxy for Geometric Optimization", SIGGRAPH 2016.
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the author to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%            Meirav Galun (http://www.wisdom.weizmann.ac.il/~/meirav/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
cd mex


%% compile all
eigenFolder = 'eigen-eigen-b30b87236a1b';

mex('-largeArrayDims',['-I',eigenFolder],'computeMeshTranformationCoeffsMex.cpp');
mex('-largeArrayDims',['-I',eigenFolder],'computeInjectiveStepSizeMex.cpp');
mex('-largeArrayDims',['-I',eigenFolder],'projectRotationMex.cpp');
mex('-largeArrayDims',['-I',eigenFolder],'computeFunctionalIsoDistMex.cpp');


%% fast routines (require libigl)
libiglFolder = 'libigl-master\include\igl';

mex('-largeArrayDims',['-I',eigenFolder],['-I',libiglFolder],'projectRotationMexFast.cpp');


%%
cd ..