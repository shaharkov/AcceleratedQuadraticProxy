%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Accelerated Quadratic Proxy for Geometric Optimization", SIGGRAPH 2016.
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the author to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%            Meirav Galun (http://www.wisdom.weizmann.ac.il/~/meirav/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% init
rng(1)
clear
initialize
close all


%% parameters
num_iter = 3000;
visualizeIterations = true;
TolX = 1e-10;
TolFun = 1e-6;


%% generate problem
% load mesh
[V, F] = read_off('data\schaeffer_cow.off'); 
V=V'; F=F'; F = F(:,[1 3 2]);
% initialize with tutte 
[V0, inds_bF] = computeTutte(F,V,true); % global scale to be as ismetric as possible
% constrain the centroid
eq_lhs_centroid = kron(eye(2),sparse(1,inds_bF,1,1,size(V,1)));
eq_rhs_centroid = [0;0];


%% setup optimization problem
optimProblem = OptimProblemIsoDist(V, F, eq_lhs_centroid, eq_rhs_centroid, V0);


%% setup solver
solver{1} = OptimSolverAcclQuadProx('AQP', optimProblem, true, true, true); solver{1}.setKappa(1000);
n_solvers = length(solver);


%% solve
for ii = 1:n_solvers
    logs{ii} = solver{ii}.solveTol(TolX, TolFun ,num_iter);
end


%% visualize
visualizeSolversForExamples;
