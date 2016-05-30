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
num_iter = 100;
visualizeIterations = true;
TolX = 1e-10;
TolFun = 1e-6;


%% load problem
load('data\jacobson_boy.mat');


%% setup optimization problem
optimProblem = OptimProblemArap(V, F, eq_lhs, eq_rhs, []);


%% setup solver
solver{1} = OptimSolverAcclQuadProx('AQP', optimProblem, true, true, false); solver{1}.setKappa(100); % no need for line search in ARAP
n_solvers = length(solver);


%% solve
for ii = 1:n_solvers
    logs{ii} = solver{ii}.solveTol(TolX, TolFun ,num_iter);
end


%% visualize
visualizeSolversForExamples;
