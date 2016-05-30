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
num_iter = 500;
visualizeIterations = false;
TolX = 1e-10;
TolFun = 1e-6;


%% generate problem
n=500;
theta = 170;

n_x = n;
n_y = ceil(n/20.0);

[x, y] = meshgrid(1:n_x, 1:n_y);
x = x(:); y = y(:);
V = [x, y];
F = delaunay(V);

anchors_fix = find(x<2);
anchor_fix_coords = V(anchors_fix,:);
anchors_move = find(x==n_x);

ctr = mean(V);
theta = theta *pi/180;
rot = [cos(theta) sin(theta); -sin(theta) cos(theta)];
anchor_move_coords = bsxfun(@plus,bsxfun(@minus,V(anchors_move,:),ctr)*rot,ctr);

% compute corresponding linear system
[eq_lhs_fix,eq_rhs_fix] = indCoordsToLinearSystem(V,anchors_fix,anchor_fix_coords);
[eq_lhs_move,eq_rhs_move] = indCoordsToLinearSystem(V,anchors_move,anchor_move_coords);


eq_lhs = [eq_lhs_fix; eq_lhs_move];
eq_rhs = [eq_rhs_fix; eq_rhs_move];


%% init - poisson
T = computeMeshTranformationCoeffsMex(F, V);
V0 = reshape(solveConstrainedLS(T, zeros(size(T,1),1), eq_lhs, eq_rhs),[],2);


%% setup optimization problem
optimProblem = OptimProblemArap(V, F, eq_lhs, eq_rhs, V0); % actually start with V0


%% setup solver
solver{1} = OptimSolverAcclQuadProx('AQP', optimProblem, true, true, false); solver{1}.setKappa(100);
solver{2} = OptimSolverAcclQuadProx('Global-Local', optimProblem, false, true, false); % AQP without acceleration reduces to global-local (see Appendix B)
n_solvers = length(solver);


%% solve
for ii = 1:n_solvers
    logs{ii} = solver{ii}.solveTol(TolX, TolFun ,num_iter);
end


%% visualize
visualizeSolversForExamples