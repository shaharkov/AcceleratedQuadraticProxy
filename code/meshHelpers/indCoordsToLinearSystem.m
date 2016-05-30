%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Accelerated Quadratic Proxy for Geometric Optimization", SIGGRAPH 2016.
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the author to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%            Meirav Galun (http://www.wisdom.weizmann.ac.il/~/meirav/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [eq_lhs,eq_rhs] = indCoordsToLinearSystem(X,ind,coords)

eq_lhs_temp = sparse(1:length(ind), ind, 1, length(ind), size(X,1));
eq_lhs = kron(eq_lhs_temp,eye(size(X,2)))*getTensorTranspose(size(X,1),size(X,2));
eq_rhs = colStack(coords');
