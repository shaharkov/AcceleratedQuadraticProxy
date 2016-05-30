%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Accelerated Quadratic Proxy for Geometric Optimization", SIGGRAPH 2016.
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the author to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%            Meirav Galun (http://www.wisdom.weizmann.ac.il/~/meirav/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = getTensorTranspose(n,m)
% Returns a matrix T such that vec(X')=T*vec(X) for any nxm matrix X.
%
% Example:
% n = 3;
% m = 4;
% X = zeros(n,m);
% X(:)=1:numel(X);
% T = getTensorTranspose(n,m);
% vec(X')
% T*vec(X)

[i,j] = meshgrid(1:n,1:m);
T = sparse(sub2ind([m n],j,i),sub2ind([n m],i,j),1);