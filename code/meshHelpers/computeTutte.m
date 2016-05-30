%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Accelerated Quadratic Proxy for Geometric Optimization", SIGGRAPH 2016.
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the author to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%            Meirav Galun (http://www.wisdom.weizmann.ac.il/~/meirav/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [V0, inds_bF] = computeTutte(F,V,scale)

% compute mesh gradient
[T,areas] = computeMeshTranformationCoeffsMex(F, V);

% constrain boundary to the circle
TR = triangulation(F, V);
temp = freeBoundary(TR);
inds_bF = temp(:,1); %boundary vertices indices
bn_vert = length(inds_bF); %number of boundary verts
tet = 0: (2*pi / bn_vert) : 2*pi;
tet(end)=[];
tet = tet';
[eq_lhs,eq_rhs] = indCoordsToLinearSystem(V(:,1:2), inds_bF, 1*[cos(tet) sin(tet)] );

% compute tutte
wT = spdiag(kron(sqrt(areas), ones(4,1)))*T;
V0 = reshape(solveConstrainedLS(wT,zeros(size(T,1),1),eq_lhs,eq_rhs),[],2);

% global scale (isometric as possible)
if scale
    V0 = ((T*V0(:))\colStack(repmat(eye(2),[1 1 size(F,1)]))) * V0;
end

