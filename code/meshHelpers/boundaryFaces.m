% Code implementing the paper "Injective and Bounded Mappings in 3D".
% Disclaimer: The code is provided as-is and without any guarantees. Please contact the author to report any bugs.
% Written by Noam Aigerman, http://www.wisdom.weizmann.ac.il/~noamaig/

function [F,inds] = boundaryFaces(T,keep)
%compute the boundary faces F of the given tetrahedral mesh 
% keep is an optional boolean vector that forces to keep all faces of 
% any tet that its cooresponding entry in "keep" is non-zero
%returns: 
%F: n x 3 array of all faces from T that are boundary or belong to a tet
% that is marked as 1 in "keep".
%inds(i)==k iff the the boundary face F(i,:) belongs to the tetrahedron 
%T(k,:) 

%alec jacobson's code with slight modification

if isempty(T) || (nargin>1 && isempty(keep));
    F=[];
    inds=[];
    return;
end
if nargin==1
    keep=zeros(size(T,1),1);
end
  % BOUNDARY_FACES
  % F = boundary_faces(T)
  % Determine boundary faces of tetrahedra stored in T
  %
  % Input:
  %  T  tetrahedron index list, m by 4, where m is the number of tetrahedra
  %
  % Output:
  %  F  list of boundary faces, n by 3, where n is the number of boundary faces
  %

  % get all faces
  allF = [ ...
    T(:,1) T(:,2) T(:,3); ...
    T(:,1) T(:,3) T(:,4); ...
    T(:,1) T(:,4) T(:,2); ...
    T(:,2) T(:,4) T(:,3)];
    keep=[keep;keep;keep;keep]~=0;
  % sort rows so that faces are reorder in ascending order of indices
  sortedF = sort(allF,2);
  % determine uniqueness of faces
  [u,~,n] = unique(sortedF,'rows');
  % determine counts for each unique face
  counts = accumarray(n(:), 1);
  % extract faces that only occurred once
  sorted_exteriorF = u(counts == 1,:);
  % find in original faces so that ordering of indices is correct
  use=ismember(sortedF,sorted_exteriorF,'rows')|keep;
  F = allF(use,:);
  
if nargout>1
    inds=1:size(T,1);
  inds=[inds inds inds inds]';

  inds=inds(use);
end
end
