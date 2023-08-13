function [x1,y1,x2,y2,eval,evec] = SegmentData(x,y,sigma_d,DistanceCutoff)
% This function applys 2D spectral data segmentation to the input data (x,
% y) using a spatial weight sigma_d.

% Tolerance for removing very small values in expotential fucntion
tol = sqrt(eps);




DistanceMatrix = pdist([x y]);
A = (squareform(DistanceMatrix));
DistanceMatrixMask = A < DistanceCutoff;
DistanceMatrix = A.*(DistanceMatrixMask);
DistanceMatrix(DistanceMatrix == 0) = 1e10;
WeightingMatrix = exp((-DistanceMatrix.^2) / (2 * sigma_d^2));

D = sparse(diag(sum(WeightingMatrix,2)));
% Form the Laplacian
L = sparse(D - WeightingMatrix);
% Find the eigenvectors and eigenvalues

% [evec,eval] = eig(L, 'vector');

% Ensure the eigenvalues and eigenvectors are sorted in ascending order
idx = find(abs(L) < tol);
L(idx) = 0;
[evec,eval] = eigs(L,2,'sm');


ClusterVar = evec(:,2) <= 0;
x1 = x(ClusterVar);
y1 = y(ClusterVar);
x2 = x(~ClusterVar);
y2 = y(~ClusterVar);