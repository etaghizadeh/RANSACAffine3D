function AffineT = findAffine (fixPoints, movPoints)
% finding the affine transformation between two sets of points. first we
% normalize the points such that they are centered at origin and second the
% average distance of points to the origin is sqrt(2);
% points are Nx3, where N is the number of points and 3 is the dimension!
N = size(movPoints, 1);

if (N<4)
    error('at least three points are necessary to find the affine transformation');
end

[movPoints, T_mov] = normalizeControlPoints(movPoints);
[fixPoints, T_fix] = normalizeControlPoints(fixPoints);

% check for devision to zero!
if(sum(isnan(movPoints(:)) | movPoints(:)==inf ) ~= 0 || sum(isnan(fixPoints(:)) | fixPoints(:)==inf ) ~= 0)
    AffineT = eye(4);
    warning('Couldn''t find an affine transformation for the given points;');
    return;
end

A = zeros(3*N, 12);
P = zeros(3*N, 1);

for ii = 1:N
    A(3*ii-2, :) = [movPoints(ii, 1) movPoints(ii, 2) movPoints(ii, 3) 1    0                0                0                0   0                0                0                0];
    A(3*ii-1, :) = [0               0               0                  0    movPoints(ii, 1) movPoints(ii, 2) movPoints(ii, 3) 1   0                0                0                0];
    A(3*ii, :)   = [0               0               0                  0    0                0                0                0   movPoints(ii, 1) movPoints(ii, 2) movPoints(ii, 3) 1];
    P(3*ii-2)    = fixPoints(ii, 1);
    P(3*ii-1)    = fixPoints(ii, 2);
    P(3*ii)      = fixPoints(ii, 3);
end

[U, D, V] = svd(A);
h         = V*pinv(D)*U'*P;
H         = [h(1) h(2)  h(3)  h(4); ...
             h(5) h(6)  h(7)  h(8); ...
             h(9) h(10) h(11) h(12); ...
             0    0     0     1];
% if the points are collinear!
if(rank(H)<4)
    AffineT = eye(4);
    warning('Couldn''t find an affine transformation for the given points;');
    return;
end
AffineT         = (T_fix'\H*T_move')';
