function AffineT = findAffine (fixPoints, movPoints)
% finding the affine transformation between two sets of points. first we
% normalize the points such that they are centered at origin and second the
% average distance of points to the origin is sqrt(2);
% points are Nx3, where N is the number of points and 3 is the dimension!
N = size(movPoints, 1);

if (N<4)
    error('at least four points are necessary to find the affine transformation');
end



[movPoints, T_mov] = normalizeControlPoints(movPoints);
[fixPoints, T_fix] = normalizeControlPoints(fixPoints);



if(sum(isnan(movPoints(:)) | movPoints(:)==inf ) ~= 0 || sum(isnan(fixPoints(:)) | fixPoints(:)==inf ) ~= 0)
    AffineT = eye(4);
%     warning('Couldn''t find an affine transformation for the given points;');
    return;
end

% fixPoints = [fixPoints ones(size(fixPoints, 1), 1)];
% movPoints = [movPoints ones(size(movPoints, 1), 1)];
% AffineT = [movPoints ones(size(movPoints, 1), 1)]\[fixPoints ones(size(fixPoints, 1), 1)];
A = zeros(3*N, 12);
P = zeros(3*N, 1);

for ii = 1:N
    t = 3*ii;
    A(t-2, :) = [movPoints(ii, 1) movPoints(ii, 2) movPoints(ii, 3) 1    0                0                0                0   0                0                0                0];
    A(t-1, :) = [0               0               0                  0    movPoints(ii, 1) movPoints(ii, 2) movPoints(ii, 3) 1   0                0                0                0];
    A(t, :)   = [0               0               0                  0    0                0                0                0   movPoints(ii, 1) movPoints(ii, 2) movPoints(ii, 3) 1];
    P(t-2)    = fixPoints(ii, 1);
    P(t-1)    = fixPoints(ii, 2);
    P(t)      = fixPoints(ii, 3);
end

[U, D, V] = svd(A);
h         = V*pinv(D)*U'*P;
AffineT   = [h(1) h(2)  h(3)  h(4); ...
             h(5) h(6)  h(7)  h(8); ...
             h(9) h(10) h(11) h(12); ...
             0    0     0     1];
AffineT         = (T_fix'\AffineT*T_mov')';
% AffineT   = H';