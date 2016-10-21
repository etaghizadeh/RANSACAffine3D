function [normPoints, T] = normalizeControlPoints(points)
% points: Nx3, while N: number of points and points are on 3 dimensional
% space.
center = mean(points);
normPoints = bsxfun(@minus, points, center);


dist = sqrt(sum(normPoints.^2,2));
scaleFactor = sqrt(2)/mean(dist);

translation = [1 0 0 -center(1); ...
               0 1 0 -center(2); ...
               0 0 1 -center(3); ...
               0 0 0 1];
scaling     = [scaleFactor 0           0           0; ...
               0           scaleFactor 0           0; ...
               0           0           scaleFactor 0; ...
               0           0           0           1];
           
T = (scaling * translation)';
normPoints = [points ones(size(points,1), 1)] * T;
normPoints = normPoints(:, 1:end-1);