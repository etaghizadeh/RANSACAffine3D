function [T, varargout] = estimateGeometricTransform3D(fix, mov, varargin)
% [T, inliers] = estimateGeometricTransform3D(fix, mov, options)
% fix: fix landmarks  Mx3
% mov: moving landmarks  Mx3 (in correspondence with fix)
% options:
%         .numIterations:  number of iterations for Ransac (default 100)
%         .thr          :  threshold of fitting (default 5)
%         .d            :  minimum number of points for affin transformation(default 3)
%         .reflection   : 'best', 0, 1 => find the optimum based on 

% [fix, T_fix] = normalizeControlPoints(fix);
% [mov, T_mov] = normalizeControlPoints(mov);

if (nargin == 2)
    options.numIterations = 1000;
	options.thr           = 5;
	options.d             = 4;
    options.reflection    = 'best';
else
    if (nargin == 3)
        options = varargin{1};
    end
end
if (~isfield(options, 'numIterations'))
    options.numIterations = 1000;
end
if (~isfield(options, 'thr'))
    options.thr = 5;
end
if (~isfield(options, 'd'))
    options.d = 4;
end
if (~isfield(options, 'reflection'))
    options.reflection    = 'best';
end

if size(fix, 1)<options.d
    error('not enought number of points are included ...');
end

if size(fix, 1) ~= size(mov, 1)
    error('number of points should be equal for the fixed and moving points.');
end

numPoints = size(fix, 1);
numInliers = 0;
Inliers = [];
distMin = inf;
indices = 1:numPoints;
for ii = 1:options.numIterations
    idx = randperm(numPoints, options.d);
    samplesFix = fix(idx, :);
    samplesMov = mov(idx, :);
    
    AffineT    = findAffine(samplesFix, samplesMov);
    movRigid   = [mov ones(size(mov, 1), 1)] * AffineT;
    movRigid   = movRigid(:, 1:3);
%     [~, ~, Transformation] = procrustes(samplesFix, samplesMov, 'reflection', options.reflection);
%     movRigid = bsxfun(@plus, Transformation.b *  mov * Transformation.T, Transformation.c(1, :));
    X = fix - movRigid;
    dist = sqrt(sum(X.^2, 2));
    ind = indices(dist<options.thr);
    if (numel(ind)>numInliers)
        numInliers = numel(ind);
        Inliers    = ind;
        distMin = mean(dist);
        if numInliers == size(fix, 1)
            break;
        end
    else if(~isempty(ind) && numel(ind) == numInliers && mean(dist) < distMin )
            distMin = mean(dist);
            Inliers    = ind;
        end
    end
end

% fix = [fix ones(size(fix, 1), 1)]*T_fix;
% fix = fix(:, 1:3);
% mov = [mov ones(size(mov, 1), 1)]*T_mov;
% mov = mov(:, 1:3);

samplesFix = fix(Inliers, :);
samplesMov = mov(Inliers, :);
[~, ~, Transformation] = procrustes(samplesFix, samplesMov, 'reflection', options.reflection);
movRigid = bsxfun(@plus, Transformation.b *  mov * Transformation.T, Transformation.c(1, :));
X = fix - movRigid;
dist = sqrt(sum(X.^2, 2));
ind = indices(dist<options.thr);
if (numel(ind)>numInliers || (numel(ind) == numInliers && mean(dist) < distMin))
    Inliers = ind;
end


samplesFix = fix(Inliers, :);
samplesMov = mov(Inliers, :);
[~, ~, T] = procrustes(samplesFix, samplesMov, 'scaling', false, 'reflection', options.reflection);
if nargout==2
    varargout{1} = Inliers;
end




% figure
% for ii = 1:size(imageEdge, 3)
%     title(num2str(ii))
%     imshow(imageEdge(:, :, ii));
%     pause(0.1)
% end