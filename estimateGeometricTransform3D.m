function [T, varargout] = estimateGeometricTransform3D(fix, mov, varargin)
% [T, inliers] = estimateGeometricTransform3D(fix, mov, options)
% fix: fix landmarks  Mx3
% mov: moving landmarks  Mx3 (in correspondence with fix)
% options:
%         .numIterations:  number of iterations for Ransac (default 100)
%         .thr          :  threshold of fitting (default 5)
%         .d            :  minimum number of points for affin transformation(default 3)
%         .reflection   : 'best', 0, 1 => find the optimum based on
%         .type        " 'random' or an 'exhaustive' search (in case of
%                        exhaustive, the .numIterations would be set to
%                        possible number of permutations. Avoid exhaustive
%                        search in case of many number of points!
%         .parallel     : true/false (default true) run it in parallel

if (nargin == 2)
    options.numIterations = 1000;
    options.thr           = 3;
    options.d             = 4;
    options.reflection    = 'best';
    options.type          = 'random';
    options.parallel      = true;
else
    if (nargin == 3)
        options = varargin{1};
    end
end
if (~isfield(options, 'numIterations'))
    options.numIterations = 1000;
end
if (~isfield(options, 'thr'))
    options.thr = 3;
end
% thr = options.thr.*sqrt(VAR)
% options.thr = sqrt(sum((options.thr./STD).^2));
if (~isfield(options, 'd'))
    options.d = 4;
end
if (~isfield(options, 'reflection'))
    options.reflection    = 'best';
end
if (~isfield(options, 'type'))
    options.type    = 'random';
end
if (~isfield(options, 'parallel'))
    options.parallel    = true;
end
if size(fix, 1)<options.d
    error('not enought number of points are included ...');
end

if size(fix, 1) ~= size(mov, 1)
    error('number of points should be equal for the fixed and moving points.');
end

numPoints = size(fix, 1);



if strcmpi(options.type, 'random')
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
    
    clusters = kmeans(fix, 4);
    fix1 = fix(clusters==1, :);
    fix2 = fix(clusters==2, :);
    fix3 = fix(clusters==3, :);
    fix4 = fix(clusters==4, :);
    mov1 = mov(clusters==1, :);
    mov2 = mov(clusters==2, :);
    mov3 = mov(clusters==3, :);
    mov4 = mov(clusters==4, :);
    
    % give a chance for points a little bit furthur. It might be helpful to
    % have points that are far from the other points in the landmark sets,
    % while they usually have lower chance in the random selection process.
    for ii = 1:100
        idx1 = randi(size(fix1, 1));
        idx2 = randi(size(fix2, 1));
        idx3 = randi(size(fix3, 1));
        idx4 = randi(size(fix4, 1));
        
        samplesFix = [fix1(idx1, :); fix2(idx2, :); fix3(idx3, :); fix4(idx4, :)];
        samplesMov = [mov1(idx1, :); mov2(idx2, :); mov3(idx3, :); mov4(idx4, :)];

        AffineT    = findAffine(samplesFix, samplesMov);
        movRigid   = [mov ones(size(mov, 1), 1)] * AffineT;
        movRigid   = movRigid(:, 1:3);
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
else if strcmpi(options.type, 'exhaustive')
        P = nchoosek(1:numPoints, options.d);
        Parallel = options.parallel;
        if Parallel
            try
                p = gcp('nocreate'); % If no pool, do not create new one.
                if isempty(p)
                    parpool();
                    p = gcp('nocreate');
                end
            catch
                Parallel = false;
            end
        end
        if Parallel
            Thr = options.thr;
            numWorkers = min(25*p.NumWorkers, size(P, 1));
            PTmp{numWorkers} = 0;
            Fix{numWorkers} = fix;
            Mov{numWorkers} = mov;
            Sz = ceil(size(P, 1)/numWorkers);
            for ii = 1:numWorkers
                PTmp{ii} = P((ii-1)*Sz+1 : min(ii*Sz, size(P, 1)), :);
                Fix{ii}  = fix;
                Mov{ii}  = mov;
            end
            IN{numWorkers} = 0;
            distMin{numWorkers} = inf;
            numInliers{numWorkers} = 0;
            parfor jj = 1:numWorkers
                indices = 1:numPoints;
                numInliers{jj} = 0;
                fix = Fix{jj};
                mov = Mov{jj};
                P = PTmp{jj};
                distMin{jj} = inf; 
                for ii = 1:size(PTmp{jj}, 1);
                    idx = P(ii, :);
                    samplesFix = fix(idx, :);
                    samplesMov = mov(idx, :);

                    AffineT    = findAffine(samplesFix, samplesMov);
                    movRigid   = [mov ones(size(mov, 1), 1)] * AffineT;
                    movRigid   = movRigid(:, 1:3);
                    X = fix - movRigid;
                    dist = sqrt(sum(X.^2, 2));
                    ind = indices(dist<Thr);
                    if (numel(ind)>numInliers{jj})
                        numInliers{jj} = numel(ind);
                        IN{jj}    = ind;
                        distMin{jj} = mean(dist);
                        if numInliers{jj} == size(fix, 1)
                            break;
                        end
                    else if(~isempty(ind) && numel(ind) == numInliers{jj} && mean(dist) < distMin{jj} )
                            distMin{jj} = mean(dist);
                            IN{jj}    = ind;
                        end
                    end
                end
            end
            clear P PTmp;
            numTmp = max(cell2mat(numInliers));
            idx = find(cell2mat(numInliers) == numTmp);
            numInliers = numTmp;
            if numel(idx)>1
                distMin = cell2mat(distMin);
                [distMin, ind] = min(distMin(idx));
                idx = idx(ind);
            else
                distMin = distMin{idx};
            end
            Inliers = IN{idx};
            fix = Fix{1};
            mov = Mov{1};
        else
            numInliers = 0;
            Inliers = [];
            distMin = inf;
            indices = 1:numPoints;
            a = axes;
            for ii = 1:size(P, 1);
                idx = P(ii, :);
                AffineT    = findAffine(fix(idx, :), mov(idx, :));
                movRigid   = [mov ones(size(mov, 1), 1)] * AffineT;
                movRigid   = movRigid(:, 1:3);
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
                else
                    if(~isempty(ind) && numel(ind) == numInliers && mean(dist) < distMin )
                        distMin = mean(dist);
                        Inliers    = ind;
                    end
                end
            end
        end
    end
end

indices = 1:numPoints;
samplesFix = fix(Inliers, :);
samplesMov = mov(Inliers, :);
AffineT    = findAffine(samplesFix, samplesMov);
movRigid   = [mov ones(size(mov, 1), 1)] * AffineT;
movRigid   = movRigid(:, 1:3);
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

