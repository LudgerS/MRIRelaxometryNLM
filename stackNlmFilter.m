function filtered =  stackNlmFilter(stack, sigma, params)
% following Manjon 2010 (or buades 2005) regarding variable names
% missing preselection and automatic sigma estimation
%
% Assumes that the third dimension of stack is the relaxometry dimension
%
% Written by Ludger Starke; Max Delbr�ck Center for Molecular Medicine in
% the Helmholtz Association, Berlin; 21-01-25
%
% License: GNU GPLv3 

if numel(size(stack)) ~= 3
    error('wrong dimension of input data')
end

nTEs = size(stack, 3);

%% parameters

n = params.centerDistance;      % block center distance
m = params.blockRadius;         % (2m + 1)^2 is the block area, needs to be m > n/2
s = params.searchRadius;        % (2s + 1)^2 is the number of searched blocks

beta = params.beta;

hSq = 2*beta*(2*m + 1)^2*nTEs*sigma^2;  % filtering parameter controlling the filter strength

%% pad image

originalDim = [size(stack, 1), size(stack, 2)];

toSymmetricGrid = n - 1 - mod(originalDim, n);
toFitBlocks = max(m - (n - 1), 0);

paddedData = padarray(stack, [toFitBlocks*[1, 1], 0], 'replicate', 'pre');
paddedData = padarray(paddedData, [toFitBlocks + toSymmetricGrid, 0], 'replicate', 'post');

indexOrigX = (1 + toFitBlocks):(toFitBlocks + originalDim(1));
indexOrigY = (1 + toFitBlocks):(toFitBlocks + originalDim(2));


%% compute block centers

dim = size(paddedData);

blockCentersX = max(m + 1, n):n:(dim(1) - m);
blockCentersY = max(m + 1, n):n:(dim(2) - m);


%% apply filter

overlapCounter = zeros(size(paddedData));
u = zeros(size(paddedData));

counter = 1;

% loop through to be filtered blocks along x
for ii = 1:numel(blockCentersX)
    
    if (ii/numel(blockCentersX) > counter/10)
        fprintf('%d%% ', counter*10)
        counter = counter + 1;
    end
    
    % loop through to be filtered blocks along y
    for jj = 1:numel(blockCentersY)
        
        blockCenter = [blockCentersX(ii), blockCentersY(jj)];
        
%         fprintf('X: %d - Y: %d\n', blockCenter(1), blockCenter(2))
        
        blockMaskB = getBlockMask(blockCenter, dim, m);
        B = paddedData(blockMaskB);
        
        temp = zeros(size(B));
        sumWeights = 0;

%         weightMatrix = zeros(dim);
%         blockCenter = [blockCentersX(10), blockCentersY(10)];

        % compare to be filtered block (B) to neighbors and compute weights
        for kk = max(1, ii - s):min(numel(blockCentersX), ii + s)
            for ll = max(1, jj - s):min(numel(blockCentersY), jj + s)
                
                blockCenter = [blockCentersX(kk), blockCentersY(ll)];
                [rangeX, rangeY] = getBlockMaskRange(blockCenter, dim, m);
                compB = paddedData(rangeX, rangeY, :);
                weight = computeWeight(B, compB(:), hSq);
                
%                 fprintf('inner Loop X: %d - Y: %d - Weight: %f\n', blockCenter(1), blockCenter(2), weight)
%                 weightMatrix(blockCenter(1), blockCenter(2)) = weight;
                
                temp = temp + weight*compB(:);
                sumWeights = sumWeights + weight;
            
            end
        end
        
%         imshow(imresize(weightMatrix, 4, 'method', 'nearest'), [])
%         drawnow
        
        temp = temp/sumWeights;
        
        u(blockMaskB) = u(blockMaskB) + temp;
        overlapCounter = overlapCounter + blockMaskB;
        
    end
end

u = u./overlapCounter;

filtered = u(indexOrigX, indexOrigY, :);

fprintf('100%%\n')





%% -- subfunctions -------------------------------------------------------

% compute block masks
function blockMask = getBlockMask(blockCenter, dim, m)

blockMask = false(dim);

rangeX = max(1, blockCenter(1) - m):min(dim(1), blockCenter(1) + m);
rangeY = max(1, blockCenter(2) - m):min(dim(2), blockCenter(2) + m);

blockMask(rangeX, rangeY, :) = true;


% compute block masks as range (faster than logical indexing)
function [rangeX, rangeY] = getBlockMaskRange(blockCenter, dim, m)

rangeX = max(1, blockCenter(1) - m):min(dim(1), blockCenter(1) + m);
rangeY = max(1, blockCenter(2) - m):min(dim(2), blockCenter(2) + m);


% compute weights
function weight = computeWeight(B, compB, hSq)
% not normalized

weight = exp(-(B - compB)'*(B - compB)/hSq);






