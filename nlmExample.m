% nlmExample.m
%
% Written by Ludger Starke; Max Delbrück Center for Molecular Medicine in
% the Helmholtz Association, Berlin; 21-01-25
%
% License: GNU GPLv3 

clear, close all

%% load build-in brain data set

load mri
reference = double(D(:,:,1,12));

%% simulate noisy data
peakSnr = 30;
sigma = max(reference(:))/peakSnr;
noisyImage = abs(reference + sigma*(randn(size(reference)) + 1i*randn(size(reference))));

% generate stack of images as artifical relaxometry data

noisyImage2 = abs(0.7*reference + sigma*(randn(size(reference)) + 1i*randn(size(reference))));
imageStack = cat(3, noisyImage, noisyImage2);


%% apply 2D filter

% set nlm parameters
params.centerDistance = 1;      % n - block center distance
params.blockRadius = 1;         % m - (2m + 1)^2 is the block area, needs to be m > n/2
params.searchRadius = 10;       % s - (2s + 1)^2 is the number of searched blocks
params.beta = 0.5;


filtered =  nlmFilter2D(noisyImage, sigma, params);
stackFiltered =  stackNlmFilter(imageStack, sigma, params);



%% show results

imshow([reference, noisyImage, filtered, stackFiltered(:,:,1)], [], 'border', 'tight')





