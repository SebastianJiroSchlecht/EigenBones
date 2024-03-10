% SimpleExample.m

clear all; close all;
randn('seed',0);

%------------------------------------------------------------------------
%  ground truth parameters
%------------------------------------------------------------------------
j = sqrt(-1);
f1 = [1.6180    0.61801];  f2 = [1.6180    0.61801*j];
lambda1 = conv(f1,conj(fliplr(f1))); lambda2 = conv(f2,conj(fliplr(f2)));

%------------------------------------------------------------------------
%  estimation from data
%------------------------------------------------------------------------
Lu=1e6;
u = (randn(2,Lu)+j*randn(2,Lu))/sqrt(2);
x = zeros(2,Lu);
X(1,:) = filter(f1,1,u(1,:)); X(2,:) = filter(f2,1,u(2,:));
Rhat1 = SpaceTimeCovMatEst(X(:,1:Lu/100),1);        % 100 sample estimate
Rhat2 = SpaceTimeCovMatEst(X(:,1:Lu),1);        % 10000 sample estimate

% save PEVDToyProblem3Estimates.mat Rhat1 Rhat2

