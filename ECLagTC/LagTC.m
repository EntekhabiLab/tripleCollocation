function [varVecRec scaleRec corrRec ratCheck] = LagTC(x)
% INPUT: X is an N X 3 vector that calculates the errors and scale parameters
% associated with triple co-location of three vectors of timeseries. The
% third is assumed to be a lagged version of the first. 
% OUTPUTS: 
%   varVec contains three elements, the covariance of the errors
%   in the first (and third) products, second product, and between the second
%   and third product, respectively.
%
%   scale returns a 1 x 3 vectors. The first elements is the
%   covariance between the truth and lagged truth, respectively. The second 
%   element is the variance of the trith, and the third element is the
%   magnitude of the scale for the second variable (the first is assumed to
%   have beta = 1). 
%
%   corr is a 1 x 2 vector containing the correlations between each
%   product and the truth. 
%
%Written by Alexandra Konings, konings@alum.mit.edu, 09/2014

C = nancov(x);

% Main calculation
beta1 = 1;
beta2 = C(2,3)/C(1,3);
TT = C(1,2)/beta1/beta2;
e11 = C(1,1) - beta1^2*TT;
e22 = C(2,2) - beta2^2*TT;
TTL = C(2,3)/beta1/beta2;

% Save to output vectors
varVecRec = [e11 e22 0];
scaleRec = [TTL TT beta2 ];
corrTruth1 = 1*sqrt(TT)/sqrt(C(1,1));
corrTruth2 = beta2*sqrt(TT)/sqrt(C(2,2));
corrRec = [corrTruth1 corrTruth2];
ratCheck = (C(3,3)-beta1.^2*TT)/e11;

