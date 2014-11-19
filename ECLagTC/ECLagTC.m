function [varVec scale corr] = ECLagTC(x)
% INPUT: X is an N X 3 vector that calculates the errors and scale parameters
% associated with triple co-location of three vectors of timeseries. The
% third is assumed to be a lagged version of the first. 
% OUTPUTS: 
%   varVecRec contains three elements, the covaraince of the errors
%   in the first (and third) products, second product, and between the second
%   and third product, respectively.
%
%   scaleRec returns a 1 x 3 vectors. The first elements is the
%   covariance between the truth and lagged truth, respectively. The second 
%   element is the variance of the trith, and the third element is the
%   magnitude of the scale for the second variable (the first is assumed to
%   have beta = 1). 
%
%   corrRec is a 1 x 2 vector containing the correlations between each
%   product and the truth. 
%
%Written by Alexandra Konings, konings@alum.mit.edu, 09/2014

C = nancov(x);

beta1 = 1;
hasData = find(~isnan(x(:,1)) & ~isnan(x(:,3)));
m31 = regYork(x(hasData,1),x(hasData,3),0,1,1);
TT = C(1,3)/m31;
e11 = C(1,1) - TT;
beta2 = C(2,3)/C(1,3);
e22 = C(2,2) - beta2^2*TT;
e12 = C(1,2) - beta2*TT;
TTL = C(1,3);

%Save for output vectors
varVec = [e11 e22 e12];
scale = [TTL TT beta2 ];
corrTruth1 = 1*sqrt(TT)/sqrt(C(1,1));
corrTruth2 = beta2*sqrt(TT)/sqrt(C(2,2));
corr = [corrTruth1 corrTruth2];

