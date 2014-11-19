function [m b] = regYork(x,y,r,Wx,Wy,maxCount)
%Williamson-York linear regression assuming (possibly correlated) errors on
%both x- and y-axes. For a good reference on the subject, see
%Cantrell, C. A. (2008), Technical note: Review of methods for lin-
%ear least-squares fitting of data and application to atmospheric chemistry
%problems. Atmospheric Chemistry and Physics, 17: 5477-5487.
%doi://10.5194/acp-8-5477-2008
%
%INPUTS: X and Y are vectors of the independent and dependent variables
%respectively. X and Y are assumed to be the same size
%        R is the cross-correlation between errors in X and Y
%        Wx and Wy are weights for the x- and y-axes, respectively. If the
%        relative sizes of the errors are known, Wx = 1/var(error_x) and
%        similarly for Wy would be a good choice. If they are not known, 
%        Wx = 1/var(X) and analogously also works well
%
%OUTPUTS: m is the slope of the regression      
%         b is the intercept of the regression        
%
%
%Written by Alexandra Konings, konings@alum.mit.edu, 09/2014

if nargin == 5
    maxCount = 100;
elseif nargin < 5
    display('Please use 5 input arguments. See help RegressYork')
end

U = x - nanmean(x);     
V = y - nanmean(y);    
alpha = sqrt(Wx.*Wy); 
m = nanmedian(y./x); %initialization

%Do iteration on regression
done = 0;
cnt = 1;
while ~done
    cnt = cnt+1;
    mOld = m;            
    Wi = Wx.*Wy./(m^2*Wy+Wx-2*m*r*alpha);   

    beta = Wi.*(U./Wy + m.*V./Wx - (m*U+V)*r./alpha);
    m = sum(Wi.*beta.*V)./sum(Wi.*beta.*U);
    if abs(m-mOld) < 1e-4
        done = 1;
    end
    
      if cnt > maxCount
          error('Regression does not converge')
          m = NaN;
          b = NaN;
      end
end  

b = nanmean(y) - m*nanmean(x);

