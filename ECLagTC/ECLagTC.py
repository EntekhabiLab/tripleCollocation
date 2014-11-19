'''
Collection of functions for correlated-error lagged triple collocation, a
statistical technique for the error characterization of two co-located 
(in time and space) timeseries representing a noisy version of the same 
variable. For a description of the method, see

Konings, A.G., K.A. McColl, S.H. Alemohammad, Chun-Hsu Su, and 
D. Entekhabi (2014). Error Characterization of Similar Products: 
Correlated-Error Triple Collocation. Submitted. 


Written by Alexandra Konings, konings@alum.mit.edu, 09/2014
'''

from __future__ import division
import numpy as np

def RegressYork(x,y,r,Wx,Wy,maxCount=100):
    '''
    Williamson-York linear regression assuming (possibly correlated) 
    errors on both x- and y-axes. For a good reference on the subject, see
    Cantrell, C. A. (2008), Technical note: Review of methods for linear
    least-squares fitting of data and application to atmospheric chemistry
    problems. Atmospheric Chemistry and Physics, 17: 5477-5487.
    doi://10.5194/acp-8-5477-2008
    
    INPUTS: X and Y are vectors of the independent and dependent variables
    respectively. X and Y are assumed to be the same size
            R is the cross-correlation between errors in X and Y
            Wx and Wy are weights for the x- and y-axes, respectively. If 
            the relative sizes of the errors are known, Wx = 1/var(error_x)
            and similarly for Wy would be a good choice. If they are not 
            known, Wx = 1/var(X) and analogously also works well
    
    OUTPUTS: m is the slope of the regression      
             b is the intercept of the regression        
    '''
    
    U = x - np.nanmean(x)   
    V = y - np.nanmean(y)    
    alpha = np.sqrt(Wx*Wy) 
    m = np.nanmedian(y/x) #initialization
    
    #Do iteration on regression
    done = False
    cnt = 1
    while not done:
        cnt = cnt+1
        mOld = m           
        Wi = Wx*Wy/( m**2*Wy+Wx - 2*m*r*alpha )
    
        beta = Wi*( U*Wy + m*V/Wx - (m*U+V)*r/alpha )
        m = np.sum(Wi*beta*V)/np.sum(Wi*beta*U)
        if np.abs(m-mOld) < 1e-4:
            done = True
        
        if cnt > maxCount:
           raise StandardError('Regression does not converge')
    
    b = np.nanmean(y) - m*np.nanmean(x)
    
    return m, b
    
def ECLagTC(x):
    C = np.cov(x.T)
    
    #Do calculation
    beta1 = 1
    hasData = np.argwhere( ~np.isnan(x[:,0]) & ~np.isnan(x[:,1]) & ~np.isnan(x[:,2]) )
    m31, b31 = RegressYork(x[hasData,0], x[hasData,2],0,1,1)
    TT = C[0,2]/m31
    e11 = C[0,0] - beta1**2*TT
    beta2 = C[1,2]/C[0,2]
    e22 = C[1,1] - beta2**2*TT
    e12 = C[0,1] - beta1*beta2*TT
    TTL = C[0,2]

    #Save for output vectors
    varVec = [e11, e22, e12]
    scale = [TTL, TT, beta2]
    corrTruth1 = 1*np.sqrt(TT)/np.sqrt(C[0,0])
    corrTruth2 = beta2*np.sqrt(TT)/np.sqrt(C[1,1])
    corrTruth = [corrTruth1, corrTruth2]


    return varVec, scale, corrTruth


def LagTC(x):
    '''
    INPUT: x is an N X 3 vector that calculates the errors and scale parameters
    associated with triple co-location of three vectors of timeseries. The
    third is assumed to be a lagged version of the first. 
    OUTPUTS: 
       varVec contains three elements, the covaraince of the errors
       in the first (and third) products, second product, and between the second
       and third product, respectively.
    
       scale returns a 1 x 3 vectors. The first elements is the
       covariance between the truth and lagged truth, respectively. The second 
       element is the variance of the trith, and the third element is the
       magnitude of the scale for the second variable (the first is assumed to
       have beta = 1). 
    
       corr is a 1 x 2 vector containing the correlations between each
       product and the truth. 
    '''

    C = np.cov(x.T)

    #Do calculation
    beta1 = 1
    beta2 = C[1,2]/C[0,2]
    TT = C[0,1]/beta1/beta2
    e11 = C[0,0] - beta1**2*TT
    e22 = C[1,1] - beta2**2*TT
    TTL = C[1,2]/beta1/beta2

    #Save for output vectors
    varVec = [e11, e22, 0]
    scale = [TTL, TT, beta2]
    corrTruth1 = 1*np.sqrt(TT)/np.sqrt(C[0,0])
    corrTruth2 = beta2*np.sqrt(TT)/np.sqrt(C[1,1])
    corrTruth = [corrTruth1, corrTruth2]


    return varVec, scale, corrTruth
    


