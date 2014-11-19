tripleCollocation
=================

A collection of implementations functions implementing variants of the (extended) triple collocation method (TC), a statistical tool to determine the RMSE and correlation coefficients associated with three coincident noisy estimates of the same true series. Functions are written in Matlab and 
(in some cases) also available with Python translation. 

These functions are freely available, with attributed required for redistribution. 

The variants currently provided are 

### 1) Error-correlated lagged triple co-location (ECLagTC)
citation:
>Konings, A.G., K.A. McColl, S.H. Alemohammad, D. Entekhabi, C-H. Su (2014). Error Characterization of Similar Data Sets: Triple Collocation with Correlated Errors, Submitted to Geophys. Res. Lett.

A variant of triple-collocation that uses a lagged version of one of the three data sets as a third dataset, so that only two independent datasets are needed. Unlike other versions of TC, the errors of different data sets are not assumed to be independent, and their cross-correlations are explicitly accounted for. 
The ECLagTC method uses an errors-in-variables regression due to York (1968). Code for performing this regression is also included

Code for an alternative version of this approach without cross-correlated errors (LagTC) is also provided . This variant of TC is based in the recognition that TC can be generalized to the statistical framework of instrumental variables. 
citation: 
>Su, C.-H., D. Ryu, W. T. Crow, and A. W. Western (2014), Beyond triple collocation: Applications to soil moisture monitoring, J. Geophys. Res. Atmos., 119, 6419–6439, doi:10.1002/2013JD021043.




