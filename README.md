# quantileFit

## Quantile regression
This function allows for quantile regression analysis on linear and non linear models

```MATLAB
% INPUTS:
%    - xdata: vector of observations x coordinates
%    - ydata: vector of observations y coordinates
%    - funMod: model to fit expressed as an anonymous function
%    - tau: percentile for the quantile regression (between 0 and 1)
%    - lb and ub: lower and upper bounds for the fitting
[fitObj,gof,o] = quantileFit(xdata,ydata,funMod,tau,lb,ub);
```

Inline-style: 
![alt text](https://github.com/PabRD/quantileFit/tree/main/quantileFit_Example.svg)
