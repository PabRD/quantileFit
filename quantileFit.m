function [fitObj,gof,o] = quantileFit(xdata,ydata,funMod,tau,lb,ub)
%QUANTILEFIT   Fit a curve with quantile regression.
%
%   FO = QUANTILEFIT(X, Y, FT, TAU) creates a fit object, FO, that encapsulates the
%   result of fitting the model specified by the fittype FT to the data X,
%   Y.
%
%   -- X must be a matrix with one (curve fitting) column. 
%
%   -- Y must be a column vector with the same number of rows as X.
%
%   -- FT is a string or a anonymous function specifying the model to fit.
%
%   If FT is a string, then it may be:
%
%       FITTYPE           DESCRIPTION
%       'linear'          Linear model
%       '2-param'         Two parameter critical power model 
%       '3-param'         Three parameter critical power model
%       'power'           Single term exponential model
%
%   To fit custom models, create an anonymous function and use this as the FT argument.
%
%
%   [FO, GOF] = QUANTILEFIT(X, Y, ...) returns appropriate goodness-of-fit measures, for
%   the given inputs, in the structure GOF. GOF includes the fields:
%       -- Pseudo R2 for linear models
%
%
%   Examples:
%
% %   To fit a linear model through x and y, specifying a percentile:
%
%      [curve, goodness] = quantileFit( x, y, 'linear', tau );
%
%   To fit a polynomial curve of degree 2 in x for the 80 percentile, using custom model:
%
%      sf = quantileFit( x, y, @(a,b,c,x) ax^2 + bx + c, 0.8 );
%
%   Remarks on fitting:
%
%   The toolbox selects default initial values for coefficients uniformly at
%   random using MultiStart for 50 iterations. As a result, the best fitting
%   is kept at the end.
%
%   See also FUNCTION_HANDLE, MULTISTART, CREATEOPTIMPROBLEM.
%   v1.2, 24/09/2023    15:00
%   v1.3, Back to fmincon
%   v1.4, options varargin Bounds A FAIRE
%
%       A FAIRE: boostrap CI, goodness of fit indicators
%       voir chap6 de "uribe 2020, Quantile regression for cross sectional and
%       time series data"

%   Copyright 2023, @PabDawan, Fit a curve with quantile regression.


if nargin < 4
    error('Number of parameter inputs should be at least 4 (xData,yData,model,Tau).')

end

if ischar(funMod)
    switch funMod
        case 'linear'
            funMod = @(a,b,x) a*x+b;
        case '2-param'
            funMod = @(Dp,Vc,t) ((Dp./t)+Vc).*t;
        case '3-param'
            funMod = @(Vi,Vc,Tau,t) (Tau.*(Vi-Vc)./(t+Tau))+Vc;
        case 'power'
            funMod = @(S,E,t) S.*t.^(E-1);
        otherwise
            error('Check the spelling of your model. List of models can be found on the function help')
    end

end
%%
[xdata,ydata] = prepareCurveData(xdata,ydata);
nbParametre = nargin(funMod) - 1;
model = @(x)quantEval(x,xdata,ydata,tau,funMod);
if nargin < 5
    lb = [];
    ub = [];
end
problem = createOptimProblem('fmincon','objective', model, ...
                                'xdata', xdata, ...
                                'ydata', ydata, ...
                                'x0',ones(1,nbParametre), ...
                                'lb', lb,'ub', ub);


%%
ms = MultiStart;
ms.Display = 'off'; % 'final', 'iter', 'off'
% rng default;
[bestx,o.fval,o.exitflag,o.output,o.solutions] = run(ms, problem, 50);
o.bestParam = bestx;

%%
strFun = func2str(funMod);
strParam = extractBetween(strFun,'@(',')');
strParam = strsplit(string(strParam),',');

paramValues = num2cell(bestx);
fitObj = cell2struct(paramValues,cellstr(strParam(1:end-1)),2);
fittedFun = @(x) funMod(struct('struc', (paramValues)).struc,x);

fitObj.funFit = fittedFun;

%% Goodness of fit
% Uribe, J. M., & Guillen, M. (2020).
% Quantile Regression for Cross-Sectional and Time Series Data.
[~, residuals] = quantEval(bestx,xdata,ydata,tau,funMod);
% poids2 = 1 - (tau - ((tau - (residuals < 0))));
% poids = 1 - (residuals < 0);

% poids = residuals.*(tau - (residuals<0))
rho = (1-tau).*(residuals<0).*abs(residuals) + tau.*(residuals<0).*abs(residuals);
SSE = sum(rho.*((ydata) - fittedFun(xdata))); % Somme des carrés des erreurs / Vhat
SST = sum(rho.*(ydata - (xdata.*bestx(end))));
R2 = 1 - (SSE / SST); % Pseudo Coefficient de détermination
gof.pseudoRsquare = R2;
gof.rho = rho;

end

function [qef, residuals] = quantEval(x,xdata,ydata,tau,funMod)
% Quantile Error Function

funMod = funMod(struct('x', num2cell(x)).x,xdata);
residuals = ydata - funMod;
qef = sum(abs(residuals).*abs(tau-(residuals<0))); % cost function for quantile regression

end



