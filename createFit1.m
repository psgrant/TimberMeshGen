function [fitresult, gof] = createFit1(ukxmax, kymax,ft)
%CREATEFIT1(UKXMAX,KYMAX)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : ukxmax
%      Y Output: kymax
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 10-Jun-2022 15:10:03


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( ukxmax, kymax );

% Set up fittype and options.
% ft = fittype( 'poly3' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );


