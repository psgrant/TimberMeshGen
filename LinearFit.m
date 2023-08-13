function [fitresult, gof] = LinearFit(xBot1, yBot1)
%CREATEFIT2(XBOT1,YBOT1)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : xBot1
%      Y Output: yBot1
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 11-Aug-2022 08:32:01


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( xBot1, yBot1 );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'yBot1 vs. xBot1', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'xBot1', 'Interpreter', 'none' );
% ylabel( 'yBot1', 'Interpreter', 'none' );
% grid on

