% ========================================================================================
%
% NAME: 	biexp_decay_fit.m
% PURPOSE: 	finds coeffs [x] of Cp=Dose*{x(1)*exp(-x(3)*(time-x(5)) + x(2)*exp(-x(4)*(time-x(5)))} 
%
% AUTHOR:   Hai-Ling Margaret Cheng
% DATE:     Sep 30, 2005
% MODIFIED: Oct 4, 2005
%
% INPUT:	x(1) and x(2): amplitudes of biexponential decay (kg/liter)
%		    x(3) and x(4): decay time constants (1/min)
%           x(5): time shift in AIF (min)
%		    xdata: time vector for Cp (min)
%           ydata: plasma concentration curve Cp (mmol/liter)
%           Dose: injection dose (mmol/kg)
%           wts: weights used in least-squares fitting
%
% OUTPUT:	err: difference between experimental and calculated data
%
% ========================================================================================

function err = biexp_decay_fit (x, xdata, ydata, Dose, wts)

%fit = Dose*(x(1)*exp(-x(3)*(xdata-x(5))) + x(2)*exp(-x(4)*(xdata-x(5))));
fit = Dose*(x(1)*exp(-x(3)*(xdata)) + x(2)*exp(-x(4)*(xdata)));
err = sqrt(wts).*(fit-ydata);

