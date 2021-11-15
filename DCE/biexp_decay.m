% ========================================================================================
%
% NAME: 	biexp_decay.m
% PURPOSE: 	finds coeffs [x] of Cp=Dose*{x(1)*exp(-x(3)*(time-x(5)) + x(2)*exp(-x(4)*(time-x(5)))} 
%
% AUTHOR:   Hai-Ling Margaret Cheng
% DATE:     Sep 30, 2005
% MODIFIED: Oct 4, 2005
%
% INPUT:	x(1) and x(2): amplitudes of biexponential decay (kg/liter)
%		    x(3) and x(4): decay time constants (1/min)
%           x(5): time shift in AIF (min)
%		    time: time vector for Cp (min)
%           Dose: injection dose (mmol/kg)
%
% OUTPUT:	Cp: plasma concentration curve (mmol/liter)
%
% ========================================================================================

function Cp = biexp_decay (x, time, Dose)

%Cp = Dose*(x(1)*exp(-x(3)*(time-x(5))) + x(2)*exp(-x(4)*(time-x(5))));
Cp = Dose*(x(1)*exp(-x(3)*(time)) + x(2)*exp(-x(4)*(time)));

