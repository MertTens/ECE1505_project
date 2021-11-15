% ========================================================================================
%
% NAME: 	fitAIF.m
% PURPOSE: 	fitting measured Cp to that derived from Ct
%
% AUTHOR:   Hai-Ling Margaret Cheng
% DATE:     April 10, 2005
% MODIFIED: April 18, 2005
%
% INPUT:	Cp: measured plasma concentration curve (mmol/liter)
%           Ct: measured tissue concentration curve (mmol/liter)
%           ve: initial estimated extravascular, extracellular volume fraction
%		    t:  time vector for Cp(t) and Ct(t) (minutes).
%
% OUTPUT:	Ktrans: estimated transfer constant (1/min)
%           ve: estimated extravascular, extracellular volume fraction
%		    chi_err: chi-squared error of curve-fitting
%
% ========================================================================================

function [Ktrans,ve,chi_err] = fitAIF(Cp,Ct,ve,t);

warning_status = warning;
warning off;

% ------------------------------------------------------------------------------------
% Initialize curve-fitting parameters.
% ------------------------------------------------------------------------------------
Ktrans = 0.1;
vars_initial=[Ktrans ve];

% Options in curve-fitting
options = optimset('MaxFunEvals',1000,'Display','off','LevenbergMarquardt','on');  

% Bounds on parameters estimated in curve fitting
lb=[0.0001 0.0001]; 
ub=[10 1];

% ------------------------------------------------------------------------------------
% Levenberg-Marquardt non-linear least-squares regression
% ------------------------------------------------------------------------------------
[vars_final,fval] = fmincon('NLS_AIF', vars_initial, [],[],[],[], ...
        lb, ub, [], options, Ct, Cp, t);

Ktrans = vars_final(1);
ve = vars_final(2);
chi_err = fval;
