% ========================================================================================
%
% NAME: 	fit_biexpCp.m
% PURPOSE: 	fits biexponential decay model to plasma curve, Cp
%
% AUTHOR:   Hai-Ling Margaret Cheng
% DATE:     Nov.2, 2005
% MODIFIED: Jul.17, 2007
%
% INPUT:	time: time vector for Cp(t) (minutes)
%           Cp: plasma concentration curve (mmol/liter)
%		    x0: initial guess fit parameters [A1 A2 m1 m2 tlag]
%           lb: lower bounds on fit parameters
%           ub: upper bounds on fit parameters
%           options: used in curve-fitting
%           Dose: dose of contrast agent (mmol/kg)
%
% OUTPUT:	x:final estimates for fit parameters [A1 A2 m1 m2 tlag]
%		    resnorm: sum of chi-squared error of curve-fitting
%
% ========================================================================================

function [x,resnorm] = fit_biexpCp(time,Cp,x0,lb,ub,options,Dose);

warning_status = warning;
warning off;

t_max = find(Cp==max(Cp));
if length(find(time>10)) <= 1
%    lb(4) = 0.009; 
%    ub(4) = 0.011; 
    [x,resnorm] = lsqnonlin(@biexp_decay_fit, x0, lb, ub, options, time(t_max:length(Cp)), Cp(t_max:length(Cp)),Dose,ones(1,length(Cp)-t_max+1));
else
    t_slow = min(find(time>10));
    [x,resnorm] = lsqnonlin(@biexp_decay_fit, x0, lb, ub, options, time(t_slow:length(Cp)), Cp(t_slow:length(Cp)),Dose,ones(1,length(Cp)-t_slow+1));
    if x(4)<x(3) && x(4)>0
        ub(4) = x(4);
    elseif x(3)<x(4) && x(3)>0
        ub(3) = x(3); 
    elseif x(4)>0
        ub(4) = x(4);
    elseif x(3)>0
        ub(3) = x(3);
    end
    wts=ones(size(time(t_max:length(Cp))));
    N = length(find(time(t_max:length(Cp)) > 5))    % find how many points, N, were acquired with NEX=4, compared to NEX=1 for DCEMRI
    wts(length(wts)-N:length(wts)) = 4*wts(length(wts)-N:length(wts));  
    [x,resnorm] = lsqnonlin(@biexp_decay_fit, x, lb, ub, options, time(t_max:length(Cp)), Cp(t_max:length(Cp)),Dose,wts);
end

