% ========================================================================================
%
% NAME: 	NLS_AIF.m
% PURPOSE: 	finds [Ktrans,ve,chi] through least-squares minimization of 
%           Cp(t) based on Tofts model with vp=0.
%
% AUTHOR:   Hai-Ling Margaret Cheng
% DATE:     April 10, 2005
% MODIFIED: April 18, 2005
%
% INPUT:	vars_to_fit: variables to fit [Ktrans Ve Vp]
%		    Ct: tissue concentration data(mM)
%		    Cp: tissue concentration data(mM)
%		    t: time vector for Cp(t) and generated Ct(t) (min)
%
% OUTPUT:	min: chi-squared difference
%
% ========================================================================================

function min = NLS_AIF (vars_to_fit, Ct, Cp, t)

Ktrans = vars_to_fit(1);
ve = vars_to_fit(2);

% find beginning of recirculation phase
idx = max(find(t <= 1.5));

% estimate Cp from Ct curve
Cp_est = Ct/ve + (1/Ktrans)*gradient(Ct,t);

% non-linear least-squares minimization of Cp after t>1.5 min.
diff = Cp-Cp_est;
diff = diff(idx:length(diff));
min=1e3*sum(diff.^2)/length(diff);
