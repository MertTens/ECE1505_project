% ========================================================================================
%
% NAME: 	NLS_AIF_vp.m
% PURPOSE: 	finds [Ktrans,ve,chi] through least-squares minimization of 
%           Cp(t) based on Tofts model with vp=0.
%
% AUTHOR:   Hai-Ling Margaret Cheng
% DATE:     April 10, 2005
% MODIFIED: April 10, 2005
%
% INPUT:	vars_to_fit: variables to fit [Ktrans Ve Vp]
%		    Y: Ct/ve + (dCt/dt)/Ktrans
%		    Cp: tissue concentration data(mM)
%		    t: time vector for Cp(t) and generated Ct(t) (min)
%
% OUTPUT:	min: chi-squared difference
%
% ========================================================================================

function min = NLS_AIF_vp (vars_to_fit, Y, Cp, t)

Ktrans = vars_to_fit(1);
vp = vars_to_fit(2);
ve = vars_to_fit(3);

del_t = t(2)-t(1);

% find beginning of recirculation phase
idx = max(find(t<1.5));

% estimate Cp from Ct curve
b = Ktrans*(1/vp + 1/ve);
Cp_est = (Ktrans/vp) * conv(exp(-b*t), Y) * del_t;
Cp_est = Cp_est(1:length(t));

% non-linear least-squares minimization of Cp after t>1.5 min.
diff = Cp-Cp_est;
diff = diff(idx:length(diff));
min=1e3*sum(diff.^2)/length(diff);
