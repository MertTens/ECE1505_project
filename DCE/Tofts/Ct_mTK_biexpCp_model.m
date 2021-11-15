% ========================================================================================
%
% NAME: 	Ct_mTK_biexpCp_model.m
% PURPOSE: 	calculates Ct from modified Tofts model assuming biexp Cp
%
% AUTHOR:   Hai-Ling Margaret Cheng
% DATE:     Jul 17, 2007
% MODIFIED: Jul 19, 2005
%
% INPUT:	t: time vector for Cp(t) and generated Ct(t) (min)
%		    Cp_coeff: [A1 A2 m1 m2 Dose] in units of kg/L, 1/min, mmol/kg
%		    respectively
%		    x: vector representing [Ktrans, vp, ve]
%
% OUTPUT:	Ct: tissue concentration (mM)
%
% NOTE: convolution is accurate only when time resolution is high.
% Therefore, in cases when temporal sampling is slow, accurate results must
% be obtained algebraicly. In the case of a biexponential AIF, there is a
% closed-form solution for Ct, as implemented in this module, which yields
% accurate Ct independent of the temporal resolution.
% ========================================================================================

function Ct = Ct_mTK_biexpCp_model (t, Cp_coeff, x)

A1 = Cp_coeff(1);
A2 = Cp_coeff(2);
m1 = Cp_coeff(3);
m2 = Cp_coeff(4);
Dose = Cp_coeff(5);

Ktrans = x(1);
vp = x(2);
ve = x(3);
tau = x(4);

% Calculate tissue concentration curve, Ct
Ct1 = vp*(A1*exp(-m1.*(t-tau)) + A2*exp(-m2.*(t-tau)));
Ct2 = Ktrans*(A1*((exp(-Ktrans.*(t-tau)/ve)-exp(-m1.*(t-tau)))/(m1-Ktrans/ve)) + ...
    A2*((exp(-Ktrans.*(t-tau)/ve)-exp(-m2.*(t-tau)))/(m2-Ktrans/ve)));
Ct = Dose*(Ct1 + Ct2);

idx = find(t<tau);
Ct(1:length(idx)) = 0;

