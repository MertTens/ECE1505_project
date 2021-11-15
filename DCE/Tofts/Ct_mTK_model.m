% ========================================================================================
%
% NAME: 	Ct_mTK_model.m
% PURPOSE: 	calculates Ct from the modified Tofts model
%
% AUTHOR:   Hai-Ling Margaret Cheng
% DATE:     March 15, 2005
% MODIFIED: Dec.16, 2005
%
% INPUT:	t: time vector for Cp(t) and generated Ct(t) (min)
%		    Cp: plasma concentration data (mM)
%		    x: vector representing [Ktrans, vp, ve]
%
% OUTPUT:	Ct: tissue concentration (mM)
%
% ========================================================================================

function Ct = Ct_mTK_model (t, Cp, x)

Ktrans = x(1);
vp = x(2);
ve = x(3);

del_t = t(2)-t(1);

% Implement convolution operation using FFTs for speed gain of 2x.
%     Ct = conv(Cp, exp(-Ktrans.*t/ve)) * del_t;
len = ceil(log(length(Cp)+length(t)) / log(2));
X = fft([Cp zeros(1,2^len-length(t))]);
Y = fft([exp(-Ktrans.*t/ve) zeros(1,2^len-length(Cp))]);
Ct = ifft(X.*Y)*del_t;

% Complete calculation of tissue concentration curve, Ct
Ct = vp*Cp + Ktrans*Ct(1:length(Cp));
