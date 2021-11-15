% ========================================================================================
%
% NAME: 	Ct_mTK_bolus_delay_model.m
% PURPOSE: 	calculates Ct from the modified Tofts model
%
% AUTHOR:   Hai-Ling Margaret Cheng
% DATE:     February 28, 2007
% MODIFIED: July 5, 2007
%
% INPUT:	t: time vector for Cp(t) and generated Ct(t) (min)
%		    Cp: plasma concentration data (mM)
%		    x: vector representing [Ktrans, vp, ve, tau]
%
% OUTPUT:	Ct: tissue concentration (mM)
%
% ========================================================================================

function Ct = Ct_mTK_bolus_delay_model (t, Cp, x)

Ktrans = x(1);
vp = x(2);
ve = x(3);
tau = x(4);

del_t = t(2)-t(1);
shift = round(tau/del_t);

% Implement convolution operation using FFTs for speed gain of 2x.
%     Ct = conv(Cp, exp(-Ktrans.*t/ve)) * del_t;
len = ceil(log(length(Cp)+length(t)) / log(2));
X = fft([Cp zeros(1,2^len-length(t))]);
Y = fft([exp(-Ktrans.*t/ve) zeros(1,2^len-length(Cp))]);
Ct = ifft(X.*Y)*del_t;

% Complete calculation of tissue concentration curve, Ct
Ct = vp*Cp + Ktrans*Ct(1:length(t));
Ct = [zeros(1,shift) Ct(1:length(t)-shift)];