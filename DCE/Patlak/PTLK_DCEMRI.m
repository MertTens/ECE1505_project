% ========================================================================================
%
% NAME: 	PTLK_DCEMRI.m
% PURPOSE: 	main program for DCE-MRI analysis using the Patlak method
%
% AUTHOR:   Hai-Ling Margaret Cheng
% DATE:     March 15, 2005
% MODIFIED: March 16, 2007
%
% INPUT:	Cp: plasma concentration curve (mmol/liter)
%           Ct: tissue concentration curve (mmol/liter)
%		    t:  time vector for Cp(t) and Ct(t) (minutes).
%
% OUTPUT:	Ktrans: transfer constant (1/min)
%		    vp: plasma volume fraction
%           Ktrans_max: maximum transfer constant (1/min)
%		    vp_max: plasma volume fraction corresponding to Ktrans_max
%
% ========================================================================================

function [Ktrans,vp,Ktrans_max,vp_max] = PTLK_DCEMRI(Cp,Ct,t);

warning_status = warning;
warning off;

del_t = t(2)-t(1);      % time resolution (min)

% ------------------------------------------------------------------------------------
% Select window from which datapoints are used for linear fit 
% ------------------------------------------------------------------------------------
idx = find(t>=0 & t<=2);    % narrow window for large Ktrans, wider for small Ktrans

% ------------------------------------------------------------------------------------
% Transform the x- and y-coordinates
% ------------------------------------------------------------------------------------
X = (cumtrapz(Cp)*del_t) ./ Cp;
Y = Ct ./ Cp;

% ------------------------------------------------------------------------------------
% Intelligent fitting of straight-line to plot of Y-vs-X
% ------------------------------------------------------------------------------------
line = robustfit(X(idx),Y(idx));    % reduce sensitivity to outliers
m = line(2);
I = find(gradient(Y,X)>0 & gradient(Y,X)<1.2*m);
startidx = I(1);
endidx = idx(length(idx));

smoothY = smooth(Y,'rlowess');
diff1 = gradient(smoothY,X);
Ktrans_max = max(diff1(startidx:endidx));
b = find(Ktrans_max == diff1);
vp_max = Y(b) - X(b)*Ktrans_max;

line = robustfit(X(startidx:endidx),Y(startidx:endidx));
Ktrans = line(2);
vp = line(1);

plot(X,Y,'.'); hold on;
plot(X(idx), X(idx)*Ktrans+vp,'r-', X(idx), X(idx)*Ktrans_max+vp_max,'k-');

