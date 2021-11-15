% ========================================================================================
%
% NAME: 	PTLK_DCEMRI.m
% PURPOSE: 	main program for DCE-MRI analysis using the Patlak method
%
% AUTHOR:   Hai-Ling Margaret Cheng
% DATE:     March 15, 2005
% MODIFIED: March 30, 2005
%
% INPUT:	Cp: plasma concentration curve (mmol/liter)
%           Ct: tissue concentration curve (mmol/liter)
%		    t:  time vector for Cp(t) and Ct(t) (minutes).
%
% OUTPUT:	Ktrans: transfer constant (1/min)
%		    vp: plasma volume fraction
%
% ========================================================================================

function [Ktrans,vp] = PTLK_DCEMRI_old(Cp,Ct,t);

warning_status = warning;
warning off;

del_t = t(2)-t(1);      % time resolution (min)
win = ceil(5/(del_t*60));        % size of window (5 sec) to look to compare slopes

% ------------------------------------------------------------------------------------
% Transform the x- and y-coordinates
% ------------------------------------------------------------------------------------
X = (cumtrapz(Cp)*del_t) ./ Cp;
Y = Ct ./ Cp;
% -------- Sort X in ascending order; re-arrange Y accordingly ---------
[X,I] = sort(X);
Y = Y(I);

% ------------------------------------------------------------------------------------
% Estimate Ktrans and vp
% ------------------------------------------------------------------------------------
% -------- Remove spurious values, esp.during onset of first-pass --------
Y = medfilt1(Y,5);

% -------- Ignore unreliable values during first-pass when both Cp and Ct are small --------
I = find(Cp==max(Cp));
J = find(Cp(1:I)<=0.1*max(Cp));
ss = max(J) + 1;
%for ss = 1:length(X)-win
%    r = corrcoef(X(ss:ss+win),Y(ss:ss+win));
%    if r(1,2)>0.9
%        break;
%    end
%end

% -------- Set up optimization linear function and settings --------
lb = []; ub = [];
x0 = [1 0];
options = optimset('MaxFunEvals', 1000,'Display','off');
straight_line = inline('x(1).*xdata + x(2)','x','xdata');

% -------- Compare slope of first 4 segments to classify uptake curve --------
stepsize=5*win;
ydata = Y(ss:ss+stepsize); xdata = X(ss:ss+stepsize);
[x,resnorm] = lsqcurvefit(straight_line, x0, xdata, ydata, lb, ub, options);
Ktrans1 = x(1);
ydata = Y(ss+stepsize:ss+2*stepsize); xdata = X(ss+stepsize:ss+2*stepsize);
[x,resnorm] = lsqcurvefit(straight_line, x0, xdata, ydata, lb, ub, options);
Ktrans2 = x(1); 
ydata = Y(ss+2*stepsize:ss+3*stepsize); xdata = X(ss+2*stepsize:ss+3*stepsize);
[x,resnorm] = lsqcurvefit(straight_line, x0, xdata, ydata, lb, ub, options);
Ktrans3 = x(1);
ydata = Y(ss+3*stepsize:ss+4*stepsize); xdata = X(ss+3*stepsize:ss+4*stepsize);
[x,resnorm] = lsqcurvefit(straight_line, x0, xdata, ydata, lb, ub, options);
Ktrans4 = x(1);
if (Ktrans1>Ktrans2) && (Ktrans2>Ktrans3) && (Ktrans3>Ktrans4)
    UPTAKE = 'rapid';
else
    UPTAKE = 'slow';
end
clear Ktrans1 Ktrans2 Ktrans3 Ktrans4

% -------- Find maximum slope over linear uptake portion --------
stepsize = win;
if (strcmp(UPTAKE,'rapid'))
    se = ss + 2*stepsize;
elseif (strcmp(UPTAKE,'slow'))
    se = ss + 10*stepsize;
end
ydata = Y(ss:se); xdata = X(ss:se);
[x,resnorm] = lsqcurvefit(straight_line, x0, xdata, ydata, lb, ub, options);
Ktrans = x(1);
slope_diff = 0;  
while (slope_diff < 1.05)
    se = se + stepsize;
    if (se > length(Y))
        break;
    end
    ydata = Y(ss:se); xdata = X(ss:se);
    [x,resnorm] = lsqcurvefit(straight_line, x0, xdata, ydata, lb, ub, options);
    slope_diff = Ktrans/x(1);      
end
se = se-stepsize;

ydata = Y(ss:se); xdata = X(ss:se);
[x,resnorm] = lsqcurvefit(straight_line, x0, xdata, ydata, lb, ub, options);
Ktrans = x(1);  % Ktrans = steepest initial slope on Y-vs-X curve    
vp = x(2);      % vp = y-intercept extrapolated from straight line fitted to steepest upslope --------
 
if vp < 0
    vp=0;
end

warning (warning_status);