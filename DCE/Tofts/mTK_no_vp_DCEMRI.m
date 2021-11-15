% ========================================================================================
%
% NAME: 	mTK_no_vp_DCEMRI.m
% PURPOSE: 	main program for DCE-MRI analysis using the modified Tofts method
%
% AUTHOR:   Hai-Ling Margaret Cheng
% DATE:     February 28, 2007
% MODIFIED: March 1, 2007
%
% INPUT:	Cp: plasma concentration curve (mmol/liter)
%           Ct: tissue concentration curve (mmol/liter)
%		    t:  time vector for Cp(t) and Ct(t) (minutes).
%
% OUTPUT:	Ktrans: transfer constant (1/min)
%		    ve: extravascular, extracellular volume fraction
%		    chi2: chi-squared error of curve-fitting
%
% ========================================================================================

function [Ktrans,ve,chi2] = mTK_no_vp_DCEMRI(Cp,Ct,t);

warning_status = warning;
warning off;

% ------------------------------------------------------------------------------------
% Define options, bounds, and initial conditions in curve-fitting.
% ------------------------------------------------------------------------------------
options = optimset('LargeScale','off','MaxFunEvals', 1000,'Display','off');
lb=[0.001 0.01]; 
ub=[10 1];
Ktrans = 0.01; ve = 0.2; 
x0=[Ktrans ve];

% ------------------------------------------------------------------------------------
% Levenberg-Marquardt non-linear least-squares regression. I discovered 
% 'lsqnonlin' is 3-4x faster than previous implementation using 'fmincon'.
% ------------------------------------------------------------------------------------
[x,resnorm,residue,flag] = lsqnonlin(@fit_mTK_no_vp_model, x0, lb, ub, options, t, Ct, Cp);
Ktrans = x(1);
ve = x(2);
chi2 = resnorm;

% Re-do fit the threshold 'margin' is exceeded. I have defined 'margin'
% such that if the mean difference between the fitted and experimental Ct
% exceeds 20% of the value of Ct, a re-fit is performed.
margin = sum(abs(residue-Ct)./Ct)/length(Ct);
if margin > 0.20
    % Try different initial conditions for Ktrans, Vp, Ve.
    x0=[0.01 0.2];         
    [x1,resnorm1] = lsqnonlin(@fit_mTK_no_vp_model, x0, lb, ub, options, t, Ct, Cp);

    x0=[0.5 0.2];      
    [x2,resnorm2] = lsqnonlin(@fit_mTK_no_vp_model, x0, lb, ub, options, t, Ct, Cp);
     
    [resnorm, I] = min([resnorm1 resnorm2 chi2]);
    if (I == 1),
	    x = x1; 
    elseif (I == 2),
	    x = x2; 
    end

    Ktrans = x(1);
    ve = x(2);
    chi2 = resnorm;
   
end

disp(['Ktrans=' num2str(Ktrans) '; ve=' num2str(ve)]);
figure; plot(t,Ct,'k.'); hold on; plot(t,Ct_mTK_no_vp_model (t, Cp, [Ktrans ve]));

warning (warning_status);
 
