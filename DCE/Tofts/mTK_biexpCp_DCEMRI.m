% ========================================================================================
%
% NAME: 	mTK_biexpCp_DCEMRI.m
% PURPOSE: 	main program for DCE-MRI analysis using the modified Tofts method
%
% AUTHOR:   Hai-Ling Margaret Cheng
% DATE:     July 17, 2007
% MODIFIED: July 22, 2007
%
% INPUT:	Cp_coeff: [A1 A2 m1 m2 Dose] in units of kg/L, 1/min, mmol/kg
%           Ct: tissue concentration curve (mmol/liter)
%		    t:  time vector for Cp(t) and Ct(t) (minutes).
%
% OUTPUT:	Ktrans: transfer constant (1/min)
%		    vp: plasma volume fraction
%		    ve: extravascular, extracellular volume fraction
%		    chi2: chi-squared error of curve-fitting
%
% ========================================================================================

function [Ktrans,vp,ve,tau,chi2] = mTK_biexpCp_DCEMRI(Cp_coeff,Ct,t);

warning_status = warning;
warning off;

% ------------------------------------------------------------------------------------
% Define options, bounds, and initial conditions in curve-fitting.
% -------------------------------------------------------------------------
lead_zeros = find(Ct<=0);
num_discard = length(lead_zeros);
t = t(num_discard:length(t));
Ct = Ct(num_discard:length(Ct));

% ------------------------------------------------------------------------------------
% Define options, bounds, and initial conditions in curve-fitting.
% ------------------------------------------------------------------------------------
options = optimset('LargeScale','on','MaxFunEvals', 1000,'Display','off');
lb=[0.001 0.001 0.01 0]; 
ub=[10 0.5 1 0.5];
Ktrans = 0.01; vp = 0.01; ve = 0.2; tau = 0.5;
x0=[Ktrans vp ve tau];

% ------------------------------------------------------------------------------------
% Levenberg-Marquardt non-linear least-squares regression. I discovered 
% 'lsqnonlin' is 3-4x faster than previous implementation using 'fmincon'.
% ------------------------------------------------------------------------------------
[x,resnorm,residue,flag] = lsqnonlin(@fit_mTK_biexpCp_model, x0, lb, ub, options, t, Ct, Cp_coeff);
Ktrans = x(1);
vp = x(2);
ve = x(3);
tau = x(4);
chi2 = resnorm;

margin = sum(abs(residue-Ct)./Ct)/length(Ct);
if margin > 0
    % Try different initial conditions for tau.
    x0=[Ktrans vp ve 0.1];         
    [x1,resnorm1] = lsqnonlin(@fit_mTK_biexpCp_model, x0, lb, ub, options, t, Ct, Cp_coeff);

    x0=[Ktrans vp ve 0.2];      
    [x2,resnorm2] = lsqnonlin(@fit_mTK_biexpCp_model, x0, lb, ub, options, t, Ct, Cp_coeff);

    x0=[Ktrans vp ve 0.3];      
    [x3,resnorm3] = lsqnonlin(@fit_mTK_biexpCp_model, x0, lb, ub, options, t, Ct, Cp_coeff);
    
    x0=[Ktrans vp ve 0.4];      
    [x4,resnorm4] = lsqnonlin(@fit_mTK_biexpCp_model, x0, lb, ub, options, t, Ct, Cp_coeff);
    
    [resnorm, I] = min([resnorm1 resnorm2 resnorm3 resnorm4 chi2]);
    if (I == 1),
	    x = x1; 
    elseif (I == 2),
	    x = x2; 
    elseif (I == 3),
	    x = x3; 
    elseif (I == 4),
	    x = x4; 
    end

    Ktrans = x(1);
    vp = x(2);
    ve = x(3);
    tau = x(4);
    chi2 = resnorm;
   
end


% Re-do fit the threshold 'margin' is exceeded. I have defined 'margin'
% such that if the mean difference between the fitted and experimental Ct
% exceeds 20% of the value of Ct, a re-fit is performed.
margin = sum(abs(residue-Ct)./Ct)/length(Ct);
if margin > 0.0
    % Try different initial conditions for Ktrans, Vp, Ve.
    x0=[0.01 0.2 0.2 tau];         
    [x1,resnorm1] = lsqnonlin(@fit_mTK_biexpCp_model, x0, lb, ub, options, t, Ct, Cp_coeff);

    x0=[0.5 0.01 0.2 tau];      
    [x2,resnorm2] = lsqnonlin(@fit_mTK_biexpCp_model, x0, lb, ub, options, t, Ct, Cp_coeff);

    x0=[0.5 0.2 0.2 tau];      
    [x3,resnorm3] = lsqnonlin(@fit_mTK_biexpCp_model, x0, lb, ub, options, t, Ct, Cp_coeff);
    
    [resnorm, I] = min([resnorm1 resnorm2 resnorm3 chi2]);
    if (I == 1),
	    x = x1; 
    elseif (I == 2),
	    x = x2; 
    elseif (I == 3),
	    x = x3; 
    end

    Ktrans = x(1);
    vp = x(2);
    ve = x(3);
    chi2 = resnorm;
   
end

%disp(['Ktrans=' num2str(Ktrans) '; vp=' num2str(vp) '; ve=' num2str(ve) '; tau=' num2str(tau)]);
%figure; plot(t,Ct,'k.'); hold on; plot(t,Ct_mTK_biexpCp_model (t, Cp_coeff, [Ktrans vp ve tau]));

warning (warning_status);
 
