% ========================================================================================
%
% NAME: 	extractAIF.m
% PURPOSE: 	main program to estimate AIF.
%
% AUTHOR:   Hai-Ling Margaret Cheng
% DATE:     April 10, 2005
% MODIFIED: April 18, 2005
%
% INPUT:	Cpref: actual plasma concentration curve (mmol/liter)
%           CtA: tissue concentration curve in tissue A (mmol/liter)
%           CtB: tissue concentration curve in tissue B (mmol/liter)
%		    t:  time vector for CtA(t) and CtB(t) (minutes).
%
% OUTPUT:	Cp_est: estimated plasma concentration curve (mmol/liter)
%
% ========================================================================================

function [Cp_est] = extractAIF;

del_t = 1/60;       % time resolution of generated Cp(t) and Ct(t) (minutes).
t = [0:del_t:7];    % time vector for generated Cp(t) and Ct(t) (minutes).
Cpref = generate_Cp('gamma variate', 0.1, t);

% ------------------------------------------------------------------------------------
% STEP 1 -- Estimate Ktrans and ve from CtA curve in tissue A, where vp is low.
%           The purpose here is to estimate Cp accurately for t>1.5 min.
%           Therefore, we obtain initial estimates for Ktrans and ve,
%           and then we refine these estimates by fitting the estimated Cp
%           to measured Cp at one or two time points after 1.5 min.
%           Even if vp is not negligible, the results won't be affected
%           since we're fitting the portion after first-pass.
% ------------------------------------------------------------------------------------
del_t = 1/60;       % time resolution of reference Cp(t) (minutes).
    
% Generate CtA, or tissue concentration curve, using one of two models 
% {'Tofts', 'Lee'}. Use the following literature values for different tissues:
%       Ktrans=0.08, vp=0.02, ve=0.14  (resting muscle)
%       Ktrans=0.20, vp=0.10, ve=0.28  (XXX)
MTT = 10/60;        % mean capillary transit time (min)
CtA = generate_Ct('Tofts', Cpref, 0.3, 0.2, 0.14, 0.01/MTT, [0:del_t:7]);

% Add random, Gaussian noise to the generated signals Cpref and Ctref
%SNR = Inf;
%[CtA,Cpref] = addnoise(CtA,Cpref,SNR);
   
% Downsample Cp and Ct according to sampling rate of measured Cp and Ct
del_ts = t(2)-t(1); % sampling time resolution of measured Cp and Ct
skip = round(del_ts/del_t); 
CtA = CtA(1:skip:length(CtA));
Cpref = Cpref(1:skip:length(Cpref));

% Increase temporal resolution of tissue curves. 
if del_ts > 1/60
    interp_factor = ceil(del_ts*60);
    CtA = interp(CtA,interp_factor);
    t = interp(t,interp_factor);
end 

% Estimate Ktrans & ve from tissue A's curve, excluding first-pass portion.  
% IN NOISE, I MUST SMOOTH Ctref BEFORE I TAKE DERIVATIVE 
[KtransA,veA,chi] = fitAIF(Cpref,CtA,0.14,t)

% Estimate first-pass portion of Cp through my method (Patlak following
% 1-st order differential solution. This will produce an estimate of vp.
[Cp_est,KtransA,vpA,veA] = fitAIF_vp(Cpref,CtA,KtransA,veA,t);
KtransA
vpA
veA
plot(t,Cpref,'r.',t,Cp_est,'k');