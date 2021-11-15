% ========================================================================================
%
% NAME: 	fitAIF_vp.m
% PURPOSE: 	fitting measured Cp to that derived from Ct
%
% AUTHOR:   Hai-Ling Margaret Cheng
% DATE:     April 10, 2005
% MODIFIED: April 20, 2005
%
% INPUT:	Cp: measured plasma concentration curve (mmol/liter)
%           CtA: measured tissue concentration curve (mmol/liter)
%           KtransA: initial estimate of transfer constant (1/min)
%           veA: initial estimate of extravascular, extracellular volume fraction
%		    t:  time vector (minutes).
%
% OUTPUT:	Cp_est: estimated plasma concentration curve (mmol/liter)
%           vpA: estimated plasma volume fraction
%
% ========================================================================================

function [Cp_est,KtransA,vpA,veA] = fitAIF_vp(Cp,CtA,KtransA,veA,t);

warning_status = warning;
warning off;

del_t = t(2)-t(1);
t90 = max(find(t<1.5));

% ------------------------------------------------------------------------------------
% Estimate vp and first-pass portion of Cp curve.
% ------------------------------------------------------------------------------------
vp = [0.005:0.005:0.2]; 
delKtrans = [];  delCp = [];
for idx=1:length(vp) 
    % estimate Cp using parameters from tissue A (my method)
    YA = CtA/veA + (1/KtransA)*gradient(CtA,t);
    b = KtransA*(1/vp(idx) + 1/veA);
    Cp_est = (KtransA/vp(idx)) * conv(exp(-b*t), YA) * del_t;
    Cp_est = Cp_est(1:length(t));
 
    % refine Cp estimate by Patlak's method
    % Patlak's vp will correspond to vpA, regardless of the accuracy of vpA.
    % However, Patlak's Ktrans will correspond to the true Ktrans only if
    % the vpA used to generate Cp_est is the true vp.
    [KtransPTLK,vpPTLK] = PTLK_DCEMRI(Cp_est,CtA,t);
    delKtrans = [delKtrans (KtransPTLK-KtransA)];
    
    chi2 = sum((Cp_est(t90:length(t))-Cp(t90:length(t))).^2) / (length(t)-t90);
    delCp = [delCp chi2];
end
range = find(abs(delCp) < 1e-2);
idx = find(abs(delKtrans(range)) == min(abs(delKtrans(range))));
vpA = vp(range(idx));

% ------------------------------------------------------------------------------------
% The permeability portion of Cp_est may not correspond to measured Cp 
% due to inaccurate ve. Thus, refine estimate of ve by fitting permeability
% phase of Cp also.
% ------------------------------------------------------------------------------------
%ve = [0.01:0.01:0.7];  delCp=[];
%for idx=1:length(ve)
%    YA = CtA/ve(idx) + (1/KtransA)*gradient(CtA,t);
%    b = KtransA*(1/vpA + 1/ve(idx));
%    Cp_est = (KtransA/vpA) * conv(exp(-b*t), YA) * del_t;
%    Cp_est = Cp_est(1:length(t));
%    delCp = [delCp sum((Cp_est - Cp).^2)/length(Cp)];
%end
%veA = ve(find(delCp==min(delCp)));

% ------------------------------------------------------------------------------------
% Refine estimate of Ktrans, using new estimate of vpA.
% ------------------------------------------------------------------------------------
%delKtrans = [];  
%Ktrans = [(KtransA*0.9):0.001:(KtransA*1.1)];
%for idx=1:length(Ktrans) 
    % estimate Cp using parameters from tissue A (my method)
%    YA = CtA/veA + (1/Ktrans(idx))*gradient(CtA,t);
%    b = Ktrans(idx)*(1/vpA + 1/veA);
%    Cp_est = (Ktrans(idx)/vpA) * conv(exp(-b*t), YA) * del_t;
%    Cp_est = Cp_est(1:length(t));
 
    % refine Cp estimate by Patlak's method
    % Patlak's vp will correspond to vpA, regardless of the accuracy of vpA.
    % However, Patlak's Ktrans will correspond to the true Ktrans only if
    % the vpA used to generate Cp_est is the true vp.
%    [KtransPTLK,vpPTLK] = PTLK_DCEMRI(Cp_est,CtA,t);
  
%    delKtrans = [delKtrans (KtransPTLK-Ktrans(idx))];
%end
%KtransA = Ktrans(find(abs(delKtrans)==min(abs(delKtrans))));

% ------------------------------------------------------------------------------------
% Calculate Cp_est using estimates of Ktrans, vp, and ve.
% ------------------------------------------------------------------------------------
YA = CtA/veA + (1/KtransA)*gradient(CtA,t);
b = KtransA*(1/vpA + 1/veA);
Cp_est = (KtransA/vpA) * conv(exp(-b*t), YA) * del_t;
Cp_est = Cp_est(1:length(t));
   
% ------------------------------------------------------------------------------------
% Refine all estimates using Tofts model for t > 0.
% ------------------------------------------------------------------------------------
vars_initial=[KtransA vpA veA];
options = optimset('MaxFunEvals',1000,'Display','off');  
lb=[0.0001 0.0001 0.0001];   ub=[10 1 1];
[vars_final,fval] = fmincon('fit_mTK_model', vars_initial, [],[],[],[], ...
    lb, ub, [], options, CtA, Cp_est, t);
KtransA = vars_final(1);
vpA = vars_final(2);
veA = vars_final(3);

% ------------------------------------------------------------------------------------
% Re-calculate Cp_est using Tofts' estimates of Ktrans, vp, and ve.
% ------------------------------------------------------------------------------------
YA = CtA/veA + (1/KtransA)*gradient(CtA,t);
b = KtransA*(1/vpA + 1/veA);
Cp_est = (KtransA/vpA) * conv(exp(-b*t), YA) * del_t;
Cp_est = Cp_est(1:length(t));