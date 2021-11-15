% ========================================================================================
%
% NAME: 	AUC_DCEMRI.m
% PURPOSE: 	main program for DCE-MRI analysis using the AUC method
%
% AUTHOR:   Hai-Ling Margaret Cheng
% DATE:     Nov.29, 2005
% MODIFIED: Dec.10, 2013
%
% INPUT:	Ct: tissue concentration curve (mmol/liter)
%		    t:  time vector for Cp(t) and Ct(t) (minutes).
%           tstart: first time point to start integration
%           interval: first X minutes over which AUC is computed
%
% OUTPUT:	AUC: integrated area over the first (interval) minutes
%
% ========================================================================================

function [AUC] = AUC_DCEMRI(Ct,t,tstart,interval);

warning_status = warning;
warning off;

% ------------------------------------------------------------------------------------
% Calculate AUCs -- no normalization against reference tissue.
% ------------------------------------------------------------------------------------
idx = find(t>=tstart & t<=(tstart+interval));
AUC = trapz(Ct(idx))*(t(2)-t(1));

%AUC = sum(Ct(idx));


if AUC < 0
    AUC=0;
end

