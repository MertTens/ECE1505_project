% ========================================================================================
%
% NAME: 	fit_mTK_biexpCp_model.m
% PURPOSE: 	finds [Ktrans,vp,ve] through least-squares minimization of modified Tofts model
%
% AUTHOR:   Hai-Ling Margaret Cheng
% DATE:     July 17, 2007
% MODIFIED: July 18, 2007
%
% INPUT:	x: variables to fit [Ktrans Vp Ve]
%		    t: time vector for Cp(t) and generated Ct(t) (min)
%		    Ct: tissue concentration data(mM)
%		    Cp_coeff: [A1 A2 m1 m2 Dose] in units of kg/L, 1/min, mmol/kg
%
% OUTPUT:	err: difference between experimental and calculated data
%
% ========================================================================================

function err = fit_mTK_biexpCp_model (x, t, Ct, Cp_coeff)

tissue_model = Ct_mTK_biexpCp_model(t, Cp_coeff, x);
err = 1e2*(tissue_model - Ct);

I = find(t < 1);
idx = length(I);
err(1:idx) = 0;
