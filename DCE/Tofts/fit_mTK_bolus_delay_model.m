% ========================================================================================
%
% NAME: 	fit_mTK_bolus_delay_model.m
% PURPOSE: 	finds [Ktrans,vp,ve,tau] through least-squares minimization of modified Tofts model
%
% AUTHOR:   Hai-Ling Margaret Cheng
% DATE:     February 28, 2007
% MODIFIED: February 28, 2007
%
% INPUT:	x: variables to fit [Ktrans Vp Ve tar]
%		    t: time vector for Cp(t) and generated Ct(t) (min)
%		    Ct: tissue concentration data(mM)
%		    Cp: tissue concentration data(mM)
%
% OUTPUT:	err: difference between experimental and calculated data
%
% ========================================================================================

function err = fit_mTK_bolus_delay_model (x, t, Ct, Cp)

tissue_model = Ct_mTK_bolus_delay_model(t, Cp, x);
err = 1e2*(tissue_model - Ct);
err(find(isnan(err))) = 1000;

end