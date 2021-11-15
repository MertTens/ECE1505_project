% ========================================================================================
%
% NAME: 	fit_mTK_model.m
% PURPOSE: 	finds [Ktrans,vp,ve] through least-squares minimization of modified Tofts model
%
% AUTHOR:   Hai-Ling Margaret Cheng
% DATE:     March 15, 2005
% MODIFIED: Oct. 24, 2005
%
% INPUT:	x: variables to fit [Ktrans Vp Ve]
%		    t: time vector for Cp(t) and generated Ct(t) (min)
%		    Ct: tissue concentration data(mM)
%		    Cp: tissue concentration data(mM)
%
% OUTPUT:	err: difference between experimental and calculated data
%
% ========================================================================================

function err = fit_mTK_model (x, t, Ct, Cp)

tissue_model = Ct_mTK_model(t, Cp, x);
err = 1e2*(tissue_model - Ct);