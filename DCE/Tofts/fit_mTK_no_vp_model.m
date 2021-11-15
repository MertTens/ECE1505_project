% ========================================================================================
%
% NAME: 	fit_mTK_no_vp_model.m
% PURPOSE: 	finds [Ktrans,ve] through least-squares minimization of modified Tofts model
%
% AUTHOR:   Hai-Ling Margaret Cheng
% DATE:     March 15, 2005
% MODIFIED: March 1, 2007
%
% INPUT:	x: variables to fit [Ktrans Ve]
%		    t: time vector for Cp(t) and generated Ct(t) (min)
%		    Ct: tissue concentration data(mM)
%		    Cp: tissue concentration data(mM)
%
% OUTPUT:	err: difference between experimental and calculated data
%
% ========================================================================================

function err = fit_mTK_no_vp_model (x, t, Ct, Cp)

tissue_model = Ct_mTK_no_vp_model(t, Cp, x);
err = 1e2*(tissue_model - Ct);