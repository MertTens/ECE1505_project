% ========================================================================================
%
% NAME: 	signal.m
% PURPOSE: 	calculates [Gd-DTPA] from signal equation for FSPGR
%
% AUTHOR:   	Hai-Ling Margaret Cheng
% DATE:     	February 26, 2004
% MODIFIED: 	February 26, 2004
%
% INPUT:	Ct: [Gd-DTPA] (mM)
%		K: gain constant in FSPGR sequence
%		theta: flip angle in FSPGR sequence (radians)
%		T10: native T1 of tissue pre-contrast (sec)
%		TR: repetition time in FSPGR sequence (sec)
%		TE: echo time in FSPGR sequence (sec)
%
% ========================================================================================

function X = signal (Ct, K, theta, T10, TR, TE)
  
% relaxation rates associated with T1 and T2
R1 = 4.5;       
R2 = 5.5;

A = K*exp(-TE*R2*Ct)*sin(theta);
B = 1-(exp(-TR*((1/T10)+(R1*Ct))));
C = 1-(cos(theta)*exp(-TR*((1/T10)+(R1*Ct))));
X =(A*B/C);
