% ========================================================================================
%
% NAME: 	unitstep.m
% PURPOSE: 	function for unit-step response.
%
% AUTHOR:   	Hai-Ling Margaret Cheng
% DATE:     	March 14, 2005
% MODIFIED: 	March 14, 2005
%
% INPUT:	t: rational number
% OUTPUT:   y=0 for t<0; 1 for t>=0
%
% ========================================================================================

function y = unitstep (t)

y = zeros(size(t));
I = find(t >= 0);
y(I) = 1;

return