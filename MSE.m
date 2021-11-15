function [mse] = MSE(ground_truth, estimates, axis, n)
%
% This function will take some ground truth data and estimates, and return
% the mse

difference = abs(ground_truth - estimates) .^ 2;

sums = sum(difference, axis);

mse = sums;

mse = mse / n;

end