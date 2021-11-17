function [A] = get_state_transition(truth);

dim = size(truth,1);
time_pts = size(truth,3);

A = zeros(dim,dim,time_pts-1);

for n = 2:time_pts
    A(n) = truth(n)./truth(n-1);
end
end
