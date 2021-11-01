function [x, fval] = compressed_sensing(x0, param);
%
% Compressed sensing recon of undersampled k-space MRI data
%
% Based on program written by Ricardo Otazo
%
% L1-norm minimzation using non-linear conjugate gradient iterations

fprintf('\n Non-linear conjugate gradient algorithm')
fprintf('\n ---------------------------------------------\n')

% Starting point
x = x0;

% Line search parameters
% Modify if needed
maxlsiter = 150;
gradToll = 1e-3;
param.l1Smooth = 1e-15;
alpha = 0.01;
beta = 0.6;
t0 = 1;
k = 0;

% Initialize the gradient
g0 = grad(x, param);
dx = -g0;

%Iterate
while(1)
    
    % Backtracking line search
    f0 = objective(x,dx,0,param);
    t =  t0;
    f1 = objective(x,dx,t,param);
    lsiter = 0;
	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 & (lsiter<maxlsiter)
		lsiter = lsiter + 1;
		t = t * beta;
		f1 = objective(x,dx,t,param);
	end

	if lsiter == maxlsiter
		disp('Error - line search ...');
		return;
	end

	% control the number of line searches by adapting the initial step search
	if lsiter > 2, t0 = t0 * beta;end 
	if lsiter<1, t0 = t0 / beta; end

    % update x
	x = (x + t*dx);

	% print some numbers	
    if param.display,
        fprintf(' ite = %d, cost = %f \n',k,f1);
    end
    
    %conjugate gradient calculation
	g1 = grad(x,param);
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
	g0 = g1;
	dx =  - g1 + bk* dx;
	k = k + 1;
	
	% stopping criteria (to be improved)
	if (k > param.nite) || (norm(dx(:)) < gradToll), break;end
end
return

function res = objective(x,dx,t,param)

% L2 norm part
% w = pagemtimes(param.E,(reshape(x+t*dx,[],1,param.time_pts))) - reshape(param.y,[],1,param.time_pts);
w = param.H.*pagemtimes(pagemtimes(param.E,x + t*dx),param.E.') - param.y;
L2Obj = w(:)' * w(:);

% L1 norm part
if param.lambda
   w = reshape((x+t*dx),[],param.time_pts)*param.W; 
   L1Obj = sum((conj(w(:)).*w(:)+param.l1Smooth).^(1/2));
else
    L1Obj=0;
end

% nuclear part
if param.nuclear
    S = svd(reshape((x+t*dx), [], param.time_pts));
    L1Nuc = sum(abs(S));
else
    L1Nuc = 0;
end

% Objective function
res=L2Obj+param.lambda*L1Obj+param.nuclear*L1Nuc;

function g = grad(x,param)

% The L2 norm part
%  test1 = 2.*pagemtimes(pagectranspose(param.R), (pagemtimes(param.R,reshape(x,[],1,param.time_pts))-reshape(param.y,[],1,param.time_pts)));
%L2Grad1 = pagemtimes(param.R,reshape(x,[],1,param.time_pts))-reshape(param.y,[],1,param.time_pts);
% L2Grad = reshape(L2Grad,[],param.time_pts);
%L2Grad1 = pagemtimes(param.R,reshape(x,[],1,param.time_pts));
L2Grad = 2.*pagemtimes(pagemtimes((param.E'),(param.H).*((param.H).*(pagemtimes(pagemtimes(param.E,x),pagetranspose(param.E))) - param.y)),pagetranspose(param.E'));



% size(test)
% imagesc(x(:,:,1))
% figure()
% imagesc(abs(ifft2(test(:,:,1))))
% rand(2,2)*rand(3,3);

%L2Grad(:)'*L2Grad(:)
%L2Grad1(:)'*L2Grad1(:)

% The L1 norm part
if param.lambda
    w = reshape(x,[],param.time_pts)*param.W;
   L1Grad = (w.*(w.*conj(w)+param.l1Smooth).^(-0.5)) * param.W';
else
    L1grad = 0;
end

% nuclear part
if param.nuclear
    [U, S, V] = svd(reshape((x), [], param.time_pts));
    S(find(S ~= 0)) = 1;
    L1NucGrad = U*S*V';
else
    L1NucGrad = zeros(param.dim*param.dim,param.time_pts);
end

g = L2Grad + param.lambda*reshape(L1Grad,param.dim,param.dim,[]) + param.nuclear*reshape(L1NucGrad,param.dim,param.dim,[]);