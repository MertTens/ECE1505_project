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

wah = x + t*dx;

% State model part
if param.stateL2
    stateL2Obj = 0;
    for n = 2:param.time_pts
        w = wah(:,:,n-1).*param.A(n-1) - wah(:,:,n);
        stateL2Obj = w(:)'*w(:);
    end
else
    stateL2Obj = 0;
end

% MLE part
if param.mle
    MLEObj = 0;
    for n = 2:param.time_pts
        w = (wah(:,:,n-1).*param.A(n-1) - wah(:,:,n)).*sqrt(param.C(:,:,n));
        MLEObj = w(:)'*w(:);
    end
else
    MLEObj = 0;
end

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
res=L2Obj+param.lambda*L1Obj+param.nuclear*L1Nuc+param.stateL2*stateL2Obj+MLEObj*param.mle;

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

% State model part
if param.stateL2
    stateL2Grad = zeros(param.dim,param.dim,param.time_pts);
    
    % First, middle last
    stateL2Grad(:,:,1) = (param.A(:,:,1).^2).*x(:,:,1)-param.A(:,:,1).*x(:,:,2);
    
    for n = 2:param.time_pts-1
        stateL2Grad(:,:,n) = (param.A(:,:,n).^2).*x(:,:,n)-param.A(:,:,n).*x(:,:,n+1) - param.A(:,:,n-1).*x(:,:,n-1) + x(:,:,n);
    end
    stateL2Grad(:,:,param.time_pts) = x(:,:,param.time_pts)-param.A(:,:,param.time_pts-1).*x(:,:,param.time_pts-1);
    stateL2Grad = stateL2Grad.*2;
else
    stateL2Grad = zeros(size(x));
end

% MLE part
if param.mle
    MLEGrad = zeros(param.dim,param.dim,param.time_pts);
    
    % First, middle last
    MLEGrad(:,:,1) = param.C(:,:,2).*(param.A(:,:,1).^2).*x(:,:,1)-param.C(:,:,2).*param.A(:,:,1).*x(:,:,2);
    
    for n = 2:param.time_pts-1
        MLEGrad(:,:,n) = param.C(:,:,n+1).*(param.A(:,:,n).^2).*x(:,:,n)-param.C(:,:,n+1).*param.A(:,:,n).*x(:,:,n+1) - param.C(:,:,n).*param.A(:,:,n-1).*x(:,:,n-1) + param.C(:,:,n).*x(:,:,n);
    end
    MLEGrad(:,:,param.time_pts) = param.C(:,:,param.time_pts).*(x(:,:,param.time_pts)-param.C(:,:,param.time_pts).*param.A(:,:,param.time_pts-1).*x(:,:,param.time_pts-1));
    MLEGrad = MLEGrad.*2;
else
    MLEGrad = zeros(size(x));
end

% nuclear part
if param.nuclear
    [U, S, V] = svd(reshape((x), [], param.time_pts));
    S(find(S ~= 0)) = 1;
    L1NucGrad = U*S*V';
else
    L1NucGrad = zeros(param.dim*param.dim,param.time_pts);
end

g = L2Grad + param.lambda*reshape(L1Grad,param.dim,param.dim,[]) + param.nuclear*reshape(L1NucGrad,param.dim,param.dim,[])+param.stateL2*stateL2Grad + param.mle*MLEGrad;