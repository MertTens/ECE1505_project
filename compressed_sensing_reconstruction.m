function [truth, recon, obsv,mse,est_cov] = compressed_sensing_reconstruction(spokes);

% Variables
recon_pts = 150;
% dim = 74;



% Add some paths
addpath('Read_data');
addpath('Reconstruction');
addpath('Simulation');
addpath('DCE');
addpath('DCE\AUC');
addpath('DCE\extractAIF');
addpath('DCE\Patlak');
addpath('DCE\Tofts');
addpath('qT1');
addpath('qT1\analysis');

% Get the data to reconstruct for hte Kalman filter
[truth, pre, obsv, est_cov, shifted_usmat,dim] = gen_recon_data(spokes);
truth = abs(ifft2(truth));

% Some useful matrices
dft2 = dftmtx(dim);
idft2 = conj(dftmtx(dim))/dim;
vec_dft2 = kron(dft2, dft2);
vec_idft2 = kron(idft2, idft2);


% Vectorize the images
vec_ground_truth = reshape(truth(:,:,1:recon_pts), dim*dim, []);
vec_observations = zeros(dim*dim, recon_pts+2);
vec_observations(:,2:recon_pts+1) = reshape(obsv(:,:,1:recon_pts), dim*dim, []);
vec_observations(:,1) = vec_dft2*vec_ground_truth(:,1);
vec_observations(:,recon_pts+2) = vec_dft2*vec_ground_truth(:,recon_pts);
vec_usmat = reshape(shifted_usmat(:,:,1:recon_pts), dim*dim, []);
vec_DCE_pre = reshape(truth(:,:,1), dim*dim, []);
vec_init = rand(dim*dim, recon_pts+2);


H = zeros(dim*dim, dim*dim, recon_pts+2);
H(:,:,1) = eye(dim*dim)*vec_dft2;
H(:,:,recon_pts + 2) = eye(dim*dim)*vec_dft2;
vec_init(:,1) = vec_ground_truth(:,1);
vec_init(:,recon_pts + 2) = vec_ground_truth(:,1);
% for n = 2:recon_pts+1
%     H(:,:,n)= diag(vec_usmat(:,n-1))*vec_dft2;
% %     H(:,:,n)= eye(dim*dim)*vec_dft2;
% %      vec_init(:,n) = vec_ground_truth(1);
% end

% Add the initialization and the end point to tie down the CS
% reconstruction
%size(H)
%size(vec_observations);
H1 = zeros(dim,dim,recon_pts + 2);
H1(:,:,1) = ones(dim,dim);
H1(:,:,recon_pts+2) = ones(dim,dim);
H1(:,:,2:recon_pts+1) = shifted_usmat(:,:,1:recon_pts);
% H1 = ones(dim,dim,recon_pts + 2);
% param.R = single(H);
param.y = reshape(vec_observations,dim,dim,[]);
param.lambda = 0.05*max(abs(vec_ground_truth(:,1)));
param.nuclear = 0.25*max(abs(vec_ground_truth(:,1)));
param.nite = 100;
param.display = 1;
param.time_pts = recon_pts + 2;
param.dim = dim;
param.H = H1;
param.E = dft2;

W = zeros(recon_pts + 2, recon_pts + 1);
for n = 1: recon_pts + 1
    W(n,n) = 1;
    W(n+1,n) = -1;
end
param.W = W;
vec_init = reshape(vec_init,dim,dim,[]);
clear W H vec_observations vec_usmat dft2 idft2
[reconstructed_data] = compressed_sensing(vec_init,  param);

for n = 1:2
    [reconstructed_data] = compressed_sensing(reconstructed_data,  param);
end
%smoothed = kalman_smoother(reconstructed_data, Ms, A, vec_ground_truth(:,recon_pts));
reconstructed_data = reconstructed_data(:,:,2:recon_pts+1);
reconstructed_data = abs(reconstructed_data);
recon = reconstructed_data;


mse = MSE(reshape(reconstructed_data,[],recon_pts), vec_ground_truth, 1, dim*dim);
figure()
plot(mse)
end