function [truth, pre, us_observations, est_cov, shifted_usmat,dim] = gen_recon_data(spokes);

%
% TODO: MAKE THE STATE TRANS ESTIMATION USE NOT THE UNIT STEP AIF
%
% Some variables to define
measurements = 1000;
cov_group_sizes = [250,500];
cov_group_size = 1;
cov_groups = measurements / cov_group_size;
recon_group_size = 20;
recon_groups = 100;
dt = 3.7;
start_angle = 15;
calc_gt_ktrans = 1;
calc_est_ktrans = 1;

% According to the variance of a noise region of the image, the noise is
% really big
% I'll make it less than what it says it is for now, but apparently its
% around 5e7 - 9e7
measurement_noise = 1e5;
process_noise = 100000;
TR = 6.1962;
FA = 0.3491;
sliceloc = 7;

[SI, pre, kt, ve, vp, t1,dim] = sim_si_images(measurements);

% Some useful matrices
dft2 = dftmtx(dim);
idft2 = conj(dftmtx(dim))/dim;
vec_dft2 = kron(dft2, dft2);
vec_idft2 = kron(idft2, idft2);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACQUIRE TIME SERIES DATA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pre = SI(:,:,1);
% cov = time_series_covariance(DCE_post);
ground_truth = SI;

k_space = fft(permute(fft(permute(ground_truth, [2,1,3])), [2,1,3]));
init_k_space = fft2(SI);


% Add some measurement noise (simulation only)
% Undersample it to get the observations
% Chagne the undersampling on every frame
GA = pi*(3-sqrt(5));
usmat = zeros(dim, dim, measurements);
shifted_usmat = zeros(dim, dim, measurements);

ang = start_angle;
for n = 1:measurements
    [usmat_single, shifted_usmat_single] = radial_undersampling(dim, spokes, ang, GA);
    ang = ang + (spokes+1)*GA;
    usmat(:,:,n) = usmat_single;
    shifted_usmat(:,:,n) = shifted_usmat_single;
end

a_factor = dim*dim / sum(sum(usmat_single));
fprintf("The acceleration factor is %f\n", a_factor)

observations = k_space.* shifted_usmat;

% true_cov = time_series_covariance(SI);

% Combine the images to get a better estimate of covariance
% est_cov = time_series_covariance(abs(ifft2(observations)));


% est_cov_imgs(:,:,cov_groups+1) = pre;
num_groups = size(cov_group_sizes,2);
est_cov = zeros(dim*dim,dim*dim);
for a = 1:num_groups
    cov_group_size = cov_group_sizes(a);
    cov_groups = measurements / cov_group_size;
    est_cov_imgs = zeros(dim, dim, cov_groups);
    for g = 1:cov_groups
        fftimg = zeros(dim, dim);
        for n = 1:cov_group_size
            idx = (g-1)*cov_group_size + n;
            obsv = observations(:,:,idx);
            fftimg(find(obsv ~= 0)) = obsv(find(obsv ~= 0));
            
        end
        est_cov_imgs(:,:,g) = abs(ifft2(fftimg));
%         figure()
%         imagesc(log(abs(fftimg)))
    end
    est_cov = est_cov + time_series_covariance(est_cov_imgs);

end
est_cov = est_cov ./ num_groups;

us_observations = observations;
truth = k_space;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VECTORIZE THE IMAGES AND GENERATE THE FT MATRICES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now, I will vectorize the images so they can be used by the kalman filter

vec_k_space = reshape(k_space, dim*dim, []);
vec_observations = reshape(observations, dim*dim, []);
vec_ground_truth = reshape(ground_truth, dim*dim, []);
vec_usmat = reshape(shifted_usmat, dim*dim, []);
vec_DCE_pre = reshape(DCE_post(:,:,1), dim*dim, []);
%vec_phase_map = reshape(phase_map(:,:,1), img_dim*img_dim, []);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALL THE KALMAN FILTER AND GET MSE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First, get some matrices ready. Then call it.


A = eye(dim*dim);
%A = state_trans;
B = eye(dim*dim);
% Cu = eye(size(covariance, 1)) .* process_noise;
Cu = cov;
%H = diag(vec_usmat)*vec_dft2*diag(vec_phase_map);
% return
H = zeros(dim*dim, dim*dim, measurements);
% H = diag(vec_usmat)*vec_dft2;
for n = 1:measurements
    H(:,:,n)= diag(vec_usmat(:,n))*vec_dft2;
end
Cw = eye(size(state_trans, 1)).*measurement_noise;
%ADD NOISE TO OBSERCVATIONS
%vec_observations = vec_noisy_us_k_space;
recon_init = vec_DCE_pre;

error_init = eye(size(state_trans, 1)) .* measurement_noise;
error_init = error_init;

reconstructed_data = kalman_filter(A, B, Cu, H, Cw, vec_observations, vec_ground_truth, recon_init, error_init);
% [reconstructed_data, Mnn] = GPU_kalman_filter(A, B, Cu, H, Cw, vec_observations, vec_ground_truth, recon_init, error_init);
% reconstructed_data = global_difference_kalman_filter(A, B, Cu, H, Cw, vec_observations, vec_ground_truth, recon_init, error_init);
% reconstructed_data = GPU_ideal_global_difference_kalman_filter(A, B, Cu, H, Cw, vec_observations, vec_ground_truth, recon_init, error_init);

reconstructed_data = abs(reconstructed_data);
recon_init = vec_ground_truth(:,measurements);
% error_init = Mnn;

cov = time_series_covariance(flipped_post);
flipped_recon = flip(reconstructed_data, 2);

[reconstructed_data_second_pass, Mnn] = GPU_kalman_filter(A, B, Cu, A, Cw, flipped_recon, vec_ground_truth, recon_init, error_init);
reconstructed_data_second_pass = abs(reconstructed_data_second_pass);
recon_init = vec_DCE_pre;
error_init = Mnn;

% cov = time_series_covariance(DCE_post);
% flipped_recon = flip(reconstructed_data_second_pass, 2);
% flipped_recon(:,1) = reconstructed_data(:,1);
% 
% [reconstructed_data_third_pass, Mnn] = GPU_kalman_filter(A, B, Cu, A, Cw, flipped_recon, vec_ground_truth, recon_init, error_init);
% reconstructed_data_third_pass = abs(reconstructed_data_third_pass);

% clims = [0, max(max(vec_ground_truth))];
for n = 1:measurements
    img = reshape(reconstructed_data(:,n), dim, dim);
    recon_series(:,:,n+1) = img;
    temp = flip(reconstructed_data_second_pass, 2);
    img2 = reshape(temp(:,n), dim, dim);
    recon_series2(:,:,n+1) = img2;
%     figure()
%     imagesc(abs(img), clims)
%     title('reecon')
%     
%     img = reshape(vec_observations(:,n), dim, dim);
%     figure()
%     imagesc(abs(ifft2(img)), clims)
%     %    imagesc(abs(ifft2(img)),'InitialMagnification','fit')
% 
%     title('observed')
%     
%     img = reshape(vec_ground_truth(:,n), dim, dim);
%     figure()
%     imagesc(abs(img), clims)
%     %    imagesc(abs(ifft2(img)),'InitialMagnification','fit')
% 
%     title('truth')
end
vec_img_obsv = zeros(size(vec_observations));
for i = 1:size(vec_img_obsv, 2)
    vec_img_obsv(:,i) = abs(vec_idft2*vec_observations(:,i));
end
recon_mse = MSE(reconstructed_data, vec_ground_truth, 1, dim*dim);
observe_mse = MSE(abs(vec_img_obsv), vec_ground_truth, 1, dim*dim);
test_mse = MSE(vec_ground_truth, vec_ground_truth, 1, dim*dim);
figure()
plot(recon_mse)
title('recon')
figure()
plot(observe_mse)
title('observed')
% figure()
% plot(test_mse)
% title('test')

if(calc_est_ktrans)
    [k_trans_est, ve_est, vp_est] = mainDCEfromMAT(T1_pre, T1_pre, recon_series, TR, FA, sliceloc);
    [k_trans_est2, ve_est2, vp_est2] = mainDCEfromMAT(T1_pre, T1_pre, recon_series2, TR, FA, sliceloc);
end


diff = abs(k_trans_est - k_trans_truth);
diff2 = abs(k_trans_est2 - k_trans_truth);
mse = sum(sum(diff.^2)) / dim*dim;
mse2 = sum(sum(diff2.^2)) / dim*dim;

% Display Purposes
% for i = 1:time_pts
%     figure()
%     imagesc(log(abs(observations(:,:,i))));
% end
% 
% 
% figure()
% imagesc(DCE_pre);
% figure()
% imagesc(DCE_settled);
% figure()
% imagesc(T1_pre);
% figure()
% imagesc(T1_post);


% Plot some contrast time curves
% x = 40;
% y = 76;
% Ct = reshape(DCE_post(x,y,:), 1, []);
% 
% plot(Ct)
end