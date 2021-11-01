function [SI, pre, kt, ve, vp, t1, dim] = sim_si_images(num_AIF_pts);

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

T1_pre_full = load('DCE_test_mat_T1_pre_contrast\T1_pre.mat').t1_map;
T1_post_full = load('DCE_test_mat_T1_post_contrast\T1_post.mat').t1_map;
DCE_pre_full = load('DCE_test_mat_pre_contrast\DCE_pre.mat').dicom_imgs;
DCE_settled_full = load('DCE_test_mat_settled_contrast\DCE_settled.mat').dicom_imgs;
DCE_post_full = load('DCE_test_mat_post_contrast\DCE_post.mat').dicom_imgs;


% Resize each image
dim = 74;
skip = 14;
% og_dim = size(DCE_pre_full, 1);
% mid = og_dim / 2;
% midrow = 80;
% midcol = 110;
% lowerrow = 1 + midrow - dim/2;
% upperrow = midrow + dim/2;
% lowercol = 1 + midcol - dim/2;
% uppercol = midcol + dim/2;
% 
% T1_pre = T1_pre_full(lowerrow:upperrow, lowercol:uppercol);
% T1_post = T1_post_full(lowerrow:upperrow, lowercol:uppercol);
% DCE_pre = DCE_pre_full(lowerrow:upperrow, lowercol:uppercol);
% DCE_post = DCE_post_full(lowerrow:upperrow, lowercol:uppercol,:);
% DCE_settled = DCE_settled_full(lowerrow:upperrow, lowercol:uppercol);
% k_trans_truth = load('ground_truth_ktrans\k_trans_truth.mat').k_trans_truth;
% ve_truth = load('ground_truth_ktrans\ve_truth.mat').ve_truth;
% vp_truth = load('ground_truth_ktrans\vp_truth.mat').vp_truth;




% T1_pre(find(T1_pre == 10000)) = 2000;
% k_trans_truth(find(k_trans_truth > 8)) = 1;

% Legend for tissue types
% nothing_p = 0;
% muscle_p = 255;
% tumor_p = 127;
% grey_matter_p = 195;

nothing = 4;
mat1 = 30;
mat2 = 51;
mat3 = 70;
mat4 = 100;
mat5 = 190;
mat6 = 256;


mat1_kt = 0.001;
mat2_kt = 10;
mat3_kt = 0.001;
mat4_kt = 10;
mat5_kt = 0.7;
mat6_kt = 0.4;


mat1_ve = 0.01;
mat2_ve = 0.01;
mat3_ve = 1;
mat4_ve = 1;
mat5_ve = 0.3;
mat6_ve = 0.7;


mat1_vp = 0.05;
mat2_vp = 0.1;
mat3_vp = 0.01;
mat4_vp = 0.3;
mat5_vp = 0.4;
mat6_vp = 0.04;


mat1_t1 = 1000;
mat2_t1 = 1000;
mat3_t1 = 1000;
mat4_t1 = 1000;
mat5_t1 = 1000;
mat6_t1 = 1000;


nothing_kt = 0;
nothing_ve = 0;
nothing_vp = 0;
nothing_t1 = 0;

% tumor_kt = 0.51;
tumor_kt = 0.7;
% tumor_ve = 0.25;
tumor_ve = 0.4;
% tumor_vp = 0.028;
tumor_vp = 0.05;
% tumor_t1 = 1875;
tumor_t1 = 1100;



muscle_kt = 0.1;
muscle_ve = 0.1;
muscle_vp = 0.01;
muscle_t1 = 1000;

grey_matter_kt = 0.02;
grey_matter_ve = 0.000001;
grey_matter_vp = 0.01;
grey_matter_t1 = 1154;

% t is in minutes here
% num_AIF_pts = 1000;
AIF_duration = 12;
t = linspace(0,AIF_duration,num_AIF_pts);
 

% SI = zeros(dim, dim, num_AIF_pts);
 SI = 0; 
 
% Take in 1000x1000 pixel png
clr = imread('1505anatomy.png');
img = rgb2gray(clr);
% img = imresize(img, dim/1000);
img = double(img);

large_k_trans = zeros(size(img));
large_ve = zeros(size(img));
large_vp = zeros(size(img));
large_t1 = zeros(size(img));


large_k_trans(find((img >= nothing) & (img < mat1))) = mat1_kt;
large_k_trans(find((img >= mat1) & (img < mat2))) = mat2_kt;
large_k_trans(find((img >= mat2) & (img < mat3))) = mat3_kt;
large_k_trans(find((img >= mat3) & (img < mat4))) = mat4_kt;
large_k_trans(find((img >= mat4) & (img < mat5))) = mat5_kt;
large_k_trans(find((img >= mat5) & (img < mat6))) = mat6_kt;

large_ve(find((img >= nothing) & (img < mat1))) = mat1_ve;
large_ve(find((img >= mat1) & (img < mat2))) = mat2_ve;
large_ve(find((img >= mat2) & (img < mat3))) = mat3_ve;
large_ve(find((img >= mat3) & (img < mat4))) = mat4_ve;
large_ve(find((img >= mat4) & (img < mat5))) = mat5_ve;
large_ve(find((img >= mat5) & (img < mat6))) = mat6_ve;

large_vp(find((img >= nothing) & (img < mat1))) = mat1_vp;
large_vp(find((img >= mat1) & (img < mat2))) = mat2_vp;
large_vp(find((img >= mat2) & (img < mat3))) = mat3_vp;
large_vp(find((img >= mat3) & (img < mat4))) = mat4_vp;
large_vp(find((img >= mat4) & (img < mat5))) = mat5_vp;
large_vp(find((img >= mat5) & (img < mat6))) = mat6_vp;

large_t1(find((img >= nothing) & (img < mat1))) = mat1_t1;
large_t1(find((img >= mat1) & (img < mat2))) = mat2_t1;
large_t1(find((img >= mat2) & (img < mat3))) = mat3_t1;
large_t1(find((img >= mat3) & (img < mat4))) = mat4_t1;
large_t1(find((img >= mat4) & (img < mat5))) = mat5_t1;
large_t1(find((img >= mat5) & (img < mat6))) = mat6_t1;




% large_k_trans(find(img == muscle_p)) = muscle_kt;
% large_k_trans(find(img == tumor_p)) = tumor_kt;
% large_k_trans(find(img == grey_matter_p)) = grey_matter_kt;
% 
% large_ve(find(img == muscle_p)) = muscle_ve;
% large_ve(find(img == tumor_p)) = tumor_ve;
% large_ve(find(img == grey_matter_p)) = grey_matter_ve;
% 
% large_vp(find(img == muscle_p)) = muscle_vp;
% large_vp(find(img == tumor_p)) = tumor_vp;
% large_vp(find(img == grey_matter_p)) = grey_matter_vp;
% 
% large_t1(find(img == muscle_p)) = muscle_t1;
% large_t1(find(img == tumor_p)) = tumor_t1;
% large_t1(find(img == grey_matter_p)) = grey_matter_t1;

% maxval = max(max(img));
% img = img ./maxval;
% kt = imresize(large_k_trans, dim/1000);
% ve = imresize(large_ve, dim/1000);
% vp = imresize(large_vp, dim/1000);
% t1 = imresize(large_t1, dim/1000);

kt = large_k_trans(1:skip:end,1:skip:end);
ve = large_ve(1:skip:end,1:skip:end);
vp = large_vp(1:skip:end,1:skip:end);
t1 = large_t1(1:skip:end,1:skip:end);
dim = size(t1, 1);

% imagesc(kt)

% First, generate Cp
AIF = 'biexponential';
peak_scale = 1;
Dose = 0.05;             % Gd dose (mmol/kg) 
r1 = 4.1/60;               % T1 relaxivity (liter/mmol/sec) -- Magnevist
TR = 6.1962;
FA = 0.3491;
Cp = generate_Cp(AIF, Dose, peak_scale, t);

CtMAT = zeros(dim, dim, num_AIF_pts);
for row = 1:dim
    for col = 1:dim
        x = zeros(4,1);
        x(1) = kt(row,col);
        x(2) = vp(row,col);
        x(3) = ve(row,col);
        x(4) = 0.0;
        temp = Ct_mTK_bolus_delay_model (t, Cp, x);
        CtMAT(row,col,:) = temp;
    end
end

R1 = ones(dim, dim, num_AIF_pts).*(1./t1) + r1.*CtMAT;
SI = (1-exp(-TR.*R1))./(1-exp(-TR.*R1)*cos(FA));
pre = (1-exp(-TR./t1))./(1-exp(-TR./t1)*cos(FA));
pre(find(t1==0)) = 0;
SI = abs(SI);
% imagesc(pre)
% figure()
% imagesc(SI(:,:,1))
SI(isnan(SI)) = 0;
end