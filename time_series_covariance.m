function [cov] = time_series_covariance(time_series)
dim = size(time_series, 1);
time_pts = size(time_series, 3);

mean_img = mean(time_series, 3);

norm_series = time_series - mean_img;
%norm_series = time_series;

% for n = 1:time_pts
%     imagesc(norm_series(:,:,n))
% %     pause(0.03)
% end

vec_norm = reshape(norm_series, dim*dim, []);

cov = zeros(dim*dim, 1);
num = 0;
for n = 1:time_pts
%     if(n/3 == floor(n/3))
%         continue
%     end
    cov = cov + vec_norm(:,n) * vec_norm(:,n)';
    num = num + 1;
end
num
cov = cov ./ num;

ident = eye(dim*dim);

vec_img = (ident .* cov) * ones(dim*dim, 1);
size(vec_img)

img = reshape(vec_img, dim, dim);
%figure()
%imagesc(img)
end