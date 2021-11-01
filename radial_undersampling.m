function [usmat, shifted_usmat] = radial_undersampling(img_dim, num_spokes, start_angle, increment)
%
% This file takes in an image, and returns the undersampled Fourier data
% according to radial undersampling with a given number of spokes and
% at a given start angle
usmat = zeros(img_dim, img_dim);


theta = start_angle;
for i = 1:num_spokes
    % Get the slope. if the slope is > 1, I should invert it and operate on
    % x vs y rather than y vs x in order to hit every pixel
    invert = false;
    m = tan(theta);
    if abs(m) > 1
        invert = true;
        m = 1/m;
    end

    % Caclulate the y intercept
    % Its gotta pass through the middle
    b = img_dim/2 - m*img_dim/2;
    
    % Maybe change ceil later
    for j = 1:(img_dim)
        x = j;
        y = round(m*x + b);

        if invert == true
            if y == 0
                y = 1;
            end
            usmat(x, y) = 1;
        else
            if y == 0
                y = 1;
            end
            usmat(y, x) = 1;
        end
    end

    % Increment theta by the increment
    theta = theta + increment;
    shifted_usmat = circshift(usmat, img_dim/2, 1);
    shifted_usmat = circshift(shifted_usmat, img_dim/2, 2);
end

end