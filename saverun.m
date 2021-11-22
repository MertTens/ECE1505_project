destdir = 'runs/2_spokes_MLE';

mkdir(destdir);

path = strcat(destdir, '/A.mat');
save(path, 'A');

path = strcat(destdir, '/mse.mat');
save(path, 'mse');

path = strcat(destdir, '/recon.mat');
save(path, 'recon');

path = strcat(destdir, '/obsv.mat');
save(path, 'obsv');

path = strcat(destdir, '/truth.mat');
save(path, 'truth');