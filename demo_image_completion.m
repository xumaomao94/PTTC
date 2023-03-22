% ------------------------------------------------------
% Run this demo to test VITTC on image completion
% ------------------------------------------------------
% Reference
% [1] Xu, L., Cheng, L., Wong, N., & Wu, Y. C. (2020). Learning tensor train representation with automatic rank determination from incomplete noisy data. arXiv preprint arXiv:2010.06564.
% [2] Xu, L., Cheng, L., Wong, N., & Wu, Y. C. (2021, December). Overfitting Avoidance in Tensor Train Factorization and Completion: Prior Analysis and Inference. In 2021 IEEE International Conference on Data Mining (ICDM) (pp. 1439-1444). IEEE.
% ------------------------------------------------------
% XU Le, 2020, Last Update: 8th Sep 2022
% ------------------------------------------------------

clear all

addpath(genpath('f_VITTC'));
addpath(genpath('f_datageneration'));
addpath(genpath('f_tensorfolding'));
addpath(genpath('f_perfevaluate'));

img_name = 'TestImages\missing_rate80_1.mat';
IfNoiseOn = true; % set true to test the noisy case
if IfNoiseOn
    % Gaussian noise with mean 0 and variance 0.1
    save_name = 'experiment results\VITTC_N10MR80jellybeans.mat';
else
    save_name = 'experiment results\VITTC_clean80jellybeans.mat';
end

%% Read the image: 'jellybeans' image, with size 256*256*3
load(img_name);
X = img.*255;  % the original image
O = S;  %  the boolean indicator tensor, showing which element is observed
if ~IfNoiseOn
    Y = X.*S;  % the observed image
else
    X_noise = imnoise(img,'gaussian',0,0.1).*255;
    Y = X_noise.*S;
end


%% Preprocessing ---- Improved ket-folding
M_separate = 2*ones(1,8); % the rows: 256 -> [2,2,2,2,2,2,2,2]
N_separate = 2*ones(1,8); % the columns: 256 -> [2,2,2,2,2,2,2,2]

% duplicate the boundaries to fold
basic_folding_block = [M_separate(1),N_separate(1),3]; % [2,2,3] in this demo
Y_dup = dup_folding(Y,basic_folding_block);
O_dup = dup_folding(O,basic_folding_block);
M_separate(1) = M_separate(1)+2; % since the basic blocks are with boundaries duplicated
N_separate(1) = N_separate(1)+2;

% Ket-folding on the tensor with duplicated folding boundaries
Y_ket = mat2tensor(Y_dup,M_separate,N_separate);
O_ket = mat2tensor(O_dup,M_separate,N_separate);

%% run VITTC
tic
[X_VITT_ket,Gcore,Lambda,Tau,~,rank_estimated] = VITTC(Y_ket,Y_ket,O_ket,...
                                                        'initmethod','randinit',...
                                                        'maxrank',64,...
                                                        'num_of_iter',30,...
                                                        'thre_stop',1e-6,...
                                                        'show_info', true);
toc
%% from Ket-folding & Dup state back to an image
X_VITT_dup = tensor2mat(X_VITT_ket,M_separate,N_separate);
X_VITT = dup_clear(X_VITT_dup,basic_folding_block);

%% evaluation
rse = rse_score(X_VITT,X);
psnr = psnr_score(X_VITT,X);
ssim = ssim_index(rgb2gray(uint8(X_VITT)),rgb2gray(uint8(X)));

%% save the result
figure;
subplot(1,3,1); imshow(X./255);
subplot(1,3,2); imshow(Y./255);
subplot(1,3,3); imshow(X_VITT./255);

% save(save_name,'X_VITT','rse','psnr','ssim')

