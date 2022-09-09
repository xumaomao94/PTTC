% ------------------------------------------------------
% Run this demo to test the accuracy of the rank estimation.
% It randomly generate noise 100 times and tests VITTC for different settings
% Under the rank pruning strategy adopted in this demo, you are supposed to
% obtain similar results as in the reference paper.
% ------------------------------------------------------
% Reference
% Xu, L., Cheng, L., Wong, N., & Wu, Y. C. (2020). Learning tensor train representation with automatic rank determination from incomplete noisy data. arXiv preprint arXiv:2010.06564.
% Section 5.1 Validation on synthetic data
% ------------------------------------------------------
% XU Le, 2020, Last Update: 8th Sep 2022
% ------------------------------------------------------
clear all;

addpath(genpath('f_VITTC'));
addpath(genpath('f_datageneration'));


dimension = [20,20,20];
savefolder = 'experiment results\';
prefix = 'exp_rank';

%% Experiment 1: TT rank = 5, SNR = 20dB, MR = 0:20:80 / 100
for tt_r = 5
    tt_rank = [1,tt_r,tt_r,1];
    [A_raw,TT_core] = tt_generate(dimension,tt_rank);
for SNR = 20
for missing_rate = 0:20:80
    
savename = sprintf('%s%s%d%s%d%s%d',savefolder,prefix,tt_r,'N',SNR,'mr',missing_rate)

rse_average = 0; % do average on the rse
rank_average = 0; % the average estimated rank
rank_accuracy = 0; % the percentage of true rank
rank_record = zeros(length(dimension)+1,100); % store the rank

for iter = 1:10
    fprintf('tt_rank: %d, snr: %d, missing rate: %d, iteration: %d\n',tt_r,SNR,missing_rate,iter);
    noise = noise_generate(A_raw,SNR);
    S = binornd(1,1-missing_rate/100,dimension);%%%%%%%%mask
    A_observed = A_raw + noise;
    A_observed = A_observed.*S;

    max_iter = 100; % the max number of iterations of the VI procedure, set to 0 to get the TTSVD result
    [A_completed,Gcore,Lambda,Tau,rse,rank_estimated] = VITTC(A_raw,A_observed,S); % A_raw is adopted to evaluate the result, if you do not have A_raw, just put A_observed in
    
    rse_average = rse_average + rse;
    rank_average = rank_average + rank_estimated;
    rank_record(:,iter) = rank_estimated;
    rank_accuracy = rank_accuracy + double(rank_estimated(:,end)' == tt_rank); % rank_estimated(:,end): is the ultimate estimation of the TTrank
end

rse_average = rse_average/iter;
rank_average = rank_average/iter;
rank_accuracy = rank_accuracy/iter;

save(savename,'rse_average','rank_average','rank_accuracy','rank_record');
end
end
end

%% Experiment 2: TT rank = 5, SNR = 0:5:20, MR = 0/ 100
for tt_r = 5
    tt_rank = [1,tt_r,tt_r,1];
    [A_raw,TT_core] = tt_generate(dimension,tt_rank);
for SNR = 0:5:20
for missing_rate = 0
    
savename = sprintf('%s%s%d%s%d%s%d',savefolder,prefix,tt_r,'N',SNR,'mr',missing_rate)

rse_average = 0; % do average on the rse
rank_average = 0; % the average estimated rank
rank_accuracy = 0; % the percentage of true rank
rank_record = zeros(length(dimension)+1,100); % store the rank

for iter = 1:10
    fprintf('tt_rank: %d, snr: %d, missing rate: %d, iteration: %d\n',tt_r,SNR,missing_rate,iter);
    noise = noise_generate(A_raw,SNR);
    S = binornd(1,1-missing_rate/100,dimension);%%%%%%%%mask
    A_observed = A_raw + noise;
    A_observed = A_observed.*S;
    
    max_iter = 100; % the max number of iterations of the VI procedure, set to 0 to get the TTSVD result
    [A_completed,Gcore,Lambda,Tau,rse,rank_estimated] = VITTC(A_raw,A_observed,S); % A_raw is adopted to evaluate the result, if you do not have A_raw, just put A_observed in
    
    rse_average = rse_average + rse;
    rank_average = rank_average + rank_estimated;
    rank_record(:,iter) = rank_estimated;
    rank_accuracy = rank_accuracy + double(rank_estimated(:,end)' == tt_rank); % rank_estimated(:,end): is the ultimate estimation of the TTrank
end

rse_average = rse_average/iter;
rank_average = rank_average/iter;
rank_accuracy = rank_accuracy/iter;

save(savename,'rse_average','rank_average','rank_accuracy','rank_record');
end
end
end

%% Experiment 3: TT rank = 10:5:20, SNR = 20dB, MR = 20 /100
for tt_r = 15:5:20
    tt_rank = [1,tt_r,tt_r,1];
    [A_raw,TT_core] = tt_generate(dimension,tt_rank);
for SNR = 20
for missing_rate = 20
savename = sprintf('%s%s%d%s%d%s%d',savefolder,prefix,tt_r,'N',SNR,'mr',missing_rate)

rse_average = 0; % do average on the rse
rank_average = 0; % the average estimated rank
rank_accuracy = 0; % the percentage of true rank
rank_record = zeros(length(dimension)+1,100); % store the rank

for iter = 1:10
    fprintf('tt_rank: %d, snr: %d, missing rate: %d, iteration: %d\n',tt_r,SNR,missing_rate,iter);
    noise = noise_generate(A_raw,SNR);
    S = binornd(1,1-missing_rate/100,dimension);%%%%%%%%mask
    A_observed = A_raw + noise;
    A_observed = A_observed.*S;

    max_iter = 100; % the max number of iterations of the VI procedure, set to 0 to get the TTSVD result
    [A_completed,Gcore,Lambda,Tau,rse,rank_estimated] = VITTC(A_raw,A_observed,S); % A_raw is adopted to evaluate the result, if you do not have A_raw, just put A_observed in
    
    rse_average = rse_average + rse;
    rank_average = rank_average + rank_estimated;
    rank_record(:,iter) = rank_estimated;
    rank_accuracy = rank_accuracy + double(rank_estimated(:,end)' == tt_rank); % rank_estimated(:,end): is the ultimate estimation of the TTrank
end

rse_average = rse_average/iter;
rank_average = rank_average/iter;
rank_accuracy = rank_accuracy/iter;

save(savename,'rse_average','rank_average','rank_accuracy','rank_record');
end
end
end



