function [A_completed,Gcore,Lambda,Tau,rse,rank_estimated,Power] = VITTC(A_raw,A_observed,Mask,varargin)
% ------------------------------------------------------
% Variational Inference for the Tensor Train completion based on the model
% proposed in the following paper
% 
% [1] Xu, L., Cheng, L., Wong, N., & Wu, Y. C. (2020). Learning tensor train representation with automatic rank determination from incomplete noisy data. arXiv preprint arXiv:2010.06564.
% [2] Xu, L., Cheng, L., Wong, N., & Wu, Y. C. (2021, December). Overfitting Avoidance in Tensor Train Factorization and Completion: Prior Analysis and Inference. In 2021 IEEE International Conference on Data Mining (ICDM) (pp. 1439-1444). IEEE.
% 
%------------------Update in 9th Sep, 2022------------------
% Some parts of the code, like initialization and rank pruning strategy are improved to
% make the algorithm stable under different testing cases. With this updating, similar
% results shall be obtained with those in the afore-listed papers. For experiments on synthetic
% data, they are supposed show better performance, e.g., under SNR=20&MR=80, according to our test.
% For experiments on image completion, we suggest to set 'method_rank_prune' as
% 'powerBased' to enable automatic rank determination and faster inference,
% which is supposed to perform a little better under noisy cases, while a
% little worse with clean data. For the same setting
% as that in [1]&[2], please set 'method_rank_prune' as 'none'.
%
% ------------------Input------------------
% A_raw
%       The original tensor data, only used for performance test, use
%       A_observed instead if there is no A_raw
% A_observed
%       The noisy, incomplete tensor data
% Mask
%       Indicating tensor,
%       -1 -> entry observed
%       -0 -> not observed
% initialmethod
%       -'svdinit'-> use Gaussian random variables to fill in empty entries
%       -'randinit'-> all entries are initialized by Gaussian variables
% num_of_iter
%       Max iteration for the VI update
% method_rank_prune
%       -'none' -> do not prune ranks during the inference
%       -'lambdaBased' -> prune ranks according to lambda
%       -'powerBased' -> (*recommended) prune ranks according to the TT core slice power, if average slice power/average core power < threshold, then discard a slice
% thre_rank_prune
%       -method_rank_prune='lambdaBased'
%           Slices with average power thres-times larger than the minimum
%           average slice power will be discarded,  slice power is obtained
%           from lambda
%           Recommended value: larger than 50, e.g., 100
%       -method_rank_prune='powerBased'
%           Slices with (average power/average core power) less than thres will be discarded
%           Recommended value: smaller than 1e-2, e.g., 1e-3, 1e-5
% thre_stop
%       stop the iteration when relative square error between current recovered tensor and
%       last update is smaller than thre_stop
% optional: maxrank
%       max rank for the TT ranks
% ------------------Output------------------
% A_completed
%       Estimated tensor data
% Gcore
%       Estimated TT cores
% Lambda
%       The mean of the Gamma varaibles
% Tau
%       The estimated noise power
% rse
%       rse between A_raw & A_completed
% rank_estimated
%       The estimated rank
% Power
%       The average power for each TT rank
% 
% ------------------------------------------------------
% XU Le, 2020, Last update: 9th Sep 2022
% ------------------------------------------------------

%% load/set the parameters for the algorithm
p = inputParser;
defaultPar.InitMethod = 'svdinit';
defaultPar.MaxIter = 100;
defaultPar.RankPruneMethod = 'powerBased';
defaultPar.RankPruneThre = 1e-2;
defaultPar.IterEndThre = 1e-12;
defaultPar.ShowInfo = false;
defaultPar.MaxRank = 64;

addRequired(p,'A_raw',@isnumeric);
addRequired(p,'A_observed',@isnumeric);
addRequired(p,'Mask',@(x) isnumeric(x)||islogical(x))
addOptional(p,'initmethod',defaultPar.InitMethod,@(x) ismember(x,{'svdinit','randinit'}));
addOptional(p,'num_of_iter',defaultPar.MaxIter,@isscalar);
addOptional(p,'method_rank_prune',defaultPar.RankPruneMethod, @(x) ismember(x,{'none','powerBased','lambdaBased'}));
addOptional(p,'thre_rank_prune',defaultPar.RankPruneThre, @isscalar);
addOptional(p,'thre_stop',defaultPar.IterEndThre,@isscalar);
addOptional(p,'show_info',defaultPar.ShowInfo,@islogical);
addOptional(p,'maxrank',defaultPar.MaxRank,@isscalar);
parse(p,A_raw,A_observed,Mask,varargin{:})

Par = p.Results;
if Par.show_info
    disp(Par)
end


%% Initialization
Size_A = size(A_observed);
indnorm = 10^(ndims(A_observed)-1)/max( abs(A_observed(:)) );
A_observed = A_observed.*indnorm;
A_Ctemp = 0;
[Lambda,Tau,Gcore] = VITTC_initialize_indpd(A_observed,Mask,Par.initmethod,Par.maxrank);

%% VI iteration
for i = 1:Par.num_of_iter
    [Gcore,V_final,W_final] = update_gcore_indpd_C_var(A_observed,Mask,Gcore,Lambda,Tau);
    Lambda = update_lambda_indpd_C_var(A_observed,Gcore,Lambda,Tau);
    Tau = update_tau_indpd_C_var(A_observed,Mask,Gcore,Lambda,Tau,V_final,W_final);
    if strcmp(Par.method_rank_prune,'lambdaBased')
        if Par.thre_rank_prune<10 && i==1
            fprintf('not the reasonable threshold for rank reduction, for lambdaBased rank pruning, threshold is better set as numbers like 100')
        end
        [Gcore,Lambda] = rank_reduce_relative_indpd_C_var(Gcore,Lambda,Par.thre_rank_prune);
    elseif strcmp(Par.method_rank_prune,'powerBased')
        if Par.thre_rank_prune>0.1 && i==1
            fprintf('not the reasonable threshold for rank reduction, for powerBased rank pruning, threshold is better set as numbers like 1e-5')
        end
        [Gcore,Lambda] =  rank_reduce_relative_ind_Power(Gcore,Lambda,Par.thre_rank_prune);
    end
    
    A_completed = tt2full(Gcore,Size_A)./indnorm;
    dist = sumsqr(A_completed(:)-A_Ctemp(:))/sumsqr(A_completed(:));
    if Par.show_info
        fprintf('the %d th iteration, rse between current and last update: %.9f \n',i,dist)
    end
    if dist < Par.thre_stop
        break
    end
    A_Ctemp = A_completed;
end

%% recover the tensor, and get the rse from the raw data
A_completed = tt2full(Gcore,Size_A);
A_completed = A_completed./indnorm;
rse = sumsqr(A_completed-A_raw)/sumsqr(A_raw);

% mean power of corresponding slices w.r.t. TT ranks
ndims_A  = ndims(A_observed);
Power_L = cell(1,ndims_A);
Power_H = cell(1,ndims_A);
for order = 1:ndims_A
    Power_L{order} = zeros(size(Lambda.mean{order}));
    Power_H{order} = zeros(size(Lambda.mean{order+1}));
    meansqr = Gcore.mean{order}(:)'*Gcore.mean{order}(:)/numel(Gcore.mean{order});
    for r = 1:length(Power_L{order})
        gcoreslice = Gcore.mean{order}(r,:,:);
        Power_L{order}(r) = gcoreslice(:)'*gcoreslice(:)/meansqr/numel(gcoreslice);
    end
    for r = 1:length(Power_H{order})
        gcoreslice = Gcore.mean{order}(:,r,:);
        Power_H{order}(r) = gcoreslice(:)'*gcoreslice(:)/meansqr/numel(gcoreslice);
    end
end
Power = cell(1,ndims_A); Power{1} = 1;
for order = 2:ndims_A
    Power{order} = Power_L{order}+Power_H{order-1};
    if Par.show_info
%         subplot(1,ndims_A,order-1); bar(Power{order});
    end
end

%% the final guessed rank
rank_estimated = ones(ndims(A_observed)+1,1);
for order = 1:ndims(A_observed)
    rank_estimated(order) = size(Gcore.mean{order},1);
end

end