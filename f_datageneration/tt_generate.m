function [A,TT_core] = tt_generate(dim,tt_rank)
% ------------------------------------------------------
% Randomly generate a tensor in the TT format, given its dimension and TT rank
% 
% ------------------Input------------------
% dim: the dimension of the tensor, e.g. [10,10,10]
% tt_rank: the tt rank, e.g. [1,5,5,1]
%         1. the length of the tt rank vector should be one more that that of
%         dim
%         2. the first and last tt rank must be 1
% 
% ------------------Output------------------
% A: the generated tt-format tensor
% TT_core: tt cores used to generate the tensor, stored as cells
% 
% XU Le, 2020
% ------------------------------------------------------
    if length(dim) ~= length(tt_rank) - 1
        msg = 'error: the length of dimension should be one less than the tt_rank';
        error(msg)
    end
    if tt_rank(1) ~= 1
        msg = 'error: the first tt_rank must be 1';
        error(msg)
    end
    if tt_rank(length(tt_rank)) ~= 1
        msg = 'error: the last tt_rank must be 1';
        error(msg)
    end
    
    % TT_cores generation
    TT_core = cell(length(dim),1);
    for i = 1:length(dim)
        TT_core{i} = randn(tt_rank(i),tt_rank(i+1),dim(i));
    end
    
    % tensor generation according to the TT_cores
    A = 1;
    for i = 1:length(dim)
        core2mat = reshape(TT_core{i},[tt_rank(i),tt_rank(i+1)*dim(i)]);
        A = A * core2mat; % I1..In-1 * r(n)In
        A = reshape(A,[prod(dim(1:i-1)),tt_rank(i+1),dim(i)]);
        A = permute(A,[1,3,2]);
        A = reshape(A,[prod(dim(1:i)),tt_rank(i+1)]);
    end
    A = reshape(A,dim);
end