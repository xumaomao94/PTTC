function A_mat = tensor2mat(A_tensor,M,N)
% ------------------------------------------------------
% For a high-order tensor, with size [M.*N,...],turn its first length(M)
% orders into 2 orders
% length(M) should be equal to length(N)
% the transformed low-order tensor is with size [prod(M),prod(N),...]
% ------------------------------------------------------
% XU Le, 2020, Last update: 9th Sep 2022
% ------------------------------------------------------

    if length(M) ~= length(N)
        msg = 'Length of the disparting images should be the same';
        error(msg);
    end
    
    partNum = length(M);
    size_A = size(A_tensor);
    if sum( M.*N ~= size_A(1:partNum) ) ~= 0
        msg = 'prod(M.*N) should be equal to the length(M) dimensions of A_tensor';
        error(msg);
    end
    
    size_A_mat = zeros(1,2*partNum);
    size_A_mat(1:2:end) = M;
    size_A_mat(2:2:end) = N;
    size_A_mat = [size_A_mat,size_A(partNum+1:end)];
    A_mat = reshape(A_tensor,size_A_mat);
    
    permu = 1:2*partNum + length(size_A) - partNum;
    permu(1:partNum) = 1:2:2*partNum;
    permu(partNum+1:2*partNum) = 2:2:2*partNum;
    
    A_mat = permute(A_mat,permu);
    A_mat = reshape(A_mat,[prod(M),prod(N),size_A(partNum+1:end)]);


end