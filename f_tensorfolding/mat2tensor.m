function A_tensor = mat2tensor(A_matrix,M,N)
% ------------------------------------------------------
% For a low-order tensor A, with size [prod(M),prod(N),...], turn its first two orders into a length(M)-order
% length(M) should be equal to length(N)
% the transformed tensor is with size [M.*N,...]
%
% ------------------------------------------------------
% XU Le, 2020, Last update: 9th Sep 2022
% ------------------------------------------------------

%% The new one
    if length(M) ~= length(N)
        msg = 'Length of the disparting images should be the same';
        error(msg);
    end
    
    if prod(M) ~= size(A_matrix,1)
        msg = 'prod(M) should be equal to the first dimension of A_matrix';
        error(msg);
    end
    
    if prod(N) ~= size(A_matrix,2)
        msg = 'prod(N) should be equal to the second dimension of A_matrix';
        error(msg);
    end
    
    partNum = length(M);
    size_A = size(A_matrix);
    
    % Split the matrix into small blocks
    size_A_tensor = [M,N,size_A(3:end)];
    A_tensor = reshape(A_matrix,size_A_tensor);
    
    % Change the permutation so that elements in a small block are in the
    % same order
    permu = 1:2*partNum+(length(A_matrix)-2);
    permu(1:2:2*partNum) = 1:partNum;
    permu(2:2:2*partNum) = partNum+1:2*partNum;
    A_tensor = permute(A_tensor,permu);
    
    size_A_tensor = [M.*N,size_A(3:end)];
    A_tensor = reshape(A_tensor,size_A_tensor);
    

end