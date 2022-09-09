function YY = dup_folding(Y,blocksize)
% ------------------------------------------------------
% Use this function to duplicate the boundaries of the folding blocks, and
% then swap it with that of the adjacent block. For more details, please
% see [1].
% Current Version only works for 3-rd order tensor.
% 
% [1] Xu, L., Cheng, L., Wong, N., & Wu, Y. C. (2020). Learning tensor train representation with automatic rank determination from incomplete noisy data. arXiv preprint arXiv:2010.06564.
% [2] Xu, L., Cheng, L., Wong, N., & Wu, Y. C. (2021, December). Overfitting Avoidance in Tensor Train Factorization and Completion: Prior Analysis and Inference. In 2021 IEEE International Conference on Data Mining (ICDM) (pp. 1439-1444). IEEE.
%
% ------------------Input------------------
% Y
%       The tensor to be folded
% blocksize
%       The basic block in the original tensor after folding
%       e.g., Y \in 256*256*3, if separations on the rows and columns are as
%       M = [2,2,2,2,2,2,2,2] and N = [2,2,2,2,2,2,2,2], respectively, then
%       the basic block to fold in the original tensor is [M(1),N(1),3].
%
% ------------------Output------------------
% YY
%       With the same order as Y, but with larger size, due to duplicated
%       boundaries along the blocks to fold
% 
% ------------------------------------------------------
% XU Le, 2020, Last update: 9th Sep 2022
% ------------------------------------------------------

Y_size = size(Y);
%% duplicate 3-rd dim
d1 = blocksize(1); dd1 = d1+2;
d2 = blocksize(2); dd2 = d2+2;
d3 = blocksize(3); dd3 = d3+2;

num = [Y_size(1)/d1; Y_size(2)/d2; Y_size(3)/d3];
if num(3) == 1
    if num(2) == 1
        YY_size = [Y_size(1)+num(1)*2; Y_size(2); Y_size(3)];
    else
        YY_size = [Y_size(1)+num(1)*2; Y_size(2)+num(2)*2; Y_size(3)];
    end
else
    YY_size = [Y_size(1)+num(1)*2; Y_size(2)+num(2)*2; Y_size(3)+num(3)*2];
end

YY = zeros(YY_size');
% duplicate the rows
YY_rowdup = zeros([YY_size(1),Y_size(2),Y_size(3)]);
if num(1) ~= 1
    for i = 1:num(1)
        % first line for a block
        if i == 1
            YY_rowdup((i-1)*dd1+1,:,:) = Y((i-1)*d1+1,:,:);
        else
            YY_rowdup((i-1)*dd1+1,:,:) = Y((i-1)*d1,:,:);
        end
        % block copy
        YY_rowdup((i-1)*dd1+(2:d1+1),:,:) = Y((i-1)*d1+(1:d1),:,:);
        % last line for a block
        if i == num(1)
            YY_rowdup((i-1)*dd1+d1+2,:,:) = Y((i-1)*d1+d1,:,:);
        else
            YY_rowdup((i-1)*dd1+d1+2,:,:) = Y(i*d1+1,:,:);
        end
    end
end

% duplicate the columns
if num(2) ~= 1
    YY_columndup = zeros([YY_size(1),YY_size(2),Y_size(3)]);
    for i = 1:num(2)
        % first line for a block
        if i == 1
            YY_columndup(:,(i-1)*dd2+1,:) = YY_rowdup(:,(i-1)*d2+1,:);
        else
            YY_columndup(:,(i-1)*dd2+1,:) = YY_rowdup(:,(i-1)*d2,:);
        end
        % block copy
        YY_columndup(:,(i-1)*dd2+(2:d2+1),:) = YY_rowdup(:,(i-1)*d2+(1:d2),:);
        % last line for a block
        if i == num(2)
            YY_columndup(:,(i-1)*dd2+d2+2,:) = YY_rowdup(:,(i-1)*d2+d1,:);
        else
            YY_columndup(:,(i-1)*dd2+d2+2,:) = YY_rowdup(:,i*d2+1,:);
        end
    end
else
    YY_columndup = YY_rowdup;
end

% % duplicate the rows
% YY_rowdup = zeros([YY_size(1),YY_size(2),Y_size(3)]);
% for i = 1:num(1)
%     % first line for a block
%     if i == 1
%         YY_rowdup((i-1)*dd1+1,:,:) = YY_columndup((i-1)*dd1+2,:,:);
%     else
%         YY_rowdup((i-1)*dd1+1,:,:) = YY_columndup((i-1)*dd1,:,:);
%     end
%     % block copy
%     YY_rowdup((i-1)*dd1+(2:d1+1),:,:) = YY_columndup((i-1)*d1+(1:d1),:,:);
%     % last line for a block
%     if i == num(1)
%         YY_rowdup((i-1)*dd1+d1+2,:,:) = YY_columndup((i-1)*d1+d1,:,:);
%     else
%         YY_rowdup((i-1)*dd1+d1+2,:,:) = YY_columndup(i*d1+1,:,:);
%     end
% end

% duplicate the 3rd dim
if num(3) ~= 1
    for i = 1:num(3)
        % first line for a block
        if i == 1
            YY(:,:,(i-1)*dd3+1) = YY_columndup(:,:,(i-1)*d3+1);
        else
            YY(:,:,(i-1)*dd3+1) = YY_columndup(:,:,(i-1)*d3);
        end
        % block copy
        YY(:,:,(i-1)*dd3+(2:d3+1)) = YY_columndup(:,:,(i-1)*d3+(1:d3));
        % last line for a block
        if i == num(3)
            YY(:,:,(i-1)*dd3+d3+2) = YY_columndup(:,:,(i-1)*d3+d3);
        else
            YY(:,:,(i-1)*dd3+d3+2) = YY_columndup(:,:,i*d3+1);
        end
    end
else
    YY = YY_columndup;
end

end