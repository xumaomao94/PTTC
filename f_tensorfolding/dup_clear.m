function Y = dup_clear(YY,blocksize)
% ------------------------------------------------------
% Use this function to clear the boundaries of the folding blocks. For more details, please
% see [1].
% Please only use this function only after you adopted the function 'dup_folding'
% Current Version only works for 3-rd order tensor.
% 
% [1] Xu, L., Cheng, L., Wong, N., & Wu, Y. C. (2020). Learning tensor train representation with automatic rank determination from incomplete noisy data. arXiv preprint arXiv:2010.06564.
% [2] Xu, L., Cheng, L., Wong, N., & Wu, Y. C. (2021, December). Overfitting Avoidance in Tensor Train Factorization and Completion: Prior Analysis and Inference. In 2021 IEEE International Conference on Data Mining (ICDM) (pp. 1439-1444). IEEE.
%
% ------------------Input------------------
% YY
%       The tensor with duplicated boundaries added by dup_folding(Y,blocksize)
% blocksize
%       The basic block in the original tensor after folding
%       e.g., Y \in 256*256*3, if separations on the rows and columns are as
%       M = [2,2,2,2,2,2,2,2] and N = [2,2,2,2,2,2,2,2], respectively, then
%       the basic block to fold in the original tensor is [M(1),N(1),3].
%
%       Remember to use the same block as in dup_folding, instead of
%       [M(1)+2,N(1)+2,3].
%
% ------------------Output------------------
% Y
%       The original tensor
% 
% ------------------------------------------------------
% XU Le, 2020, Last update: 9th Sep 2022
% ------------------------------------------------------
YY_size = size(YY);

d1 = blocksize(1); dd1 = d1+2;
d2 = blocksize(2); dd2 = d2+2;
d3 = blocksize(3); dd3 = d3+2;

num = [YY_size(1)/dd1; YY_size(2)/dd2; YY_size(3)/dd3];
if num(3) <= 1
    if num(2) <= 1
        Y_size = [YY_size(1)-num(1)*2; YY_size(2); YY_size(3)];
    else
        Y_size = [YY_size(1)-num(1)*2; YY_size(2)-num(2)*2; YY_size(3)];
    end
else
    Y_size = [YY_size(1)-num(1)*2; YY_size(2)-num(2)*2; YY_size(3)-num(3)*2];
end

% do average on the duplicated pixels
% YY(2+dd1:dd1:end,:,:) = (YY(2+dd1:dd1:end,:,:) + YY(dd1:dd1:end-1,:,:))/2;
% YY(dd1-1:dd1:end-2,:,:) = (YY(dd1-1:dd1:end-2,:,:) + YY(dd1+1:dd1:end,:,:))/2;
% YY(:,2+dd2:dd2:end,:) = (YY(:,2+dd2:dd2:end,:) + YY(:,dd2:dd2:end-1,:))/2;
% YY(:,dd2-1:dd2:end-2,:) = (YY(:,dd2-1:dd2:end-2,:) + YY(:,dd2+1:dd2:end,:))/2;
% if num(3) > 1
%     YY(:,:,dd3+2:dd3:end) = (YY(:,:,dd3+2:dd3:end) + YY(:,:,dd3:dd3:end-1))/2;
%     YY(:,:,dd3-1:dd3:end-2) = (YY(:,:,dd3-1:dd3:end-2) + YY(:,:,dd3+1:dd3:end))/2;
% end




YY(2+dd1:dd1:end,:,:) = (YY(2+dd1:dd1:end,:,:) );
YY(dd1-1:dd1:end-2,:,:) = (YY(dd1-1:dd1:end-2,:,:) );
YY(:,2+dd2:dd2:end,:) = (YY(:,2+dd2:dd2:end,:) );
YY(:,dd2-1:dd2:end-2,:) = (YY(:,dd2-1:dd2:end-2,:) );
if num(3) > 1
    YY(:,:,dd3+2:dd3:end) = (YY(:,:,dd3+2:dd3:end) );
    YY(:,:,dd3-1:dd3:end-2) = (YY(:,:,dd3-1:dd3:end-2) );
end

v1 = 1:YY_size(1);
v1([1:dd1:end,dd1:dd1:end]) = [];
v2 = 1:YY_size(2);
if num(2) > 1
    v2([1:dd2:end,dd2:dd2:end]) = [];
end
v3 = 1:YY_size(3);
if num(3) > 1
    v3([1:dd3:end,dd3:dd3:end]) = [];
end
Y = YY(v1,v2,v3);

end