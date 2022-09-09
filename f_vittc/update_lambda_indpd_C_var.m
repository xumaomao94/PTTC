function Lambda_update = update_lambda_ind_var(A_observed,Gcore,Lambda,Tau)
    ndims_A = ndims(A_observed);
    size_A = size(A_observed);
    R = zeros(1,ndims_A+1);
    for i = 1:ndims_A+1
        R(i) = length(Lambda.mean{i});
    end
    
%   when i = 1
%     Gcoreselfcor_sum_after = sum(Gcore.mean{1}.^2, 3) + size_A(1)*Gcore.var{1};
%     alpha = 1/2 * size_A(1) * R(2) + Lambda.a;
%     beta = 1/2 * Gcoreselfcor_sum_after * Lambda.mean{2}' + Lambda.b;
%     Lambda.mean{1} = alpha/beta;
%     Lambda.var{1} = alpha/(beta.^2);
    
    for i = 2:ndims_A % Lambda{1} and Lambda{ndims_A} are set to 1
%         alpha = zeros(1,R(i));
%         beta = zeros(1,R(i));
%         Gcorecor_before = Gcore.cor{i-1};
%         Gcoreselfcor_sum_before = reshape( diag(sum(Gcorecor_before,3)),R(i-1),R(i));
%         Gcorecor_after = Gcore.cor{i};
%         Gcoreselfcor_sum_after = reshape( diag(sum(Gcorecor_after,3)),R(i),R(i+1));
        
        Gcoreselfcor_sum_before = sum( Gcore.mean{i-1}.^2, 3 ) + sum(Gcore.var{i-1},3);
        Gcoreselfcor_sum_after = sum(Gcore.mean{i}.^2, 3) + sum(Gcore.var{i},3);
%         for r = 1:R(i)
            % ************************** need to check whether is still true when i = 1 and i = R(N)
%             alpha(r) = 1/2 * size_A(i-1) * R(i-1) + 1/2 * size_A(i) * R(i+1) + Lambda.a;
%             beta(r) = 1/2 * Gcoreselfcor_sum_before(:,r)' * Lambda.mean{i-1}' + ...
%                 1/2 * Gcoreselfcor_sum_after(r,:) * Lambda.mean{i+1}' + Lambda.b;
            
        alpha = 1/2 * size_A(i-1) * R(i-1) + 1/2 * size_A(i) * R(i+1) + Lambda.a;
        beta = 1/2 * Gcoreselfcor_sum_before' * Lambda.mean{i-1}' + 1/2 * Gcoreselfcor_sum_after * Lambda.mean{i+1}' + Lambda.b;
%         end
%         Lambda.mean{i} = alpha./beta;
%         Lambda.var{i} = alpha./(beta.^2);

        Lambda.mean{i} = (alpha./beta)';
        Lambda.var{i} = (alpha./(beta.^2))';

    end
    
%   when i = N+1
%     Gcoreselfcor_sum_before = sum( Gcore.mean{ndims_A}.^2, 3 ) + size_A(ndims_A)*Gcore.var{ndims_A};
%     alpha = 1/2 * size_A(ndims_A) * R(ndims_A) + Lambda.a;
%     beta = 1/2 * Gcoreselfcor_sum_before' * Lambda.mean{ndims_A}' + Lambda.b;
%     Lambda.mean{ndims_A+1} = alpha/beta;
%     Lambda.var{ndims_A+1} = alpha/(beta.^2);
    
    Lambda_update = Lambda;
end