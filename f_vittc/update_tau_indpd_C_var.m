function Tau_update = update_tau_indpd_C_var(A_observed,Mask,Gcore,Lambda,Tau,V_final,W_final)
    % the same as the fix in update_Gcore_fix4, replace enumerations by matrix product, fixed at 17:31, June 12th
    ndims_A = ndims(A_observed);
    size_A = size(A_observed);
    R = zeros(1,ndims_A+1);
    for i = 1:ndims_A+1
        R(i) = length(Lambda.mean{i});
    end
    
    alpha = sum(reshape(Mask,1,[]))/2 + Tau.a;
    beta = Tau.b;
    
%     gmean = tt2full(Gcore,size_A);
%     beta = beta + 1/2*sum(reshape(gmean.^2.*Mask,[],1)) + ...
%         1/2*sum(reshape(A_observed.^2.*Mask,[],1)) - ...
%         sum(reshape(A_observed.*gmean.*Mask,[],1));
    

    gmean = W_final;
    gcor = V_final;
    M_2vec = Mask(:);
    A_2vec = A_observed(:);
    beta = beta + 1/2*M_2vec'*(A_2vec.*A_2vec+gcor-2*A_2vec.*gmean);
    
    
    Tau.mean = alpha/beta;
    Tau.var = alpha/(beta^2);
    Tau_update = Tau;
end