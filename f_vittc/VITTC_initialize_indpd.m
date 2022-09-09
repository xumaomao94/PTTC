function [lambda,tau,Gcore] = VITTC_initialize_indpd(A_observed,Mask,initmethod,maxrank)
    Size_A = size(A_observed);
    ndims_A = ndims(A_observed);
    
    % Gcore mean, determine R at the same time
    R = zeros(ndims_A+1,1);
    R(1) = 1; R(end) = 1;
    Gcore.mean = cell(1,ndims_A);
    Gcore.var = cell(1,ndims_A);
    
    
    Mask_index = find(Mask~=0);
    if strcmp(initmethod,'svdinit')
        A_guess = A_observed.*Mask+(1-Mask).*normrnd(mean(A_observed(Mask_index)),std(A_observed(Mask_index)),Size_A);
    elseif strcmp(initmethod,'randinit')
        A_guess = normrnd(mean(A_observed(Mask_index)),std(A_observed(Mask_index)),Size_A);% normal random, with A_mean and A_var
    end
    
    Aremain_M = reshape(A_guess,Size_A(1),[]); % use A_observed + random guess as initialization
    for i = 1:ndims_A-1
        [U,S,V] = svd(Aremain_M,'econ');
        R(i+1) = min(size(S,1),maxrank);
        
        U = U*S.^(1/ndims_A);
        S = S.^((ndims_A-1)/ndims_A);
        Gcore.mean{i} = permute(reshape(U(:,1:R(i+1)),[R(i),Size_A(i),R(i+1)]),[1,3,2]);
        Aremain_M = reshape(S(1:R(i+1),1:R(i+1))*V(:,1:R(i+1))',[R(i+1)*Size_A(i+1),prod(Size_A(i+1:end))/Size_A(i+1)]);
    end
    Gcore.mean{ndims_A} = reshape(Aremain_M,[R(ndims_A),R(ndims_A+1),Size_A(ndims_A)]);
    
    % lambda
    lambda.a = 10^(-6);
    lambda.b = 10^(-6);
    lambda.mean = cell(1,ndims_A+1); % for conviency, we set lambda(0) = 1 and lambda(N) = 1. The real lambda(i) is denoted as lambda(i+1) here!
    lambda.var = cell(1,ndims_A+1);
    for i = 1:(ndims_A+1)
        lambda.mean(i) = {lambda.a/lambda.b * ones(1,R(i))};
        lambda.var(i) = {lambda.a/(lambda.b^2) * ones(1,R(i))};
    end
    
    % tau
    tau.a = 10^(-6);
    tau.b = 10^(-6);
    tau.mean = tau.a/tau.b;
    tau.var = tau.a/(tau.b^2);
    
    
    for i = 1:ndims_A
        for d = 1:Size_A(i)
            Gcore.var{i}(:,:,d) = 1./(lambda.mean{i}'*lambda.mean{i+1});
        end
    end
end