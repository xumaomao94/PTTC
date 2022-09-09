function [Gcore_update,V_final,W_final] = update_gcore_indpd_C_var(A_observed,Mask,Gcore,Lambda,Tau)

    size_A = size(A_observed);
    ndims_A = ndims(A_observed);
    R = zeros(1,ndims_A+1);
    for i = 1:ndims_A+1
        R(i) = length(Lambda.mean{i});
    end
   
    % the update for the order-th TTcore
    V = 1; % to store the kroncor's product before the current core
    W = 1; % to store the mean's product before the current core
    
    U = cell(1,ndims_A); % to store the kroncor's product after the current core
    U{ndims_A} = 1;
    X = cell(1,ndims_A); % to store the mean's product after the current core
    X{ndims_A} = 1;
    
    for order = 1:ndims_A-1
        i_locate = ndims_A-order; % from ndims_A-1 to 1
        i_after = i_locate+1; % index of the core after i_locate
        
        % update U
        Gcorevar_kronform = spalloc(R(i_after)^2,R(i_after+1)^2,R(i_after)*R(i_after+1));%zeros(R(i_after)^2,R(i_after+1)^2);
        for d = 1:size_A(i_after)
            Gcorevar_kronform(1:R(i_after)+1:R(i_after)^2,1:R(i_after+1)+1:R(i_after+1)^2) = Gcore.var{i_after}(:,:,d);
            Gcorecor_kronform = kron(Gcore.mean{i_after}(:,:,d),Gcore.mean{i_after}(:,:,d)) + Gcorevar_kronform;
            U{i_locate}((d-1)*R(i_after)^2+1:d*R(i_after)^2,:) = Gcorecor_kronform * U{i_after};
        end
        U{i_locate} = reshape(U{i_locate},[R(i_after)^2,prod(size_A(i_after:end))]);

        % update X
        X{i_locate} = reshape(permute(Gcore.mean{i_after},[1,3,2]),[R(i_after)*size_A(i_after),R(i_after+1)]) * X{i_after};
        X{i_locate} = reshape(X{i_locate},[R(i_after),prod(size_A(i_after:end))]);
    end
    
    for order = 1:ndims_A
        % update the precision matrix
        
        if order ~= 1
            Gcorevar_kronform = spalloc(R(order-1)^2,R(order)^2,R(order-1)*R(order)); %zeros(R(order-1)^2,R(order)^2);
            
            V_temp = V;
            V = zeros(prod(size_A(1:order-1))/size_A(order-1),R(order)^2*size_A(order-1));
            for d = 1:size_A(order-1)
                Gcorevar_kronform(1:R(order-1)+1:R(order-1)^2,1:R(order)+1:R(order)^2) = Gcore.var{order-1}(:,:,d);
                Gcorecor_kronform = kron(Gcore.mean{order-1}(:,:,d),Gcore.mean{order-1}(:,:,d)) + Gcorevar_kronform;
                V(:,(d-1)*R(order)^2+1:d*R(order)^2) = V_temp*Gcorecor_kronform;
            end
            V = reshape( permute( reshape(V,[prod(size_A(1:order-1))/size_A(order-1),R(order)^2,size_A(order-1)] ),[1,3,2] ), [prod(size_A(1:order-1)),R(order)^2] );
        end
        
        M_2matrix = reshape(Mask,[prod(size_A(1:order-1)),size_A(order),prod(size_A(order+1:end))]);
        M_2matrix = permute(M_2matrix,[1,3,2]);
     
        if order ~= 1
            W = W * reshape(Gcore.mean{order-1},R(order-1),R(order)*size_A(order-1));
            W = reshape(permute(reshape(W,[prod(size_A(1:order-2)),R(order),size_A(order-1)]),[1,3,2]),[prod(size_A(1:order-1)),R(order)]);
        end
        
        A_2matrix = reshape(A_observed,[prod(size_A(1:order-1)),size_A(order),prod(size_A(order+1:end))]);
        A_2matrix = permute(A_2matrix,[1,3,2]);
        
        Gcoreprecision_Lambdapart = Lambda.mean{order}'*Lambda.mean{order+1};
        
        Tau_times_M = Tau.mean * M_2matrix;
        for d = 1:size_A(order)
            Gcoreprecision_VUpart = V' * Tau_times_M(:,:,d) * U{order}';
            Gcore_mean_times_precision =  W' * (A_2matrix(:,:,d).*Tau_times_M(:,:,d)) * X{order}';
            for p = 1:R(order)
                for q = 1:R(order+1)
                    Gcore.var{order}(p,q,d) = 1/(Gcoreprecision_Lambdapart(p,q)+Gcoreprecision_VUpart((p-1)*R(order)+p,(q-1)*R(order+1)+q));
                    Gcore.mean{order}(p,q,d) = (  - ( sum(sum(Gcoreprecision_VUpart((p-1)*R(order)+1:p*R(order),(q-1)*R(order+1)+1:q*R(order+1)).*Gcore.mean{order}(:,:,d))) - Gcoreprecision_VUpart((p-1)*R(order)+p,(q-1)*R(order+1)+q)*Gcore.mean{order}(p,q,d)) ...
                        + Gcore_mean_times_precision(p,q)  ) * Gcore.var{order}(p,q,d);
                end
            end
        end
    end
    
    % saved for tau's use
    Gcorevar_kronform = spalloc(R(ndims_A)^2,R(ndims_A+1)^2,R(ndims_A)*R(ndims_A+1));%zeros(R(ndims_A)^2,R(ndims_A+1)^2);
    
    V_temp = V;
    V = zeros(prod(size_A(1:ndims_A))/size_A(ndims_A),R(ndims_A+1)^2*size_A(ndims_A));
    for d = 1:size_A(ndims_A)
        Gcorevar_kronform(1:R(ndims_A)+1:R(ndims_A)^2,1:R(ndims_A+1)+1:R(ndims_A+1)^2) = Gcore.var{ndims_A}(:,:,d);
        Gcorecor_kronform = kron(Gcore.mean{ndims_A}(:,:,d),Gcore.mean{ndims_A}(:,:,d)) + Gcorevar_kronform;
        V(:,(d-1)*R(ndims_A+1)^2+1:d*R(ndims_A+1)^2) = V_temp*Gcorecor_kronform;
    end
    V = reshape( permute( reshape(V,[prod(size_A(1:ndims_A))/size_A(ndims_A),R(ndims_A+1)^2,size_A(ndims_A)] ),[1,3,2] ), [prod(size_A(1:ndims_A)),R(ndims_A+1)^2] );
    V_final = V;
    
    W = W * reshape(Gcore.mean{ndims_A},R(ndims_A),R(ndims_A+1)*size_A(ndims_A));
    W = reshape(permute(reshape(W,[prod(size_A(1:ndims_A-1)),R(ndims_A+1),size_A(ndims_A)]),[1,3,2]),[prod(size_A(1:ndims_A)),R(ndims_A+1)]);
    W_final = W;
    
    Gcore_update = Gcore;
end