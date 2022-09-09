function [Gcore_reduced,Lambda_reduced] = rank_reduce_relative_indpd_C_var(Gcore,Lambda,threshold)
% discard the rows: lambda > threshold*min(lambda)
    ndims_A = length(Gcore.mean);
    for i = 1:ndims_A
        Gcore_mean_temp = Gcore.mean{i};
        Gcore_var_temp = Gcore.var{i};

%         Gcore_cor_temp = Gcore.cor{i};
%         Gcore_kron_temp = Gcore.kron{i};
        lambda_row_temp = Lambda.mean{i};
        lambda_column_temp = Lambda.mean{i+1};
        lambda_var_row_temp = Lambda.var{i};
        
        row_discarded = lambda_row_temp > threshold*min(lambda_row_temp);
        ind_row_discarded = find(row_discarded);
        column_discarded = lambda_column_temp > threshold*min(lambda_column_temp);
        ind_column_discarded = find(column_discarded);
        
        ind_rowkron_discarded = find(~kron(~row_discarded,~row_discarded));
        ind_colkron_discarded = find(~kron(~column_discarded,~column_discarded));
%         r_row = length(lambda_row_temp);
%         r_column = length(lambda_column_temp);
        
        % eliminate rows and columns from the Gcore
        % -- mean --
        Gcore_mean_temp(ind_row_discarded,:,:) = [];
        Gcore_mean_temp(:,ind_column_discarded,:) = [];
        Gcore_var_temp(ind_row_discarded,:,:) = [];
        Gcore_var_temp(:,ind_column_discarded,:) = [];

        
        % -- var & cor --
%         varrow_discarded = kron(ones(1,r_column),row_discarded);
%         varcolumn_discarded = kron(column_discarded,ones(1,r_row));
%         ind_var_discarded = find(varrow_discarded | varcolumn_discarded);
        
%         ind_kron_row_discarded = find(kron(row_discarded,ones(1,r_row)) | kron(ones(1,r_row),row_discarded));
%         ind_kron_column_discarded = find(kron(column_discarded,ones(1,r_column)) | kron(ones(1,r_column),column_discarded));
        
%         Gcore_var_temp(ind_var_discarded,:,:) = [];
%         Gcore_var_temp(:,ind_var_discarded,:) = [];
%         Gcore_cor_temp(ind_var_discarded,:,:) = [];
%         Gcore_cor_temp(:,ind_var_discarded,:) = [];
%         Gcore_kron_temp(ind_kron_row_discarded,:,:) = [];
%         Gcore_kron_temp(:,ind_kron_column_discarded,:) = [];
       
        % eliminate elements in Lambda, only operate on the lambda_row, since the lambda_column are needed in the next rank reduction
        % -- mean --
        lambda_row_temp(ind_row_discarded) = [];
        % -- var, actually var of lambda makes no sense in this program --
        lambda_var_row_temp(ind_row_discarded) = [];
        
        % substitute the structure with temp
        Gcore.mean{i} = Gcore_mean_temp;
        Gcore.var{i} = Gcore_var_temp;

%         Gcore.cor{i} = Gcore_cor_temp;
%         Gcore.kron{i} = Gcore_kron_temp;
        Lambda.mean{i} = lambda_row_temp;
        Lambda.var{i} = lambda_var_row_temp;
    end
    Gcore_reduced = Gcore;
    Lambda_reduced = Lambda;
end