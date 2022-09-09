function A = tt2full(Gcore,Size_A)
%% a much faster one adopting the idea of TTSVD
A = 1;
for i = 1:length(Size_A)
%     gcore_temp = Gcore.mean{i}; % of size r(i)*n(i)*r(i+1)
    Size_gcore = size(Gcore.mean{i});
    A = A*reshape(permute(Gcore.mean{i},[1,3,2]),Size_gcore(1),[]);
    A = reshape(A,size(A,1)*Size_gcore(3),[]);
end
A = reshape(A,Size_A);
end