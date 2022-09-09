function rse = rse_score(X,Y)
    % Get the relative square error between X and Y
    % X is the estimated tensor, and Y is the original tensor
    % rse = ||X-Y||_F  /  ||Y||_F
    if size(X)~=size(Y)
        error('The compared tensor should be with the same size')
    end
    
    rse = sqrt(sum( (X(:) - Y(:)).^2 )/sum(Y(:).^2));

end
