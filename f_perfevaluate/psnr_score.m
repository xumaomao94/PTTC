function psnr = psnr_score(X,Y)
    % Get the PSNR between two images, X and Y
    % X is the estimated image, and Y is the original image
    % psnr = 10 log_10{max(Y)^2 / [||X-Y||_F^2  /  numel(Y)]}
    if size(X)~=size(Y)
        error('The compared tensor should be with the same size')
    end
    
    psnr = 10 * log10(  max(Y(:)).^2 / mean( (X(:) - Y(:)).^2 )  );

end