function noise = noise_generate(Signal,SNR)
% ------------------------------------------------------
% Randomly generate noise one the given Signal, with preset SNR
% 
% ------------------Input------------------
% Signal: clean data
% SNR: noise power, defined as SNR = 10 log10(S_power/N_power)
% 
% ------------------Output------------------
% noise: generated 0-mean Gaussian noise, with the same size of Signal
% 
% XU Le, 2020
% ------------------------------------------------------

    S_power = sumsqr(Signal)/numel(Signal); % power of the signal
    N_power = S_power/(10^(SNR/10)); % power of the noise
    noise = sqrt(N_power)*randn(size(Signal));
end