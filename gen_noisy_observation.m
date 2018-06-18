function [yy,x] =gen_noisy_observation(SNR,mode,k,N,idx_y)  
    % Create a time domain signal. Stream of diracs
    x = zeros(N, 1);
    if mode == 1
        x(k) = randn(size(k));
    elseif mode == 2
        x(k) = ones(size(k));
    else
        msg = 'Error - mode_spikes: Must be either 1 or 2';
        error(msg);
    end
    % Take its fourier transform    
    X_k = fft(x);
    % Noiseless samples
    y   = X_k(idx_y) / sqrt(N);
    
    % Add noise to the measured vector
    P_y   = y' * y / length(y);
    e     = randn(size(y)) + 1j * randn(size(y));
    P_e   = e' * e / length(e);
    sigma = sqrt(10^(-SNR/10) * P_y / P_e);
    e     = sigma * e;
    yy    = y + e;    
end