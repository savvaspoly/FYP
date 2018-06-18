function [x] = fri_extrapolate_fourier_tls_real(y, N, cadzow_it, K, idx)

M = length(y);
y = y * sqrt(N);

% Obtain the extrapolated Fourier transform

if nargin == 2
    
    cadzow_it  = 0;
    idx        = 1 : M;
    [y_ext, K] = extrapolate_fourier_tls(y, N, cadzow_it);
    
elseif nargin == 3
    
    idx       = 1 : M;
    y_ext     = extrapolate_fourier_tls(y, N, cadzow_it);
    
elseif nargin == 4
    
    idx       = 1 : M;
    y_ext     = extrapolate_fourier_tls(y, N, cadzow_it, K);
    
elseif nargin == 5
    
    y_ext     = extrapolate_fourier_tls(y, N, cadzow_it, K, idx);
    
end

% Construct the sampling matrix
F_M = fft(eye(N));
D   = F_M(idx,:);

% Extrapolate FT from y
xx = real(ifft(y_ext * sqrt(N)));

% Get peaks
a_x = abs(xx);
a_x = [0; eps; a_x; eps; 0];
peaks = find((a_x(3:end-2) > a_x(4:end-1)) & (a_x(3:end-2) > a_x(2:end-3)) ...
             & (a_x(2:end-3) >= a_x(1:end-4)) & (a_x(4:end-1) >= a_x(5:end)));

% Compute the sparse vector
x        = zeros(N, 1);
x(peaks) = real(D(:,peaks) \ y);
[~,idx]  = sort(abs(x));
idx_z    = idx(1:end-K);
x(idx_z) = 0;
