function [X_k, K] = extrapolate_fourier_tls(y, N, cadzow_it, K, idx)

M          = length(y);
L          = N - M;
estimate_k = false;

if nargin < 5
    idx = 1 : M;
    if nargin < 4
        estimate_k = true;
        if nargin < 3
            cadzow_it = 0;
        end
    end
end

% Params
sv_thresh = 0.05;
idx_m     = setdiff(1:N, idx);

% Construct a square (or almost square) Toeplitz matrix
mid   = floor(M/2);
Y_toe = toeplitz(y(mid:end), y(mid:-1:1));

% Estimate K
if estimate_k
    sing_vals = svd(Y_toe); 
    sing_vals = sing_vals / sing_vals(1);
    K         = sum(sing_vals > sv_thresh);
end

% Denoise the samples
if cadzow_it > 0
    Y_toe = cadzow(Y_toe, K, cadzow_it);
    y     = [Y_toe(1,end:-1:1).'; Y_toe(2:end,1)];
end

% Get the annihilating filter and 
Y_toe   = toeplitz(y(K+1:end), y(K+1:-1:1));
[~,~,V] = svd(Y_toe);
h       = V(:,end);

% Construct the system
y1      = zeros(N, 1);
y1(idx) = y;
Y1      = toeplitz(y1, [y1(1); y1(end:-1:end-K+1)]);
E       = zeros(N, L);
for ith = 1 : L
    e_i             = zeros(N, 1);
    e_i(idx_m(ith)) = 1;
    E_i             = toeplitz(e_i, [e_i(1); e_i(end:-1:end-K+1)]);
    E(:,ith)        = E_i * h;
end
E_ext   = [Y1*h E];
[~,~,V] = svd(E_ext, 0);
ym      = V(2:end,end);

% Construct the Fourier vector
X_k = zeros(N, 1);
X_k(idx)   = y;
X_k(idx_m) = ym;

% Denoise the samples
X_toe = toeplitz(X_k(mid:end), X_k(mid:-1:1));
X_toe = cadzow(X_toe, K, cadzow_it);
X_k   = [X_toe(1,end:-1:1).'; X_toe(2:end,1)];

X_k(idx) = y;

end




function Y = cadzow(X, k, iter)
[m,n] = size(X);

Y = X;
for it = 1 : iter
    % Force the rank of the matrix to be 'k' (svd function returns the
    % singular values in decreasing order)
    [U,S,V] = svd(Y);
    Y       = U(:,1:k) * S(1:k,1:k) * V(:,1:k)';
    
    % Average the diagonal elements of Y in order to obtain a Toeplitz
    % matrix
    c = zeros(m, 1);
    r = zeros(1, n);
    c(1) = mean(diag(Y));
    r(1) = c(1);
    for ith_row = 2 : m
        c(ith_row) = mean(diag(Y,1-ith_row));
    end
    for ith_col = 2 : n
        r(ith_col) = mean(diag(Y,ith_col-1));
    end
    Y = toeplitz(c, r);
end

end
