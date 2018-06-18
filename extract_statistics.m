%Functions needed BPDN_s, new_algorithm_s, compare_merge_s, add_noise_to_y_s
function [comp]=extract_statistics(SNR,K,spacing,reps,mode_s, mode_d,M)
% Parameters
N         = 256;
cadzow_it = 5;

% Acquisition matrix
idx_0 = 1;
idx_F = M + idx_0 - 1;
idx_y = idx_0:idx_F;
F_M   = fft(eye(N));
D_n   = F_M(idx_y,:) / sqrt(N);

% Generate support sets
if mode_d == 1
    locs_d=gen_spread_support_sets(K,reps,spacing,N);          %equally spaced+-1 randomnsess
elseif mode_d == 2
    locs_d=gen_rand_support_sets(K,reps,N);              %just random
else
    msg = 'Error - mode_dist: Must be either 1 or 2';
    error(msg);
end

comp=zeros(3,3);
for ith_k = 1 : reps
    
    disp(['SNR = ' num2str(SNR) ', K = ' num2str(K) ', iteration = ' num2str(ith_k)])
    % Sparse vector and measured samples
    k = locs_d{ith_k};
    k = k(:); %support_set
    %Create the noisy sampled signal in F domain yy. x is the original to
    %be recontructed.
    [yy,x] =gen_noisy_observation(SNR,mode_s,k,N,idx_y);
    %=============================FRINA_algorithm=========================
    % Extrapolate entire FT vector - TLS
    x_tls = fri_extrapolate_fourier_tls_real(yy, N, cadzow_it, K, idx_y);
    %=============================BPDN_algorithm=========================
    LAM=0.005;
    [x_bpdn,x_bpdn1] = BPDN_s(yy,D_n,K,LAM);
    %=============================NEW_algorithm=========================
    x_new=new_algorithm(D_n,yy,x_tls,x_bpdn,x_bpdn1,K,idx_y);  
    
    E_x = x' * x;
    %RMSE, success rates and MSE statistics
    comp(1,1)= comp(1,1)+sqrt(mean(abs((yy - D_n*x_bpdn).^2)));
    comp(1,2)= comp(1,2)+sqrt(mean(abs((yy - D_n*x_tls).^2)));
    comp(1,3)= comp(1,3)+sqrt(mean(abs((yy - D_n*x_new).^2)));

    comp(2,1)= comp(2,1)+suc_rate(x,x_bpdn);
    comp(2,2)= comp(2,2)+suc_rate(x,x_tls);
    comp(2,3)= comp(2,3)+suc_rate(x,x_new);
  
    comp(3,1)= comp(3,1)+(x-x_bpdn)' * (x-x_bpdn) / E_x;
    comp(3,2)= comp(3,2)+(x-x_tls)' * (x-x_tls) / E_x;
    comp(3,3)= comp(3,3)+(x-x_new)' * (x-x_new) / E_x;
    
end
comp=comp/reps;
end
