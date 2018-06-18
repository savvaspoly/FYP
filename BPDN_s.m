function [x_bpdn,x_bpdn1] = BPDN_s(yy,D_n,K,LAM) 
        N=size(D_n,2);
        cvx_begin quiet
            variable x_bpdn(N) % declares x to be an optimization variable of dimension n.
            minimize( 0.5*norm(yy-D_n*x_bpdn)+ LAM*norm(x_bpdn,1) );
        cvx_end 
        x_bpdn1=x_bpdn;
        [~,idx] = sort(abs(x_bpdn));
        zero_i  = idx(1:end-K);
        x_bpdn(zero_i) = 0;  %The x_bpnd with K sparsity
        
        [~,idx] = sort(abs(x_bpdn1));
        zero_i  = idx(1:end-(K+2)); 
        x_bpdn1(zero_i) = 0; %The x_bpnd with K+2 sparsity
end