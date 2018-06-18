function [x_fi,D1]=fusion_method2(D_n,yy,x_tls,x_bpdn,K,w)
    x_fi=zeros(size(x_bpdn));
    same_loc = find ( (x_tls~=0) & (x_bpdn~=0) );  %common support set
    x_fi(same_loc)= x_tls(same_loc)*w + x_bpdn(same_loc)*(1-w); %weighted average according to w=between 0 and 1. Usually w=0.5
    itr=K-sum((same_loc~=0)); %The remaining sparsity to recover
    
    dif_loc = find (bitxor ( (x_tls~=0) , (x_bpdn~=0) ));  %joint support set esxluding the common support set
    D_nf=D_n(:,dif_loc); %Sub-matrix based on dif_loc indices
    select_vals = x_tls(dif_loc)+ x_bpdn(dif_loc); 
    n=size(D_nf,2);
    yyy=yy-D_n*x_fi; %substract from the observation the common locations' constribution
    cvx_begin quiet
        variable x_m(n) % declares x to be an optimization variable of dimension n.
        minimize( norm(yyy-D_nf*x_m) );
    cvx_end 
    %Create estimate x for the Least RMSE approach
    [~,I] = sort(abs(x_m),'descend');
    x_fi(dif_loc(I(1:itr)))=select_vals(I(1:itr));   %selecting itr remaining magnitudes
    %The estimated support set of this fusion method
    D1=same_loc;
    D1=[D1;dif_loc(I(1:itr))];
end