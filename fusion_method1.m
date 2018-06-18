function [x_fi,D2]=fusion_method1(D_n,yy,x_tls,x_bpdn,K,w)
    diff_locs = find (bitxor ( (x_tls~=0) , (x_bpdn~=0) )); %joint support set esxluding the common support set
    same_locs = find (bitand ( (x_tls~=0) , (x_bpdn~=0) )); %common support set
    loc=[same_locs;diff_locs];
    D_nf=D_n(:,same_locs);
    D_nf=[D_nf,D_n(:,diff_locs)]; %sub-matrix based on same and different locations' indices
    select_val = x_tls(same_locs)*w + x_bpdn(same_locs)*(1-w); %weighted average. w=0.5 usually
    select_val=[select_val; x_tls(diff_locs) + x_bpdn(diff_locs)];
    n=size(D_nf,2);
    cvx_begin quiet
        variable x_m(n) % declares x to be an optimization variable of dimension n.
        minimize( norm(yy-D_nf*x_m) );
    cvx_end 
    [~,I] = sort(abs(x_m),'descend');
    %Create estimate x for the Least RMSE approach
    x_fi=zeros(size(x_tls));
    x_fi(loc(I(1:K)))=select_val(I(1:K)); %selecting the larger K magnitudes
    %The estimated support set of this fusion method
    D2=loc(I(1:K));
    
end