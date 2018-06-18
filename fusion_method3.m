function [x_fi,D3]=fusion_method3(D_n,yy,x_tls,x_bpdn,K,w)
    x_fi=zeros(size(x_bpdn));
    same_loc = find ( (x_tls~=0) & (x_bpdn~=0) ); %common support set
    x_fi(same_loc)= x_tls(same_loc)*w + x_bpdn(same_loc)*(1-w);%weighted average. w=0.5 usually
    itr=K-sum((same_loc~=0)); %The remaining sparsity to recover
    
    LL=sum (x_fi~=0);
    dif_loc = find (bitxor ( (x_tls~=0) , (x_bpdn~=0) )); %joint support set esxluding the common support set
    D_nf=D_n(:,same_loc);
    D_nf=[D_nf,D_n(:,dif_loc)]; %sub matrix based on the given indices
    select_vals =x_tls(dif_loc)+ x_bpdn(dif_loc);
    n=size(D_nf,2);
    cvx_begin quiet
        variable x_m(n) % declares x to be an optimization variable of dimension n.
        minimize( norm(yy-D_nf*x_m) );
    cvx_end 
    x_m(1:LL)=[];
    %Create estimate x for the Least RMSE approachs
    [~,I] = sort(abs(x_m),'descend');
    x_fi(dif_loc(I(1:itr)))=select_vals(I(1:itr)); %selecting itr remaining magnitudes
    %The estimated support set of this fusion method
    D3=same_loc;
    D3=[D3;dif_loc(I(1:itr))];
end