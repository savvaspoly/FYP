function [x_fi,x_fi2]=fusion_method1_ext(D_n,yy,A,C,K)
    x_fi=zeros(size(C,1),1);
    all_locs = unique(A(:));
    same_locs=all_locs((size(A,2)-1)<histc(A(:),all_locs));
    diff_locs = setdiff(all_locs,same_locs);
    
    B=(C&1)>0;
    x_w=sum(B,2);
    x_sum=sum(C,2);   
    loc=[same_locs;diff_locs];
    D_nf=D_n(:,same_locs);
    D_nf=[D_nf,D_n(:,diff_locs)];
    select_vals=x_sum(same_locs)./x_w(same_locs);
    select_vals=[select_vals;x_sum(diff_locs)./x_w(diff_locs)];
    n=size(D_nf,2);
    cvx_begin quiet
        %cvx_precision high
        variable x_m(n) % declares x to be an optimization variable of dimension n.
        minimize( norm(yy-D_nf*x_m) );
    cvx_end
   [~,I] = sort(abs(x_m),'descend');
   x_fi2=x_fi;
   x_fi(loc(I(1:(K))))=select_vals(I(1:(K)));
   x_fi2(loc(I(1:K)))=x_m(I(1:K));
end