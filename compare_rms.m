% Compare x_tls and x_bpdn and add the same locations
function [x_b]=compare_rms(x,x_tls,x_bpdn,N,w)

    x_b=zeros(N,1);
    x1= ((x~=0) & (x_bpdn~=0) & (x_tls~=0) ); %same locs
    x2= ((x~=0) & (x_bpdn~=0) ) | ((x~=0) & (x_tls~=0) ); %desired ones
    diff_loc =find (bitxor ( (x1~=0) , (x2~=0) ));
    same_loc =find (x1);

    x_b(same_loc)=x_bpdn(same_loc)*(1-w) + x(same_loc)*w;
    x_b(diff_loc)=x_bpdn(diff_loc) + x_tls(diff_loc);
end