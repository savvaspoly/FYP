function [prcntg] = suc_rate(x,x_prd)
     p=sum(x~=0 & x_prd~=0);
     k=sum(x~=0);
     prcntg=p/k*100;
end

