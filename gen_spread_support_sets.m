function locs_d1=gen_spread_support_sets(K,len,space,N)
    shift_limit=N-((space+1)*K-space-1);
    if shift_limit<=0
        msg = 'Error - out of N range: Change space and K';
        error(msg)
    end
    rem=len-shift_limit;
    locs_d1=cell(shift_limit,1);
    ii=1;
    while ii<=shift_limit
         arr=[];
         for i=1:K
            arr=[arr,i*(space+1) - space + (ii-1)];
         end 
         arr(2:end-1)=randi([-1 1],1,K-2)+arr(2:end-1);
         locs_d1{ii}=arr;
         ii=ii+1; 
    end

if rem<shift_limit
   locs_d1 = [locs_d1;locs_d1(1:rem)];
else
   while length(locs_d1)<len
        locs_d1 = [locs_d1;locs_d1];
   end
end
locs_d1=locs_d1(1:len);
