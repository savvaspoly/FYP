function locs_d1=gen_rand_support_sets(K,len,N)
    locs_d1=cell(len,1);
    ii=1;
    while ii<=len
         arr=[];
         arr=sort(randperm(N,K),'ascend');
         locs_d1{ii}=arr;
         ii=ii+1; 
    end
end
