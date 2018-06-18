function per=extract_graphs(mode_s,mode_d,reps,SNR)
%Effect of varying K
variable_test=2:2:32;
size_p=3;
per=zeros(3,size_p,length(variable_test));
for i=1:length(variable_test)
    max_d=floor((256/variable_test(i))-1);
    per(:,:,i)=extract_statistics(SNR,variable_test(i),max_d,reps,mode_s,mode_d,64);
end

colours =["b","r","--k","--c","--m","--y","--g"];
set=[1,2,3];
size_p=length(set);
names ={'BSS-Spread','BSS','GSS-Spread','GSS'};
leg_names ={'x-BPDN','x-FRINA','x-NEW'};
if ((mode_s==2)&&(mode_d==1))
    k=1;
elseif(mode_s==2) &&(mode_d==2)
    k=2;
elseif(mode_s==1) &&(mode_d==1)
    k=3;
elseif(mode_s==1) &&(mode_d==2)
    k=4;
else
    msg = 'Error - mode_dist and mode_spikes: Must be either 1 or 2';
    error(msg);
end

%Plot the average RMSE in relation to observation  for each K
figure;
for i=1:size_p
    hold on
    plot(variable_test,squeeze(per(1,set(i),:)),colours(i), 'Marker', 'x','MarkerEdgeColor','k');
    grid on;
    grid minor;
end
legend(leg_names,'Location','southeast');
title (['RMSE vs K, SNR=',num2str(SNR),', ' ,names{k}]);
ylabel('RMSE');
xlabel('K');

%Plot the average % success rates  for each K
figure; 
for i=1:(size_p)
    hold on
    plot(variable_test,squeeze(per(2,set(i),:)),colours(i), 'Marker', 'x','MarkerEdgeColor','k');
    grid on;
    grid minor;
end
legend(leg_names,'Location','southwest');
title (['Success vs K, SNR=',num2str(SNR),', ' ,names{k}]);
ylabel('Success%');
xlabel('K');

%Plot the average MSE for each K
figure;
for i=1:size_p
    hold on
    plot(variable_test,squeeze(per(3,set(i),:)),colours(i), 'Marker', 'x','MarkerEdgeColor','k');
    grid on;
    grid minor;
end
legend(leg_names,'Location','southwest');
title (['MSE vs K, SNR=',num2str(SNR),', ' ,names{k}]);
ylabel('MSE');
xlabel('K');

end