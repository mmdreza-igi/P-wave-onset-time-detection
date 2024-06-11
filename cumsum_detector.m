clc
clear
tic

[D,T0,H] = rdsac('2020-02-15-2144-00SSNQR__001_BH_Z_SNQR__BH_Z_.SAC');

%D=D(1:10000);
DN = cumsum(D);
Window_s=50;
OverLap=5;
DNn=[];
for i=1:(length(DN)-Window_s)/OverLap
    a(i,:)=DN(i*OverLap:i*OverLap+Window_s);
end
for k=1:i-1
    b(k,:)=a(k+1,:)-a(k,:);
    f(k)=mean(b(k,:));
end
B=b';

%% findpeaks

% for j=1:i-1
%     y=findpeaks(b(j,:));
%     num_p(j)=length(y);
% end

%% diffrential
interval=5;
sum_x=zeros(length(b),1);
for j=1:i-1
    for k=1:interval:Window_s
        X=b(j,k:k+interval-1);
        min_x=min(X);
        max_x=max(X);
        diffren=abs(max_x-min_x)/interval;
        sum_x(j,1)=sum_x(j,1)+diffren;
    end
end
        

%% picking P_wave

%% 1
% for j=2:i-7
%     if sum_x(j+1)>max(sum_x(1:j)) && sum_x(j+2)>max(sum_x(1:j+1)) && sum_x(j+3)>max(sum_x(1:j+2)) &&  sum_x(j+4)>max(sum_x(1:j+3))   
%         break
%     end
% end

%% 2
% vatar_m = sqrt(mean(sum_x(1:end))^2 + length(b)^2);
% zavie_m = asind(mean(sum_x(1:end))/vatar_m);
% for j=1:i-1
%      vatar=sqrt(sum_x(j)^2+length(b)^2);
%     zavie(j)=asind(sum_x(j)/vatar);
% end
% for j=1:i
%     if zavie(j) > 0.7*zavie_m
%         break
%     end
% end
% first_pick = j-1;
% second_pick = first_pick*OverLap + Window_s;

%% 3

lta_w = 150;
sta_w = 20;

for j=1:i-lta_w-1
    sta = mean(sum_x(j+lta_w-sta_w:j+lta_w));
    lta = mean(sum_x(j:j+lta_w));
    sta_lta(j+lta_w-sta_w) = sta/lta;
end

pick = (find(sta_lta==max(sta_lta))-sta_w/2)*OverLap+Window_s;


%% plot

figure;plot(sum_x)
hold on

%xline(find(sta_lta==max(sta_lta)),'--r')

figure;plot(D)
hold on
lg1=xline(pick,'r');
hold on 
lg2=xline(5411,'k');
lgd=legend([lg1,lg2],{'CUMSUM-STA/LTA(5370)','manual(5411)'});
toc