clc
clear
tic
%% input data
[D1,T0,H]=rdsac('2020-02-15-2144-00SSNQR__001_BH_Z_SNQR__BH_Z_.SAC');
%D1=importdata('2022-06-19 0608 18_MDJ.txt');
P=7500;
%% some constent parameters
OrderNum=4;
n0=P-1000;
n1=P+1000;
t=0;
p=10;
tic
%% main processing part

for window=n0:p:n1    %% moving in a window that we think p onset is in
    
    if window==n0
        for j=1:window      %% building first AR model for noise part
            for k=1:OrderNum
                x(j,k)=D1(OrderNum+j-k);
            end
            y(j,1)=D1(OrderNum+j);
        end
        x(:,OrderNum+1)=y;
        [~,R]=qr(x);
        for i=1:OrderNum
            er(i)=(1/(n0-OrderNum))*sum(R(i:end,end).^2);
            aic(i)=(n0-OrderNum)*log10(er(i))+2*(i+1);
        end
    else
        for j=window-p+1:window      %% building first AR model for noise part
            for k=1:OrderNum
                x1(j-window+p,k)=D1(OrderNum+j-k);
            end
            y1(j-window+p,1)=D1(OrderNum+j);
        end
        x1(:,OrderNum+1)=y1;
        
        X=[R;x1];
        [~,R]=qr(X);
        x1=[];
        y1=[];
        
        for i=1:OrderNum
            er(i)=(1/(n0-OrderNum+p*t))*sum(R(i:end,end).^2);
            aic(i)=(n0-OrderNum+p*t)*log10(er(i))+2*(i+1);
        end
    end
    t=t+1;
    samp(window,1)=t-1;
    AIC(t)=min(aic);
    x1=[];
    y1=[];
end
t=0;
for window=n1:-p:n0
    if window==n1
        for j1=window:n1+1500 %%n1+1000 %% building second AR model for signal part
            for k1=1:OrderNum
                x2(j1-window+1,k1)=D1(j1-k1);
            end
            y2(j1-window+1,1)=D1(j1);
        end
        x2(:,OrderNum+1)=y2;
        [~,R1]=qr(x2);
        
        for i1=1:OrderNum
            er1(i1)=(1/(n0-OrderNum))*sum(R1(i1:end,end).^2);
            aic1(i1)=(n0-OrderNum)*log10(er1(i1))+2*(i1+1);
        end
    else
        for j1=window+1:window+p
            for k1=1:OrderNum
                x3(-j1+window+p+1,k1)=D1(j1-k1);
            end
            y3(-j1+window+p+1,1)=D1(j1);
        end
        x3(:,OrderNum+1)=y3;
        X3=[R;x3];
        [~,R1]=qr(X3);
        
        x3=[];
        y3=[];
        
        for i1=1:OrderNum
            er1(i1)=(1/(n0-OrderNum+p*t))*sum(R1(i1:end,end).^2);
            aic1(i1)=(n0-OrderNum+p*t)*log10(er1(i1))+2*(i1+1);
        end
    end
    t=t+1;
    x3=[];
    y3=[];
    AIC1(t)=min(aic1);
end
AIC1(1)=AIC1(2);
AICC=flip(AIC1)+AIC;
figure;plot(AICC)

%% pick p wave
min_AIC=min(AICC);
A=find(AICC==min_AIC);
B=find(samp==A);
figure;plot(D1)
toc
hold on
lg1=xline(B,'r--');
hold on
lg2=xline(7830 ,'k--');
lgd=legend([lg1,lg2],{'AR-AIC(6760)','manual(7830)'});