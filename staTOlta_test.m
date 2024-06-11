clc
clear
tic
%opening sac data 
[D,T0,H]=rdsac('1997-02-28-1702-14S-AZR___001_SP_Z_AZR___SP_Z_.SAC'); 

t = T0 + (H.B + (0:H.DELTA:(H.NPTS - 1)*H.DELTA)')/86400;
% D(length(D)+1)=0;
% t(length(D))=t(length(D)-1)+0.020;
%D=bandpass(D,[1 6],100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% remove DC bias (r) and calculation characteristic function (e)

c1=2;c2=1;
fd=zeros(1,length(D));
r=zeros(1,length(D));
dr=zeros(1,length(D));
e=zeros(1,length(D));
for i=1:length(D)
    if i==1
        fd(i)=D(1);
        r(i)=c1*D(1)+fd(i);
        dr(i)=c2*fd(i);
        e(i)=r(i)^2+dr(i)^2;
    else
        fd(i)=D(i)-D(i-1);
        r(i)=c1*D(i-1)+fd(i);
        dr(i)=c2*fd(i);
        e(i)=r(i)^2+dr(i)^2;
    end
    
end
r=r';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot waveform and calculate and plot zer-crossing 
% figure;plot(t,r);
% xlim = [min(t),max(t)];
% set(gca,'XLim',xlim)
% datetick('x','keeplimits');
% xlabel(sprintf('%s to %s',datestr(xlim(1)),datestr(xlim(2))))
% yline(sum(r)/length(r),'linewidth',1);

%  hold on

dc=sum(r)/length(r);
%zc=zeros(length(r),1);
ti=1;
for i=2:length(r)
    if r(i-1)*r(i)==abs(r(i-1)*r(i))
    else
        alfa=(r(i)-r(i-1))/(t(i)-t(i-1));
        beta=r(i)-alfa*t(i);
        pic(ti)=i;
        zc(ti,1)=(dc-beta)/alfa;
        zc(ti,2)=t(i);
%         plot(zc(ti),0,"r+")
        ti=ti+1;
    end
end

% hold on

 for i=2:length(pic)
     mr=r(pic(i-1));
     for j=pic(i-1):pic(i)
         if pic(i)==length(D)
             j=pic(i-1);
         end
             
         if abs(mr)<=abs(r(j+1))
             mr=r(j+1);
             mt=t(j+1);
         else 
             mr=r(j+1);
             mt=t(pic(i-1));
             
         end
     end
     maxr(i)=mr;
     maxt(i)=mt;
     
 end
 
% plot(maxt,maxr,'gO')
      
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot characteristic function of waveform

% figure;plot(e')
% xlim = [min(t),max(t)];
% set(gca,'XLim',xlim)
% datetick('x','keeplimits');
% xlabel(sprintf('%s to %s',datestr(xlim(1)),datestr(xlim(2))))
% 
% rt=datestr(t,'MM:ss:fff');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate sta and lta an find where sta/lta exceed preset tereshold and
% save data we need
e=e';
ti=1;
sta_window_length=10;
lta_window_length=100;
sta=zeros(length(D)-sta_window_length+1,1);
lta=zeros(length(D)-lta_window_length+1,1);
for i=1:length(e)
    if i<sta_window_length+1
        sta(1)=e(i)/sta_window_length+sta(1);
    else
        ti=ti+1;
        sta(ti)=sta(ti-1)+0.7*(e(i)-sta(ti-1));
    end
end
ti=1;
for i=1:length(D)
    if i<=lta_window_length
        lta(1)=e(i)/lta_window_length+lta(1);
    else
        ti=ti+1;
        lta(ti)=lta(ti-1)+0.05*(e(i)-lta(ti-1));
    end
end
%data=zeros(30000,5);
%data_t=cell(30000,1);
ti=1;
for i=1:length(sta)-lta_window_length
    gm(i)=lta(i)*5;
    if sta(lta_window_length-1+i)/lta(i)>40      
        
        
        gama(ti)=sta(lta_window_length-1+i)/lta(i);
        tgama(ti)=t(lta_window_length+9+i);
      
       % nsta(ti)=sta(49+i);
       % data(i,1)=t(104+i,1);
       % data_t(i,1)={convertCharsToStrings(rt(104+i,:))};
       % data(i,2)=r(104+i);
       % data(i,3)=dr(104+i);
        ti=1+ti;
    else
        %data(i,1)=0;
        %data(i,2)=0;
        %data(i,3)=0;
        %data(i,4)=0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;plot(r) 
yline(0)
hold on
plot(tgama(1),0,'ro')

%%%% simpling frequently sta

for i=2:length(tgama)
    if tgama(i)-tgama(i-1)<t(6)-t(1)
        tgama(i)=tgama(i-1);
    end
end
%%%% tereshold for zero-crssing
for i=2:length(tgama)
    npick(i)=1;
    for j=1:length(zc)
        if zc(j,2)>tgama(i-1) && zc(j,2) <tgama(i)
            npick(i)=1+npick(i);
        end

    end
    if npick(i)-npick(i-1)<15
        tgama(i)=tgama(i-1);
    end
end

for i=1:length(t)
    if tgama(1)==t(i)
        t_main=i;
    else
    end
end
    p_pick=tgama(1);
for i=t_main:-1:t_main-100
   if abs(D(i-1))>abs(D(i))
        p_pick=t(i);
        P_sample=i;
        break
   else
    end
end
% figure;plot(t,r) 
% yline(sum(r)/length(r))
% hold on
% plot(tgama,0,'ko')
% hold on
% % PP=datestr(p_pick,'MM:ss:fff');
% % PP=convertCharsToStrings(PP);
% % PP=PP+' -> mm:ss:ms';
% xline(p_pick)
% %text(p_pick,-5.5*10^3,PP)
figure;plot(D);
hold on
plot(P_sample,0,'ro')
hold on
lg1=xline(P_sample,'r');
hold on 
lg2=xline(6243,'k');
lgd=legend([lg1,lg2],{'sta/lta(3353)','manual(6243)'});
toc