clear;clc;close all;
%%%%%%%%  MMSE ZF MRC %%%%%%%
w=[400 900 1600 2500 3600 4900 6400 8100 10000];
b=load('MMSEtwenty.mat');
bb=load('ZFtwenty.mat');
bbb=load('MRCtwenty.mat');
g=load('MMSEseventy.mat');
gg=load('ZFseventy.mat');
ggg=load('MRCseventy.mat');

figure(1)
muRateMMSE=zeros(1,length(w));
for s=1:length(w)
     M=w(s);
     muRateMMSE(s)=log2(1+0.14*M);
end
plot(w,muRateMMSE,'r-','DisplayName','upper bound','LineWidth',1);
hold on;
%title('Ergodic rates of an LIS-based system with randomly located devices');
xlabel('Number of antennas on each LIS unit (M)');
ylabel('Ergodic Capacity (bps/Hz)');
plot(w, b.ergodicR_k,'b--o','DisplayName','K=20 MMSE','LineWidth',1);
plot(w, bb.ergodicR_k,'g--o','DisplayName','K=20 ZF','LineWidth',1);
plot(w, bbb.ergodicR_k,'m-.o','DisplayName','K=20 MRC','LineWidth',1);
plot(w, g.ergodicR_k,'b--pentagram','DisplayName','K=70 MMSE','LineWidth',1);
plot(w, gg.ergodicR_k,'g--pentagram','DisplayName','K=70 ZF','LineWidth',1);
plot(w, ggg.ergodicR_k,'m-.pentagram','DisplayName','K=70 MRC','LineWidth',1);
grid on;
grid minor; 
hold off; 
legend show;



