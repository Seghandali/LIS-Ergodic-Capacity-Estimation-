clear;clc;close all;
%%%%%%%%  MMSE ZF %%%%%%%
w=[400 900 1600 2500 3600 4900 6400 8100 10000];
g=load('MMSEseventy.mat');
gg=load('ZFseventy.mat');
h=load('MMSEeighty.mat');
hh=load('ZFeighty.mat');

figure(1)
for s=1:length(w)
     M=w(s);
     muRateMMSE(s)=log2(1+0.14*M);
end
plot(w,muRateMMSE,'r-*','DisplayName','upper bound','LineWidth',1);
hold on;
%title('Ergodic Capacity of an LIS-based system with randomly located devices');
xlabel('Number of antennas on each LIS unit (M)');
ylabel('Ergodic Capacity (bps/Hz)');
plot(w, g.ergodicR_k,'b--^','DisplayName','K=70 MMSE','LineWidth',1);
plot(w, gg.ergodicR_k,'g--^','DisplayName','K=70 ZF','LineWidth',1);
plot(w, h.ergodicR_k,'b--o','DisplayName','K=80 MMSE','LineWidth',1);
plot(w, hh.ergodicR_k,'g--o','DisplayName','K=80 ZF','LineWidth',1);
grid on;
grid minor; 
hold off; 
legend show;



