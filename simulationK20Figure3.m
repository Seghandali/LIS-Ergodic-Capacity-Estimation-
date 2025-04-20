%  Ergodic Rate, Monte Carlo Simulations compare with Analytical results-ZF
% The execution time depends on the value of K,for example with iteration t=100,for K=10, it is approximately equal to 12 hours,and for K=70, it is approximately equal to 90 hours.
%%%%% K=20, randomly located devices, combination of LOS and NLOS interference  %%%%%
clear; clc; close all;
%tic
w=[400 900 1600 2500 3600 4900 6400 8100 10000];
ergodicR_k=zeros(1,length(w));
varianc=zeros(1,length(w));
K=20;     % Number of devices   
L=0.25;   % Length of LIS unit : 2L=0.5m
lambda=0.1; % Carrier frequency : 3GHz
snr_dB = 3; % Uplink target SNR
snr = 10.^(snr_dB./10);
p_k =snr; p_j =snr; 
B_PL=3.7;  %  path loss exponent
for s=1:length(w)
M=w(s); % Number of antennas in each LIS unit
deltaL=(2*L)/sqrt(M); % antenna spacing
P=M/2; % represents the number of dominant paths among all NLOS paths 
N=100;
gamma_k=zeros(1,N);
R_k=zeros(1,N);
variance=zeros(1,N);
for t=1:N % for each iteration of MONTE CARLO simulation
% Finding location of device j & k(target device) : x device & y device , z device 
Xdevice=zeros(1,K);
Ydevice=zeros(1,K);
Zdevice=zeros(1,K);
x=unifrnd(-2,2,1,K-1);
y=unifrnd(-2,2,1,K-1);
z=unifrnd(0.5,2,1,K-1);
for j=1:K-1 % randomly and uniformly deploy the devices in a 4m*4m*2m space.
    Xdevice(1,j)=x(1,j);
    Ydevice(1,j)=y(1,j);
    Zdevice(1,j)=z(1,j);
end
% target device is located at (0,0,1)
Xdevice(1,K)=0;
Ydevice(1,K)=0;
Zdevice(1,K)=1;
% Finding location of antenna m of LIS unit k : xLIS_km & yLIS_km , zLIS_km = 0
n=sqrt(M);
XLIS=zeros(1,M);
YLIS=zeros(1,M);
for m=0:M-1 % Each LIS unit has M antennas, antenna spacing of deltaL in a rectangular lattice form with an area limited to 2L x 2L
    %[q, r] = quorem(sym(i), sym(n));
    q=floor(m/n);
    r=m-n*q;
    y_i=(q-(n-1)/2)*(-1)*deltaL;
    x_i=(r-(n-1)/2)*deltaL;
    XLIS(1,m+1)=x_i;
    YLIS(1,m+1)=y_i;
end
% ***************************Channel Model*******************************
   H=zeros(M,K);  % The channel H is the M × K matrix between device k or j and antenna m of LIS unit k
% The desired channel h_kk between device k and LIS unit k:LOS path, target device is located at (0,0,1): j=K
h_kk=zeros(M,1);
alfaL_kkm=zeros(M,1);
lL_kkm=zeros(M,1);
h_kkm=zeros(M,1);
betaL_kkm=zeros(M,1);
for m=1:M
    d_kkm= sqrt((0-XLIS(1,m))^2+(0-YLIS(1,m))^2+Zdevice(1,K)^2);% distance between device k(target device) and antenna m of LIS unit k
    alfaL_kkm(m,1)=sqrt(Zdevice(1,K)/d_kkm); % The antenna gain 
    lL_kkm(m,1)=1/sqrt(4*pi*d_kkm^2); % The free space path loss attenuation
    h_kkm(m,1)=exp((-1)*1i*2*pi*d_kkm/lambda); % h_kkm is the LOS channel state between device k and antenna m of LIS unit k
    betaL_kkm(m,1)=alfaL_kkm(m,1)*lL_kkm(m,1); % LOS channel gain between device k and antenna m of LIS unit k
    h_kk(m,1)=betaL_kkm(m,1)*h_kkm(m,1); % desired channel h_kk(m,1) between device k and antenna m of LIS unit k
    H(m,K)= h_kk(m,1);
end
%%%%% The interference channel h_jk between device j and LIS unit k:combination of LOS and NLOS
% deterministic LOS component from device j to LIS unit k given by hL_jk
hL_jk=zeros(M,K-1);
alfaL_jkm=zeros(M,K-1);
lL_jkm=zeros(M,K-1);
h_jkm=zeros(M,K-1);
betaL_jm=zeros(M,K-1);
for j=1:K-1 % j = 1: K except target device,td
    for m=1:M
        d_jkm= sqrt((Xdevice(1,j)-XLIS(1,m))^2+(Ydevice(1,j)-YLIS(1,m))^2+(Zdevice(1,j))^2);% distance between device j and antenna m of LIS unit k
        alfaL_jkm(m,j)=sqrt(Zdevice(1,j)/d_jkm); % The antenna gain 
        lL_jkm(m,j)=1/sqrt(4*pi*d_jkm^2); % The free space path loss attenuation
        h_jkm(m,j)=exp((-1)*1i*2*pi*d_jkm/lambda); % LOS channel state between device j and antenna m of LIS unit k 
        betaL_jm(m,j)=alfaL_jkm(m,j)*lL_jkm(m,j); % LOS channel gain between device j and antenna m of LIS unit k
        hL_jk(m,j)=betaL_jm(m,j)*h_jkm(m,j); % deterministic LOS component of interference channel from device j to antenna m of LIS unit k
    end
end
% correlated NLOS component from device j to LIS unit k given by hNL_jk
hNL_jk=zeros(M,K-1);

for j=1:K-1 % j = 1: K except target device,td
    lNL_jkm=zeros(M,1);
    dv=zeros(sqrt(M),1);
    dh=zeros(sqrt(M),1);
    d=zeros(M,1);
    D_jk=zeros(M,P);
    R_jk=zeros(M,P);
    for m=1:M
        djkm= sqrt((Xdevice(1,j)-XLIS(1,m))^2+(Ydevice(1,j)-YLIS(1,m))^2+(Zdevice(1,j))^2);% distance between device j and antenna m of LIS unit k
        lNL_jkm(m,1)=djkm^(-B_PL/2); % The path loss attenuation factors 
    end
    LNL_jk=diag(lNL_jkm);
    g_jk=sqrt(1./2).*(randn(P,1) + 1i.*randn(P,1)); % An independent fast-fading channel vector
    for p=1:P 
        tetav_jkp=unifrnd(-pi/2,pi/2);
        tetah_jkp=unifrnd(-pi/2,pi/2);
        alfaNL_jkp=sqrt(cos(tetav_jkp)*cos(tetah_jkp));
        fiv_jkp=sin(tetav_jkp);
        fih_jkp=sin(tetah_jkp)*cos(tetav_jkp);
        for m=1:sqrt(M)
            dv(m,1)=exp((m-1)*1i*(2*pi*deltaL/lambda)*fiv_jkp); % (3)
            dh(m,1)=exp((m-1)*1i*(2*pi*deltaL/lambda)*fih_jkp); % (4)
        end
        d=(1/sqrt(M))*(kron(dv,dh)); % (2)
        for m=1:M
            D_jk(m,p)=alfaNL_jkp*d(m,1);
            R_jk(m,p)=LNL_jk(m,m)*D_jk(m,p);  % (R_jk)^1/2 : the deterministic correlation matrix from device j to LIS unit k
        end
    end
    Hn=R_jk*g_jk;
    for m=1:M
        hNL_jk(m,j)=Hn(m,1); % The correlated NLOS component of interference channel between device j and LIS unit k
    end
end
%toc
% The interference channel h_jk between device j and LIS unit k:combination of LOS and NLOS
h_jk=zeros(M,K-1);
d_c=10;
for j=1:K-1
    d_jk= sqrt((Xdevice(1,j))^2+(Ydevice(1,j))^2+(Zdevice(1,j))^2);% distance between device j and the center of LIS unit k
    ric_dB=13 - 0.03 * d_jk; % Rician factor (k[dB])
    ric=10.^(ric_dB./10);
    % probability of LOS path
    if (0<d_jk)&&(d_jk<10)
        PLOS_jk=(d_c-d_jk)/d_c; % (48)
    else
        PLOS_jk=0;        
    end
    a=rand;
    for m=1:M
        if a<PLOS_jk
            h_jk(m,j)=sqrt(ric/(ric+1))* hL_jk(m,j);
            H(m,j)= h_jk(m,j);
        else
            h_jk(m,j)=sqrt(1/(ric+1))* hNL_jk(m,j);
            H(m,j)= h_jk(m,j);
        end
    end
end
 zf=inv(H'*H);
 I_k=zf(K,K);
 gamma_k(t)=p_k/I_k; 
% uplink data rate 
 R_k(t)=log2(1+gamma_k(t)); 
 variance(t)= var(R_k(t));
end
% averaging for monte-carlos
ergodicR_k(s)=mean(R_k);
varianc(s)=var(R_k);
end

% Distance from Favorable Propagation
sumC=0; 
for j=1:K
    sumC=sumC+log2(1+snr*(norm(H(:,j))^2));
end
I=eye(K);
upperB=log2(det(I+snr*H'*H));
deltaC=(sumC-upperB)/upperB 
mu=mean(gamma_k);
v=var(gamma_k);  
muRate=log2(1+mu)-v/(2*(1+mu)^2); 
mean(gamma_k) /M;

%plots
figure(1)
muRateZF=zeros(1,length(w));
for s=1:length(w)
    M=w(s);
    %muRateZF(s)=log2(1+0.14*M);
    gamma_k_estimation = p_k*((M/(4*pi*L^2))*atan(L^2/(Zdevice(1,K)*sqrt(2*L^2+Zdevice(1,K)^2))));
    gamma_k_estimation / M;
    muRateZF(s)=log2(1+gamma_k_estimation);
end
plot(w,muRateZF,'r-','DisplayName','LIS subjects to the channel hardening effect (upper bound)');
grid on;
plot(w,ergodicR_k,'b-O','DisplayName','K=20, ZF simulation');
%title('Ergodic Capacity of an LIS-based system with randomly located devices');
xlabel('Number of antennas on each LIS unit(M)');
ylabel('Ergodic Capacity (bps/Hz)');
legend show;








    
    
    






    
    
    