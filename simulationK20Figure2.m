%Ergodic Rate, Monte Carlo Simulations compare with Analytical results-MRC
%%%%% K=20, randomly located devices, combination of LOS and NLOS interference  %%%%%
clear; clc; close all;
w=[400 900 1600 2500 3600 4900 6400 8100 10000];
ergodicR_k=zeros(1,length(w));
varianc=zeros(1,length(w));
K=20;     % Number of devices   
L=0.25;   % Length of LIS unit : 2L=0.5m
lambda=0.1; % Carrier frequency : 3GHz
snr_dB = 3; % Uplink target SNR
snr = 10.^(snr_dB./10);
p_k =snr; p_j =snr; 
T_k = 0.5; % Channel imperfectness
B_PL=3.7;  %  path loss exponent
delta = 1; % Hardware impairments
for s=1:length(w)
M=w(s);
deltaL=(2*L)/sqrt(M); % antenna spacing
P=M/2; % represents the number of dominant paths among all NLOS paths 
N=100;
R_k=zeros(1,N);
variance=zeros(1,N);
for t=1:N % for each iteration of MONTE CARLO simulation    
e_km_var=sqrt(1./2).*(randn(M,1) + 1i.*randn(M,1));% Estimation error vector independent random variables
c_k=sqrt(1./2).*(randn(M,1) + 1i.*randn(M,1)); % Hardware impairments at LIS unit k 

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
    q=floor(m/n);
    r=m-n*q;
    y_i=(q-(n-1)/2)*(-1)*deltaL;
    x_i=(r-(n-1)/2)*deltaL;
    XLIS(1,m+1)=x_i;
    YLIS(1,m+1)=y_i;
end
% The desired channel h_kk between device k and LIS unit k:LOS path, target device is located at (0,0,1): j=K
h_kk=zeros(M,1);
e_k=zeros(M,1);
f_k=zeros(M,1);
alfaL_kkm=zeros(M,1);
lL_kkm=zeros(M,1);
h_kkm=zeros(M,1);
betaL_kkm=zeros(M,1);
for m=1:M
    d_kkm= sqrt((0-XLIS(1,m))^2+(0-YLIS(1,m))^2+1);% distance between device k(target device) and antenna m of LIS unit k
    alfaL_kkm(m,1)=sqrt(1/d_kkm); % The antenna gain 
    lL_kkm(m,1)=1/sqrt(4*pi*d_kkm^2); % The free space path loss attenuation
    h_kkm(m,1)=exp((-1)*1i*2*pi*d_kkm/lambda); % h_kkm is the LOS channel state between device k and antenna m of LIS unit k
    betaL_kkm(m,1)=alfaL_kkm(m,1)*lL_kkm(m,1); % LOS channel gain between device k and antenna m of LIS unit k
    h_kk(m,1)=betaL_kkm(m,1)*h_kkm(m,1); % desired channel h_kk(m,1) between device k and antenna m of LIS unit k
    e_k(m,1)=betaL_kkm(m,1)*e_km_var(m,1); % Estimation error vector
    f_k(m,1)=sqrt(1-T_k^2)*h_kk(m,1)+T_k*e_k(m,1); % MF receiver Under the imperfect CSI results from an least square estimator
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
d_jk=zeros(1,K-1);
ric_dB=zeros(1,K-1);
ric=zeros(1,K-1);
hNL_jk=zeros(M,K-1);
S28=zeros(1,K-1);
S29=zeros(1,K-1);
s751=zeros(M,K-1);
tic
for j=1:K-1 % j = 1: K except target device,td
    lNL_jkm=zeros(M,1);
    dv=zeros(sqrt(M),1);
    dh=zeros(sqrt(M),1);
    d=zeros(M,1);
    D_jk=zeros(M,P);
    R_jk=zeros(M,P);
    d_jk(1,j)= sqrt((Xdevice(1,j))^2+(Ydevice(1,j))^2+(Zdevice(1,j))^2);% distance between device j and the center of LIS unit k
    ric_dB(1,j)=13 - 0.03 * d_jk(1,j); % Rician factor (k[dB])
    ric(1,j)=10.^(ric_dB(1,j)./10);
    for m=1:M
        djkm= sqrt((Xdevice(1,j)-XLIS(1,m))^2+(Ydevice(1,j)-YLIS(1,m))^2+(Zdevice(1,j))^2);% distance between device j and antenna m of LIS unit k
        lNL_jkm(m,1)=djkm^(-B_PL/2); % The path loss attenuation factors 
    end
    LNL_jk=diag(lNL_jkm);
    g_jk=sqrt(1./2).*(randn(P,1) + 1i.*randn(P,1)); % An independent fast-fading channel vector
    s28=0;
    s29=0;
    for p=1:P 
        tetav_jkp=unifrnd(-pi/2,pi/2);
        tetah_jkp=unifrnd(-pi/2,pi/2);
        alfaNL_jkp=sqrt(cos(tetav_jkp)*cos(tetah_jkp));
        fiv_jkp=sin(tetav_jkp);
        fih_jkp=sin(tetah_jkp)*cos(tetav_jkp);
        for m=1:sqrt(M)
            dv(m,1)=exp((m-1)*1i*(2*pi*deltaL/lambda)*fiv_jkp); 
            dh(m,1)=exp((m-1)*1i*(2*pi*deltaL/lambda)*fih_jkp); 
        end
        d=(1/sqrt(M))*(kron(dv,dh)); 
        for m=1:M
            D_jk(m,p)=alfaNL_jkp*d(m,1);
            R_jk(m,p)=LNL_jk(m,m)*D_jk(m,p);  % (R_jk)^1/2 : the deterministic correlation matrix from device j to LIS unit k
        end
        s281=norm(h_kk'*R_jk(:,p))^2;
        s28=s28+s281;
        s291=0;
        for m=1:M
            s2911=(alfaNL_jkp*betaL_kkm(m,1)*lNL_jkm(m,1))^2/M;
            s291=s291+s2911;
        end
        s29=s29+s291;
    end
    S28(j)=((1-T_k^2)/(ric(1,j)+1))*s28; 
    S29(j)=(T_k^2/(ric(1,j)+1))*s29; 
    H=R_jk*g_jk;
    for m=1:M
        hNL_jk(m,j)=H(m,1); % The correlated NLOS component of interference channel between device j and LIS unit k
        s751(m,j)=(norm(R_jk(m,:)))^2;
    end
end
toc
% The interference channel h_jk between device j and LIS unit k:combination of LOS and NLOS
h_jk=zeros(M,K-1);
d_c=10;
for j=1:K-1
    % probability of LOS path
    if (0<d_jk(1,j))&&(d_jk(1,j)<10)
        PLOS_jk=(d_c-d_jk(1,j))/d_c; 
    else
        PLOS_jk=0;        
    end
    a=rand;
    for m=1:M
        %h_jk(m,j)=PLOS_jk*sqrt(ric(1,j)/(ric(1,j)+1))* hL_jk(m,j)+(1-PLOS_jk)*sqrt(1/(ric(1,j)+1))* hNL_jk(m,j);
        if a<PLOS_jk
            h_jk(m,j)=sqrt(ric(1,j)/(ric(1,j)+1))* hL_jk(m,j);
        else
            h_jk(m,j)=sqrt(1/(ric(1,j)+1))* hNL_jk(m,j);
        end
    end
end

% pk : (17)
pk=atan(L^2/(sqrt(2*L^2+1)));
pk_=pk^2/(pi^2*deltaL^4);
% qk : (32)
qk=L^2/((L^2+1)*(2*L^2+1))+((L*(2*L^2+3))/((1*(L^2+1))^(3/2)))*atan(L/sqrt(L^2+1));
qk_=qk/(16*pi^2*deltaL^2);
% mean and variance of Yjk 
muL_jk=zeros(M,K-1);
muY_jk=zeros(1,K-1);
varY_jk=zeros(1,K-1);
A=0;
C=0;
for j=1:K-1 % k = 1: K except target device, td
    s26=0;
    s27=0;
    for m=1:M
        s261= h_kk(m,1)'*hL_jk(m,j); 
        s26=s26+s261;
        s271=(betaL_kkm(m,1)*betaL_jm(m,j))^2; 
        s27=s27+s271;
    end
    S26=(sqrt(ric(1,j)*(1-T_k^2)/(ric(1,j)+1)))*s26;
    S27=(ric(1,j)*(T_k^2)/(ric(1,j)+1))*s27;
    muY_jk(j)=S27+S28(j)+S29(j)+norm(S26)^2;
    varY_jk(j)=(S27+S28(j)+S29(j))^2+2*norm(S26)^2*(S27+S28(j)+S29(j)); 
    
    A1=p_j*muY_jk(j);
    A=A+A1;  
    C1=p_j^2*varY_jk(j);
    C=C+C1; 
end
% mean and variance of Zk 
for m=1:M
    varzwL_km(m)=delta*p_k*betaL_kkm(m,1)^4; 
end
sum30=0;
for m=1:M
    su72=0;
    sm74=0;
    s75=0;
    su76=0;
    for j=1:K-1 % k = 1: K except td=(K-1)/2+1 target device
        su721=(sqrt(p_j*ric(1,j)/(ric(1,j)+1))*betaL_kkm(m,1)*betaL_jm(m,j)*conj(h_kk(m,1))*h_jk(m,j));
        su72=su72+su721;    
        sm741=(sqrt(p_j*ric(1,j)/(ric(1,j)+1))*betaL_jm(m,j)*h_jk(m,j));
        sm74=sm74+sm741;  
        s75=s75+(p_j/(ric(1,j)+1))*s751(m,j); 
        su761=sqrt(p_j*ric(1,j)/(ric(1,j)+1))*betaL_jm(m,j)*h_kkm(m,1)*conj(h_jk(m,j));
        su76=su76+su761; 
    end
    firstp(m)=(1-T_k^2)*norm(su72)^2; 
    secp(m)=betaL_kkm(m,1)^2*(s75+T_k^2*norm(sm74)^2);
    varzwR_km(m)=delta*(firstp(m)+secp(m)); 
    omegawLR_km(m)=delta*sqrt(p_k)*betaL_kkm(m,1)^3*su76; 
    sum301=varzwL_km(m)+varzwR_km(m)+2*real(omegawLR_km(m));
    sum30=sum30+sum301;
end
muzw_k=sum30; 
varzw_k=sum30^2; %
% mean and variance of Ik
muc_tk=zeros(1,K-1);
mua_tmk=zeros(M,K-1);
for j=1:K-1 % k = 1: K except td=(K-1)/2+1 target device
    tk=0;
    for m=1:M
        tk1=h_kk(m,1)'*hL_jk(m,j);
        tk= tk+tk1;
        mua_tmk(m,j)=sqrt(T_k^2*ric(1,j)/(ric(1,j)+1))*betaL_kkm(m,1)*betaL_jm(m,j)*h_jk(m,j); 
   end
   muc_tk(j)=sqrt(ric(1,j)*(1-T_k^2)/(ric(1,j)+1))*tk;
end
f=0;
for i=1:K-1
     for j=1:K-1
         if i~=j
            W=0;
            for m=1:M
                 w1=conj(mua_tmk(m,i))* mua_tmk(m,j);
                 W=W+w1;
            end
            w_ijk(i,j)=2*real(muc_tk(i)*conj(muc_tk(j))*W); 
            f1=p_j*p_j*w_ijk(i,j);
            f=f+f1;
         end
     end
 end  
muI_k=p_k*T_k^2*qk_+A+sqrt(pk_)+muzw_k;
varI_k=p_k^2*T_k^4*qk_^2+T_k^2*(2-T_k^2)*qk_+C+f+varzw_k; 
% mean and variance of the asymptotic uplink SINR 
mug_k=p_k*pk_*(1-T_k^2)*((1/muI_k)+(varI_k/(muI_k^3))); 
varg_k=p_k^2*pk_^2*(1-T_k^2)^2*((varI_k/muI_k^4)-(varI_k^2/muI_k^6));
% asymptotic mean and variance of Rk 
muR_k=log2(1+mug_k)-varg_k/(2*(1+mug_k)^2); %  estimate the ergodic rate 
varR_k=varg_k/(1+mug_k)^2-varg_k^2/(4*(1+mug_k)^4); %  verify system reliability

R_k(t)= muR_k;
variance(t)= varR_k;
end
% averaging
ergodicR_k(s)=mean(R_k);
varianc(s)=mean(variance);
end
% plots
figure(1)
plot(w,ergodicR_k,'b-O','DisplayName','K=20, MRC');
%title('Ergodic Capacity of an LIS-based system with randomly located devices');
xlabel('Number of antennas on each LIS unit (M)');
ylabel('Ergodic Capacity (bps/Hz)');
grid on;
legend show;
ylim([2 6.5]);

    
    