%2020.03.18
%use two linear models to approxmiate the nonlinear model
clear,clc,
load('E:\back up flash\20191230\20170814NNIdentifymodel\lib\20190605_RNN_linearization\github\net12.mat')
load('E:\back up flash\20191230\20170814NNIdentifymodel\lib\20190605_RNN_linearization\github\xstar_np.mat')
global W1
global B1
global B2
global W2
global B3
global xkpc
global ukpc
global refk %col vector
global Nh
global Nc
global Abar
global Bbar
global Ae
global Be 
global Gall
global Call
global Cbar
global scale

scale=10.0;

%extract rnn from matlab
IW=net.IW;
LW=net.LW;
rnnb=net.b;%all the biaes values

%nonlinear state space equation
%    xk=tanh( W1*xk+B2+B1*uk(k) );
%    y(k)=W2*xk+B3;
B1=cell2mat(IW(1));
W1=cell2mat(LW(1));
B2=cell2mat(rnnb(1)); 
W2=cell2mat(LW(2));
B3=cell2mat(rnnb(2));

Ns=2;
Ae=reshape(xstar(1:Ns^2),Ns,Ns);
Be=xstar(Ns^2+1:Ns^2+Ns,1);
Ge=xstar(Ns^2+Ns+1:Ns^2+Ns*2,1);
Ce=(xstar(Ns^2+Ns*2+1:Ns^2+Ns*3,1))';
xstarpara=xstar;
%
% testU=(y1range(2)-y1range(1))/(x1range(2)-x1range(1))*(Xtr-x1range(1))+y1range(1);
% testU=Xtr*scale;
% testY=Ttr*scale;% to be predicted
fs=10000;
Ts=1/fs;
t=(0:Ts:5)';
amp=1.2;
testU=amp*(sin(2*pi*250*t+1.5*pi)+1);

Lseq=length(testU);
ekall=zeros(Lseq,1);
yPredTest=zeros(Lseq,1);%modeling error
% plot(testU);hold on,plot(testY);
Nrnn=20;

simN=3000;
%LPF state space equation
Nf=3;%order of filter
[b,a]=butter(Nf,0.2);
[Abar,Bbar,Cbar,Dbar]=tf2ss(b,a);
bxn=zeros(Nf+Nrnn+Ns,1);%state [beta(k); x(k); eta(k)]

ek=0;

Gall=[zeros(Nrnn+Nf,1);Ge];
Call=[zeros(1,Nrnn+Nf) Ce];
refsig=testU;%reference signal
yk=zeros(simN,1);
uk=testU;
betak=zeros(Nf,1);
xk=zeros(Nrnn,1);
etak=zeros(Ns,1);
%parameter for UKF
Nall=Ns+Nrnn+Nf;
xkhat=zeros(Nall,1);%ukf estimated states
Rv=eye(Nall)*0.000002;%process noise
Rn=0.005;%measurement noise, scalar!!!
Pk=Rv;
%sigma points parameters alpha beta gamma
alpha=0.00001;
kappa=3/Nall/1000;%can be set to zero
beta=2;%for guassian 2 is optimal
lambda=alpha^2*(Nall+kappa)-Nall;
gamma=sqrt(Nall+lambda);
W0m=lambda/(Nall+lambda);
W0c=lambda/(Nall+lambda)+(1-alpha^2+beta);
Wim=1/2/(Nall+lambda);
Wic=Wim;

y_est=zeros(simN,1);

WimMat=diag([W0m ones(1,2*Nall)*Wim]);
WicMat=diag([W0c ones(1,2*Nall)*Wic]);
Pk_=zeros(Nall);

%noise
measn=normrnd(0,0.01,simN,1);
prn=normrnd(0,0.0002,Nall,simN);

%predictive control parameter
xkpc=zeros(Nall,1);
ukpc=0;
Nh=14;%control horizon
Nc=8;
refk=zeros(Nh,1);

rpc=0.8;
%
options=optimoptions(@fmincon,'Display','off');

uoptGD=zeros(simN,1);
opterror=zeros(simN,2);
ek1=0;
%

%add a weight matrix We for the tracking error
dwe=1-(1:Nh)/Nh+0.1;
We=diag(dwe.^2);
We=eye(Nh);
for k=1:simN
    %system dynamics
    betak=bxn(1:Nf);
    xk=bxn(Nf+1:Nf+Nrnn);
    etak=bxn(Nf+Nrnn+1:end);
    
%     uk(k)=(y1range(2)-y1range(1))/(x1range(2)-x1range(1))*(uk(k)-x1range(1))+y1range(1);
    tmp=Cbar*betak/scale;
    bxn=[Abar*betak+Bbar*uk(k);
        tanh(W1*xk+B2+B1*tmp);
         Ae*etak+scale*Be*(W2*xk+B3)]+Gall*ek1;    
    yk(k)=Call*bxn+ek1;
    ek1=measn(k)*0;
    
    %state estimation
    %Sigma points
    sqrtPk=chol(Pk,'lower');%cholesky factorization
    xi=[xkhat repmat(xkhat,1,Nall)+gamma*sqrtPk repmat(xkhat,1,Nall)-gamma*sqrtPk];
    %Time update
    uk_1=uk(k);
    betaki=xi(1:Nf,:);
    xki=xi(Nf+1:Nf+Nrnn,:);
    etaki=xi(Nf+Nrnn+1:end,:);
    tmp=Cbar*betaki/scale;
    xikk_1=[Abar*betaki+repmat(Bbar*uk_1,1,2*Nall+1);
            tanh(W1*xki+repmat(B2,1,2*Nall+1)+B1*tmp);
            Ae*etaki+scale*Be*(W2*xki+B3)]+repmat(Gall*ek,1,2*Nall+1);
    xkhat_=sum(xikk_1*WimMat,2);
    Pk_=Rv;
    for k1=1:2*Nall+1
        Pk_=Pk_+WicMat(k1,k1)*( (xikk_1(:,k1)-xkhat_)*(xikk_1(:,k1)-xkhat_)' );
    end
    ykk_1=Call*xikk_1;%row vector
    ykhat_=sum(ykk_1*WimMat);%scalar
    %Measurement update    
    Pykayka=Rn;%yka means yk with ~ above
    for k1=1:2*Nall+1
        Pykayka=Pykayka+WicMat(k1,k1)*( (ykk_1(:,k1)-ykhat_)*(ykk_1(:,k1)-ykhat_)' );
    end
    Pxkyk=zeros(Nall,1);
    for k1=1:2*Nall+1
        Pxkyk=Pxkyk+WicMat(k1,k1)*( (xikk_1(:,k1)-xkhat_)*(ykk_1(:,k1)-ykhat_)' );
    end
    Kk=Pxkyk/Pykayka;%Pxkyk*inv(Pykayka)
    xkhat=xkhat_+Kk*(yk(k)-ykhat_);
    ek=yk(k)-ykhat_;
    Pk=Pk_-Kk*Pykayka*Kk';
    y_est(k)=Call*xkhat; 
    %predictive control: fmincon-------------------------------------------
    xkpc=xkhat;
%     xkpc=bxn;
    ukpc=uk(k);
    refk(1:Nc)=(refsig(k+1:k+Nc,1)); %col vector
    if(Nh>Nc)
        refk(Nc+1:Nh)=refk(Nc);
    end    
%     lb=-ones(4,1);
%     ub=ones(4,1);
    lb=[];ub=[];
    x0=ones(Nc,1)*uk(k);    
%     xstar=fmincon(@preError_1_V3,x0,[],[],[],[],lb,ub,[],options);%dont forget @  !!
%     uk(k+1)=xstar(1);   
   
%second method:linearization at each point  20190606
%current state is xkpc, current input is ukpc
%convert to M_k+1=Apc*M_k+Bpc*uk+Mc  where Mc is constant
%output:    y_k=Call*M_k
Npc=Nf+Nrnn+Ns;%order of the model
Apc=zeros(Npc,Npc);
Bpc=zeros(Npc,1);
Cpc=Call;
Apc(1:Nf,1:Nf)=Abar;
Apc(Nf+Nrnn+1:end,Nf+1:end)=[scale*Be*W2 Ae];
Bpc(1:Nf)=Bbar;
%compute f0 where the linearization happens
Mc=zeros(Npc,1);
betak_0=xkpc(1:Nf);
xk_0=xkpc(Nf+1:Nf+Nrnn);
f0=W1*xk_0+B2+B1*Cbar*betak_0/scale;
tanhf0=tanh(f0);
df0=1-tanhf0.^2;
tmpB1CW=[B1*Cbar/scale W1];
dfdbetax=repmat(df0,1,size(tmpB1CW,2)).*tmpB1CW;
Mc(Nf+1:Nf+Nrnn)=-dfdbetax*xkpc(1:Nf+Nrnn)+tanhf0;
Mc(Nf+Nrnn+1:end)=scale*Be*B3;
Apc(Nf+1:Nf+Nrnn,1:Nf+Nrnn)=dfdbetax;
%predictive control setup
Gpc=zeros(size(Cpc,1)*Nh,size(Apc,2));
Hpc=zeros(Nh,Nh);
Fpc=zeros(Nh,1);
hpc=zeros(Nh+1,1);
hpc(Nh+1)=Cpc*Apc^(Nh-1)*Bpc;
Lpc=zeros(Nh,1);

Dpc=eye(Nh);
Dpc(Nh+1:Nh+1:(Nh-1)*Nh+Nh-1)=-1;
Dpc(end)=0;
[V_c,D_c] = eig(Apc);
% V_c=real(V_c);
% D_c=real(D_c);
%MPC parameters setup
    for k2=1:Nh
        if(k2==1)
            hpc(k2)=0;
            Lpc(1)=Cpc*Mc;
        else
            hpc(k2)=Cpc*real( V_c*D_c^(k2-2)/V_c )*Bpc;  %Apc^(k2-2)
            Lpc(k2)=Lpc(k2-1)+Cpc*real(V_c*D_c^(k2-1)/V_c)*Mc;
        end
        Gpc( (k2-1)*size(Cpc,1)+1:k2*size(Cpc,1),: )=Cpc*real(V_c*D_c^k2/V_c);
        Hpc(k2,1:k2)=(hpc(k2:-1:1,1))';
    end

Fpc=hpc(2:end);
Qpc=rpc*(Dpc'*We*Dpc)+Hpc'*We*Hpc;
fpc=Hpc'*We*(Gpc*xkpc+Fpc*ukpc-refsig(k+1:k+Nh,1)+Lpc);
xopt=-Qpc\fpc;
uk(k+1)=xopt(1); 
%set up new linearization for better approximation
%compute linearization point
xkSL=xkpc;
lenSL=5;
ukSL=[ukpc; xopt(1:lenSL-1)];
for kSL=1:lenSL  %nonlinear forward
    % xkSL=Apc*xkSL+Bpc*ukSL(kSL)+Mc;
    betakSL=xkSL(1:Nf);
    xk_sl=xkSL(Nf+1:Nf+Nrnn);
    etakSL=xkSL(Nf+Nrnn+1:end);
    tmpSL=Cbar*betakSL/scale;
    xkSL=[Abar*betakSL+Bbar*ukSL(kSL);
        tanh(W1*xk_sl+B2+B1*tmpSL);
         Ae*etakSL+scale*Be*(W2*xk_sl+B3)];   
end
%linearize at xkSL
%new linear model
Apc2=Apc;
%compute f0 where the linearization happens
Mc2=Mc;
betak2_0=xkSL(1:Nf);
xk2_0=xkSL(Nf+1:Nf+Nrnn);
f20=W1*xk2_0+B2+B1*Cbar*betak2_0/scale;
tanhf20=tanh(f20);
df20=1-tanhf20.^2;
tmpB1CW=[B1*Cbar/scale W1];
dfdbetax2=repmat(df20,1,size(tmpB1CW,2)).*tmpB1CW;
Mc2(Nf+1:Nf+Nrnn)=-dfdbetax2*xkSL(1:Nf+Nrnn)+tanhf20;
Apc2(Nf+1:Nf+Nrnn,1:Nf+Nrnn)=dfdbetax2;
lenSL=7;
 for k2=lenSL+1:Nh
%         hpc(k2)=Cpc*real( V_c*D_c^(lenSL-2)/V_c )*Apc2^(k2-lenSL)*Bpc;  %Apc^(k2-2)
%         Lpc(k2)=Lpc(k2-1)+Cpc*real(V_c*D_c^(lenSL-1)/V_c)*Apc2^(k2-lenSL)*Mc2;
%        
%         Gpc( (k2-1)*size(Cpc,1)+1:k2*size(Cpc,1),: )=Cpc*real(V_c*D_c^lenSL/V_c)*Apc2^(k2-lenSL);
        
        hpc(k2)=Cpc*Apc2^(k2-lenSL)*real( V_c*D_c^(lenSL-2)/V_c )*Bpc;  %Apc^(k2-2)
        Lpc(k2)=Lpc(k2-1)+Cpc*Apc2^(k2-lenSL)*real(V_c*D_c^(lenSL-1)/V_c)*Mc2;
       
        Gpc( (k2-1)*size(Cpc,1)+1:k2*size(Cpc,1),: )=Cpc*Apc2^(k2-lenSL)*real(V_c*D_c^lenSL/V_c);
        Hpc(k2,1:k2)=(hpc(k2:-1:1,1))';
 end
hpc(Nh+1)=Cpc*Apc^(lenSL-2)*Apc2^(Nh+1-lenSL)*Bpc;
Fpc=hpc(2:end);
Qpc=rpc*(Dpc'*Dpc)+Hpc'*Hpc;
fpc=Hpc'*(Gpc*xkpc+Fpc*ukpc-refsig(k+1:k+Nh,1)+Lpc);
xopt=-Qpc\fpc;
uk(k+1)=xopt(1);   

end

%
%tracking results and tracking errors
subplot(1,2,1),
plot(yk(1:simN)),hold on,plot(refsig(1:simN));
subplot(1,2,2),
plot(yk(1:simN)-refsig(1:simN));













