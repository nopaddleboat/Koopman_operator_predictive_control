%2020.03.18
%linearize the nonlinear mode with Koopman operators and design predictive
%controller for output tracking.
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

Vnh=12;

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

fs=10000;
Ts=1/fs;
t=(0:Ts:5)';
amp=1.5;
testU=amp*(sin(2*pi*50*t+1.5*pi)+1);
% testU=rand(length(t),1)*amp;
%
Lseq=length(testU);
ekall=zeros(Lseq,1);
yPredTest=zeros(Lseq,1);%modeling error
% plot(testU);hold on,plot(testY);
Nrnn=20;

%
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
alpha=0.0001;
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
Nh=Vnh;%control horizon

Nc=8;
refk=zeros(Nh,1);

rpc=0.1;
%
options=optimoptions(@fmincon,'Display','off');

uoptGD=zeros(simN,1);
opterror=zeros(simN,2);
ek1=0;
%
%simN=900;
%add a weight matrix We for the tracking error
dwe=1-(1:Nh)/Nh+0.1;
We=diag(dwe.^2);
We=eye(Nh);

cnt=0;
ypre_TL=zeros(simN,1);
flag=1;
xkTL=xk;
uk=uk/scale;
%generate other m-1 observalbes
mobs=8;%m observales
%measurement: Wm(i,:)*tanh(W1*xk+Bm(:,i))
Wm=(rand(mobs-1,size(W2,2))-0.5)*2;
Bm=(rand(size(B1,1),mobs-1)-0.5)*2;
B3m=(rand(mobs,1)-0.5)*2;
%Ks number of samples in each epoch
Ks=Vnh;
Kep=50; %number of epoch
ukSam=(rand(Kep,Ks)-0.5)*2;
ukSam(:,1)=0;

Xtmp=zeros(mobs,Ks*Kep); %used for constructing Xlift and Xlift_p
ypre_Koop=zeros(simN,1);

y_t1=zeros(simN,1);
y_t2=zeros(simN,1);
%
xk=zeros(Nrnn,1);
% Ks=Vnh;
% Kep=20; %number of epoch
Xlift=zeros(mobs+1,(Ks-1)*Kep);
Xlift_p=zeros(mobs,(Ks-1)*Kep);
for k=1:simN-2      
    %Taylor linearization
    %xk=Apc*xk+Bpc*uk+B_c  : B_c constant term
    %yk=W2*xk+B3
       %Koopman operator       
        for k_koop=1:Kep %Kep epochs
            xkKoop=xk;
            for k_koop1=1:Ks
                Xtmp(1,k_koop1+Ks*(k_koop-1))=(W2*xkKoop+B3)*scale;
                tmp=tanh( repmat(W1*xkKoop,1,mobs-1)+Bm);
                for k_koop2=2:mobs
                    Xtmp(k_koop2,k_koop1+Ks*(k_koop-1))=(Wm(k_koop2-1,:)*tmp(:,k_koop2-1)+B3)*scale;
                end         
                rng=max(refsig(k:k+Nh+5))-min(refsig(k:k+Nh+5));
                Radp=rng/10/3;
                xkKoop=tanh(W1*xkKoop+B2+B1*(ukSam(k_koop,k_koop1)*Radp+uk(k)));
            end            
        end
        %construct Xlift and Xlift_p

        for k_koop=1:Kep
            Xlift(:,1+(Ks-1)*(k_koop-1):(Ks-1)*k_koop)=[Xtmp(:,1+Ks*(k_koop-1):-1+Ks*k_koop);ukSam(k_koop,1:Ks-1)*Radp+uk(k)];
            Xlift_p(:,1+(Ks-1)*(k_koop-1):(Ks-1)*k_koop)=Xtmp(:,2+Ks*(k_koop-1):Ks*k_koop);
        end
        
        AB=Xlift_p/Xlift; 
        %parameters for the lifted linear system
        ApcK=AB(:,1:mobs);
        BpcK=AB(:,end);
        CpcK=[1 zeros(1,mobs-1)]; 
%         disp(['loop:' num2str(k)]);
%         disp( num2str(AB*Xlift(:,1)-Xlift_p(:,1)));
        xkKoop=zeros(mobs,1);
        xkKoop(1)=(W2*xk+B3)*scale;
        tmp=tanh( repmat(W1*xk,1,mobs-1)+Bm );
        for k_koop2=2:mobs
           xkKoop(k_koop2)=(Wm(k_koop2-1,:)*tmp(:,k_koop2-1)+B3)*scale;
        end    
        
%         y_t1(k)=CpcK*xkKoop;         
%         xkKooptmp=ApcK*xkKoop+BpcK*uk(k);         
%         y_t2(k+1)=CpcK*xkKooptmp; %        
    %linearization based on Koopman operator  
    
    %predictive control  with
    %ypre_Koop(k)=CpcK*xkKoop; 
    %xkKoop=ApcK*xkKoop+BpcK*uk(k);
    Apc=ApcK;
    Bpc=BpcK;
    Cpc=CpcK;    
    
    xkpc=xkKoop;
    ukpc=uk(k);
    %predictive control setup    
    Gpc=zeros(size(Cpc,1)*Nh,size(Apc,2));
    Hpc=zeros(Nh,Nh);
    Fpc=zeros(Nh,1);
    hpc=zeros(Nh+1,1);
    hpc(Nh+1)=Cpc*Apc^(Nh-1)*Bpc;

    Dpc=eye(Nh);
    Dpc(Nh+1:Nh+1:(Nh-1)*Nh+Nh-1)=-1;
    Dpc(end)=0;

    %MPC parameters setup
        for k2=1:Nh
            if(k2==1)
                hpc(k2)=0;
            else
                hpc(k2)=Cpc*Apc^(k2-2)*Bpc;  %Apc^(k2-2)
            end
            Gpc( (k2-1)*size(Cpc,1)+1:k2*size(Cpc,1),: )=Cpc*Apc^k2;
            Hpc(k2,1:k2)=(hpc(k2:-1:1,1))';
        end

    Fpc=hpc(2:end);
    Qpc=rpc*(Dpc'*Dpc)+Hpc'*Hpc;
    fpc=Hpc'*(Gpc*xkpc+Fpc*ukpc-refsig(k+1:k+Nh,1));
    xopt=-Qpc\fpc;
    uk(k+1)=xopt(1);     
    
    %system dynamics
    yk(k)=(W2*xk+B3)*scale;
    xk=tanh(W1*xk+B2+B1*uk(k));   
end
disp('finished!');
%
%plot the tracking results and the tracking errors
subplot(1,2,1),
plot(yk(1:simN)),hold on,plot(refsig(1:simN));
subplot(1,2,2),
plot(yk(1:simN)-refsig(1:simN));











