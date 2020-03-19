%2020.03.18 
% compare the prediction accuracy of linearization with Taylor series and
%Koopman operator
clear,clc,
%
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
% %

fs=10000;
Ts=1/fs;
t=(0:Ts:5)';
amp=1.322;
testU=amp*(sin(2*pi*242.6*t+1.5*pi)+1);

Lseq=length(testU);
ekall=zeros(Lseq,1);
yPredTest=zeros(Lseq,1);%modeling error
% plot(testU);hold on,plot(testY);
Nrnn=20;

%
simN=1500;
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
Nh=8;%control horizon
Nc=8;
refk=zeros(Nh,1);

rpc=0.8;
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
%
cnt=0;
ypre_TL=zeros(simN,1);
flag=1;
xkTL=xk;
uk=uk/scale;
%generate other m-1 observalbes
mobs=20;%m observales
%measurement: Wm(i,:)*tanh(W1*xk+Bm(:,i))
Wm=(rand(mobs-1,size(W2,2))-0.5)*2;
Bm=(rand(size(B1,1),mobs-1)-0.5)*2;
B3m=(rand(mobs,1)-0.5)*5;
%Ks number of samples in each epoch
Ks=20;  %20
Kep=60; %number of epoch  20 when range is large make Kep large 
Kep1=0;
%
[b,a]=butter(6,0.15);
ukSamx=(rand(Kep*Ks,1)-0.5)*4;
ukSamy=filtfilt(b,a,ukSamx);
% plot(y),hold on, %plot(x);
% ukSam=(reshape(ukSamy,Ks,Kep))';
ukSam=(rand(Kep-Kep1,Ks)-0.5)*2;
for k=1:Kep
    tmp=ukSam(k,1);
    ukSam(k,:)=ukSam(k,:)-tmp;
end

Np=12;
%
Xtmp=zeros(mobs,Ks*Kep); %used for constructing Xlift and Xlift_p
ypre_Koop=zeros(simN,1);

Xlift=zeros(mobs+1,(Ks-1)*Kep);
Xlift_p=zeros(mobs,(Ks-1)*Kep);

y_t1=zeros(simN,1);
y_t2=zeros(simN,1);

Wg=W1;
Wslope=(rand(1,mobs-1)-0.5)*0+1;

flag=1;
for k=1:simN-2      
    %Taylor linearization
    %xk=Apc*xk+Bpc*uk+B_c  : B_c constant term
    %yk=W2*xk+B3
    if(flag==1)
        %update coefficients
        f0=W1*xk+B2+B1*uk(k);
        tanhf0=tanh(f0);
        df0=1-tanhf0.^2;
        Apc=W1.*repmat(df0,1,size(W1,2));
        Bpc=df0.*B1;
        B_c=tanhf0-Apc*xk-Bpc*uk(k);  
        flag=0;
        xkTL=xk;
        
        %Koopman operator       
        for k_koop=1:Kep %Kep epochs
            xkKoop=xk;
            for k_koop1=1:Ks
                Xtmp(1,k_koop1+Ks*(k_koop-1))=(W2*xkKoop+B3)*scale;
                tmp=tanh( repmat(Wslope,20,1).*repmat(Wg*xkKoop,1,mobs-1)+Bm);
%                 tmp=( repmat(Wg*xkKoop,1,mobs-1)+Bm);
                for k_koop2=2:mobs
                    Xtmp(k_koop2,k_koop1+Ks*(k_koop-1))=(Wm(k_koop2-1,:)*tmp(:,k_koop2-1)+B3m(k_koop2))*scale;
                end          
                %Radp: range adapter
                rng=max(uk(k:k+Np))-min(uk(k:k+Np));
                Radp=rng/2;
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
        tmp=tanh( repmat(Wslope,20,1).*repmat(Wg*xk,1,mobs-1)+Bm );
%         tmp=( repmat(Wg*xk,1,mobs-1)+Bm );
        for k_koop2=2:mobs
           xkKoop(k_koop2)=(Wm(k_koop2-1,:)*tmp(:,k_koop2-1)+B3m(k_koop2))*scale;
        end    
        
        y_t1(k)=CpcK*xkKoop; 
        
        xkKooptmp=ApcK*xkKoop+BpcK*uk(k);
         
        y_t2(k+1)=CpcK*xkKooptmp; 
        cnt=0;
    else
        if(cnt==Np)  %prediction horizon
            flag=1;
%             cnt=0;
        end
    end    
    %linearized system dynamics: Taylor series
    ypre_TL(k)=(W2*xkTL+B3)*scale;
    xkTL=Apc*xkTL+Bpc*uk(k)+B_c;
    
    cnt=cnt+1;
    
    %linearization based on Koopman operator  
    ypre_Koop(k)=CpcK*xkKoop; 
    xkKoop=ApcK*xkKoop+BpcK*uk(k);
       
    
    %system dynamics
    yk(k)=(W2*xk+B3)*scale;
    xk=tanh(W1*xk+B2+B1*uk(k));
    
    
   
end
%
%prediction results
ind1=100;
% subplot(1,3,1),
% plot(yk(ind1:simN)),hold on,plot(ypre_TL(ind1:simN));
% subplot(1,3,2),
% plot(yk(ind1:simN)-ypre_TL(ind1:simN));
% 
% subplot(1,3,3),
subplot(1,2,1),
plot(yk(ind1:simN),'red','linewidth',2),
hold on,plot(ypre_Koop(ind1:simN),'green--','linewidth',2); 
plot(ypre_TL(ind1:simN),'blue--','linewidth',2)% plot(y_t1(ind1:simN),'*'); plot(y_t2(ind1:simN),'x');
legend('ref','Koopman','Taylor');
subplot(1,2,2),
plot(yk(ind1:simN)-ypre_Koop(ind1:simN),'red--','linewidth',2),hold on,
plot(yk(ind1:simN)-ypre_TL(ind1:simN),'blue--','linewidth',2);
legend('Koopman','Taylor');

e1=yk(ind1:simN)-ypre_Koop(ind1:simN);

e2=yk(ind1:simN)-ypre_TL(ind1:simN);

disp([norm(e1) norm(e2)]);
disp([norm(e1,inf) norm(e2,inf)]);
disp([var(e1) var(e2)]);
disp('finished!');












