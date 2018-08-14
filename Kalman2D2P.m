clc
clear all
close all


% here is the code for simple kalman filter

Ts=0.1; 
F_x=[1 0;0 1]; 
F_u=[Ts 0;0 Ts]; 
H=[1 0; 0 1]; 

segmaVx=.5;
segmaVy=.5; 

Q=[segmaVx.^2 0  ; 0 segmaVy.^2];

t=0:Ts:10;
sys=ss(F_x,F_u, eye(2),[],Ts);

vx=[zeros(1,30) .25*ones(1,20) -.20*ones(1,20) .15*ones(1,length(t)-70)]+normrnd(0,segmaVx,1,length(t));
vy=[zeros(1,10) .60*ones(1,60) -.20*ones(1,length(t)-70)]+normrnd(0,segmaVy,1,length(t));
v=[vx;vy];

% generating the true data for linear system
x0 = [0 0]';
Xtrue=lsim(sys,v,t,x0);
xtrue=Xtrue(:,1);
ytrue=Xtrue(:,2);

% generating measurement data by adding noise to the true data: 
segmaY1 = 0.1;
xm1=xtrue+normrnd(0,segmaY1,length(xtrue),1);
ym1=ytrue+normrnd(0,segmaY1,length(ytrue),1);
R1=[segmaY1^2 0;0 segmaY1^2];

segmaY2 = 0.2;
xm2=xtrue+normrnd(0,segmaY2,length(xtrue),1);
ym2=ytrue+normrnd(0,segmaY2,length(ytrue),1);
R2=[segmaY2^2 0;0 segmaY2^2];

% here we are initializing estimation
Xkfest1=zeros(2,length(t));
Pkf11=zeros(2,length(t));
KK1=zeros(2,length(t));

Xkfest1(:,1)=x0;
Pkf1=F_u*Q*F_u';
Pkf11(:,1)= diag(Pkf1);

Xkfest2=zeros(2,length(t));
Pkf22=zeros(2,length(t));
KK2=zeros(2,length(t));

Xkfest2(:,1)=x0;
Pkf2=F_u*Q*F_u';
Pkf11(:,1)= diag(Pkf2);

for i=2:length(t)
    Pkf1=F_x*Pkf1*F_x'+F_u*Q*F_u'; %predicting P
    Xkfest1(:,i) = F_x*Xkfest1(:,i-1)+F_u*v(:,i-1); %Predicitng the state
    K1=Pkf1*H'/(H*Pkf1*H'+R1); % Kalman gains  
    Xkfest1(:,i)=Xkfest1(:,i)+K1*([xm1(i); ym1(i)]-H*Xkfest1(:,i)); %Correcting: estimating the state
    Pkf1=(eye(2)-K1*H)*Pkf1; %Correcting: estimating P
    Pkf11(:,i) = diag(Pkf1);
    
    
    Pkf2=F_x*Pkf2*F_x'+F_u*Q*F_u'; %predicting P
    Xkfest2(:,i) = F_x*Xkfest2(:,i-1)+F_u*v(:,i-1); %Predicitng the state
    K2=Pkf2*H'/(H*Pkf2*H'+R2); % Kalman gains  
    Xkfest2(:,i)=Xkfest2(:,i)+K1*([xm2(i); ym2(i)]-H*Xkfest2(:,i)); %Correcting: estimating the state
    Pkf2=(eye(2)-K2*H)*Pkf2; %Correcting: estimating P 
    Pkf22(:,i) = diag(Pkf2);
end

figure
plot(t,Xkfest1(1,:),'r',t,Xkfest2(1,:),'--',t,xtrue,'g')
legend('xEstimated1','xEstimated2','xTrue')
xlabel('time [sec]');
ylabel('displacementX [m/s]');
title('displacementX');
figure
plot(t,Xkfest1(2,:),'r',t,Xkfest2(2,:),'--',t,ytrue,'g')
legend('yEstimated','yEstimated2','yTrue')
xlabel('time [sec]');
ylabel('displacementY [m/s]');
title('displacementY');



