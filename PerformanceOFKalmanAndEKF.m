clc
clear all
close all

%% here is the code for simple kalman filter

% Defining Parameter 
Ts=0.1; 
A=[1 0;0 1]; 
B=[Ts 0;0 Ts]; 
C=[1 0; 0 1]; 
x0=[0;0]; 
sys=ss(A,B, eye(2),[],Ts); 
t=0:Ts:3;
xsm=10*cos(2*pi*t);
ysm=10*sin(2*pi*t);

segmavx=.5;
segmavy=.5; 
% here the values are from experiments

vx=[xsm]+normrnd(0,segmavx,1,length(t));
vy=[ysm]+normrnd(0,segmavy,1,length(t));
v=[vx;vy];

%generating the true data
Xtrue=lsim(sys,v,t,x0);
xtrue=Xtrue(:,1);
ytrue=Xtrue(:,2);
measurementsV=[0.2.^2 0;0 0.2.^2];

% generating measurement data by adding noise to the true data:
xm=xtrue+normrnd(0,0.2,length(xtrue),1);
ym=ytrue+normrnd(0,0.2,length(ytrue),1);

% here we are initializing estimation
% here we are initializing estimation
Xkfest=zeros(2,length(t));
Xkfest(:,1)=x0;


Xekfest=zeros(2,length(t));
Xekfest(:,1)=x0;

%% here are the noise
R=measurementsV*C*C';
Q=[segmavx.^2 0  ; 0 segmavy.^2];

%Initializing P
Pkf=B*Q*B';
Pekf=B*Q*B';
Jk=[1 0; 0 1];
Vk=[1 0; 0 1];

Fk=[1 0;0 1]; % partial devt wrt states
Wk=[Ts 0;0 Ts]; %% partial devt wrt inputs


for i=2:1:length(t)
    
Pkf=A*Pkf*A'+B*Q*B'; %predicting P for KF
Xkfest(:,i)=A*Xkfest(:,i-1)+B*v(:,i-1); %Predicitng the state for KF
K=Pkf*C'/(C*Pkf*C'+R); % Kalman gains for KF
Xkfest(:,i)=Xkfest(:,i)+K*([xm(i); ym(i)]-C*Xkfest(:,i)); %Correcting: estimating the state for KF
Pkf=(eye(2)-K*C)*Pkf; %Correcting: estimating P  for KF 
 
    
    
    
Pekf=Fk*Pekf*Fk'+Wk*Q*Wk'; %predicting P for EKF
Xekfest(:,i)=A*Xekfest(:,i-1)+B*v(:,i-1); %Predicitng the state for EKF
K=Pekf*Jk'*inv( Jk*Pekf*Jk'+Vk*R*Vk') ; % Kalman gains for EKF
Xekfest(:,i)=Xekfest(:,i)+K*([xm(i); ym(i)]-C*Xekfest(:,i)); %Correcting: estimating the state  EKF
Pekf=(eye(2)-K*Jk)*Pekf; %Correcting: estimating P for EKF
end


figure
plot(t,xtrue,'b')
hold on
plot(t,Xkfest(1,:),'r')
plot(t,Xekfest(1,:),'m')
plot(t,xm,'g')
hold off;
legend('True', 'KF', 'EKF')
xlabel('time [sec]');
ylabel('displacementX [m/s]');
title('displacementX');

figure
plot(t,ytrue,'b')
hold on
plot(t,Xkfest(2,:),'r')
plot(t,Xekfest(2,:),'m')
plot(t,ym,'g')
hold off;
legend('True', 'KF', 'EKF')
xlabel('time [sec]');
ylabel('displacementY [m/s]');
title('displacementY');


figure
plot(xtrue,ytrue,'b')
hold on
plot(Xkfest(1,:),Xkfest(2,:),'r')
plot(Xekfest(1,:),Xekfest(2,:),'m')
plot(xm,ym,'g')
hold off;
xlabel('time [sec]');
ylabel('displacementY [m/s]');
title('displacementY');


