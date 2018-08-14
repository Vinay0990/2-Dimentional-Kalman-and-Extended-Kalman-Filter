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
t=0:Ts:15;

segmavx=.5;
segmavy=.5; 
% here the values are from experiments

vx=[zeros(1,30) .25*ones(1,20) -.20*ones(1,20) .15*ones(1,length(t)-70)]+normrnd(0,segmavx,1,length(t));
vy=[zeros(1,10) .60*ones(1,60) -.20*ones(1,length(t)-70)]+normrnd(0,segmavy,1,length(t));
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
Xkfest=zeros(2,length(t));
Xkfest(:,1)=x0;

%% here are the noise
R=measurementsV*C*C';
Q=[segmavx.^2 0  ; 0 segmavy.^2];

%Initializing P

Pkf=B*Q*B';

for(i=2:1:length(t))
Pkf=A*Pkf*A'+B*Q*B'; %predicting P
Xkfest(:,i)=A*Xkfest(:,i-1)+B*v(:,i-1); %Predicitng the state
K=Pkf*C'/(C*Pkf*C'+R); % Kalman gains
Xkfest(:,i)=Xkfest(:,i)+K*([xm(i); ym(i)]-C*Xkfest(:,i)); %Correcting: estimating the state
Pkf=(eye(2)-K*C)*Pkf; %Correcting: estimating P
end

figure
plot(t,Xkfest(1,:),'r',t,xm,'--',t,xtrue,'g')
legend('xEstimated','xAttenuated','xTrue')
xlabel('time [sec]');
ylabel('displacementX [m/s]');
title('displacementX');
figure
plot(t,Xkfest(2,:),'r',t,ym,'--',t,ytrue,'g')
legend('xEstimated','xAttenuated','xTrue')
xlabel('time [sec]');
ylabel('displacementY [m/s]');
title('displacementY');
