clc; clear all;
syms phi delta phi_p delta_p T_phi T_delta

M = [4.09537, 0.31541; 0.31541, 0.16216]
C1 = [0, 4.11551; -0.58511, 0.51821]
K0 = [-61.76730 , -4.80033 ; -4.80033, -1.31521]
K2 = [0, 6.69605 ;0, 0.57886]
v = 3.6 %4.012 %3.19
M_inv = inv(M)
Vector_Estados = [phi ; delta; phi_p; delta_p]
Vector_Entrada = [T_phi; T_delta]

A = [zeros(2),eye(2);-M_inv * (K0 + (v^2) * K2),-M_inv * (v * C1)]
B = [zeros(2);M_inv*[1;0],M_inv*[0;1]]
C = eye(4)
D = [0, 0; 0 ,0; 0 ,0; 0 ,0]


tmax=10;
Aa=[A zeros(4,1);-[0,0,0,1] 0]
Ba=[B;[0,0]]
Ca = eye(5)
Da = [0, 0; 0 ,0; 0 ,0; 0 ,0; 0,0]
% 
Qa=[6.7 0 0 0 0; 0 0.3 0 0 0 ; 0 0 0.2 0 0;0 0 0 0.2 0;0 0 0 0 0.1];
R=0.4;
% 
 sys=ss(Aa,Ba,Ca,Da);
 
 Kalqr=lqr(sys,Qa,R)
 ke=Kalqr(:,1:4);
 ki=Kalqr(:,2:5);
 eabkilqr=eigs(Aa-Ba*Kalqr)
 x0=[0.087 ,0.34 ,0.5, 0.5];
 [var]=sim('rctalineal_integrador_bicicleta',[0 tmax])
 x=var.x;
 u=var.u;
 t=var.t;
% 
 figure
% 
 subplot(5,1,1)
 plot(t,x(:,1),'b')
 hold on
 line([0 10],[0 0],'Color','red')
 ylabel('Angulo Inclinación')
% 
 subplot(5,1,2)
 plot(t,x(:,2),'b')
 hold on
 line([0 10],[0 0],'Color','red')
 ylabel('Angulo Dirección')
% 
 subplot(5,1,3)
 plot(t,x(:,3),'b')
 hold on
 line([0 10],[0 0],'Color','red')
 ylabel('Vel Inclinación')
% 
subplot(5,1,4)
plot(t,x(:,3),'b')
hold on
line([0 10],[0 0],'Color','red')
ylabel('Vel Direccion')
subplot(5,1,5)
% 
plot(t,u)
ylabel('Entradas')
xlabel('t')
% 
figure
plot(x(:,1),x(:,2))
ylabel('Angulo Dirección')
xlabel('Angulo inclinación')


Aa=[A , zeros(4,1);-[0 1 0 0] ,0]
eig(Aa)
Pa = [-4 ,-4 ,-6 ,-3 ,-1]
%[-4.23 ,-4.23 ,-6.69 ,-3.6 ,-1.3]
Ka=place(Aa,Ba,Pa)
eabki=eigs(Aa-Ba*Ka)
ke=Ka(:,1:4);
ki=Ka(:,5);
x0=[0.087 ,0.34 ,0.5, 0.5];
[var]=sim('rctalineal_integrador_bicicleta',[0 tmax])
x=var.x;
u=var.u;
t=var.t;
figure

subplot(5,1,1)
plot(t,x(:,1),'b')
hold on
line([0 10],[0 0],'Color','red')
ylabel('Angulo Inclinación')

subplot(5,1,2)
plot(t,x(:,2),'b')
hold on
line([0 10],[0 0],'Color','red')
ylabel('Angulo Dirección')

subplot(5,1,3)
plot(t,x(:,3),'b')
hold on
line([0 10],[0 0],'Color','red')
ylabel('Vel Inclinación')

subplot(5,1,4)
plot(t,x(:,3),'b')
hold on
line([0 10],[0 0],'Color','red')
ylabel('Vel Direccion')
subplot(5,1,5)

plot(t,u)
ylabel('Entradas')
xlabel('t')

figure
plot(x(:,1),x(:,2))
ylabel('Angulo Dirección')
xlabel('Angulo inclinación')

x=var.x;
% x1(:,1)=x(1,1,:);
% x2(:,1)=x(2,1,:);
% x = [x1,x2];
u=var.u;
t=var.t;
Cep=var.Cep;
Tep=var.Tep;

figure('name','CONTROL ROBUSTO (integrador) Col Polos vs Control óptimo')
subplot(2,2,1)
plot(t,x(:,1),'b')
hold on
line([0 tmax],[xeq(1,1) xeq(1,1)],'Color','red')
hold on
line([0 tmax],[xeq(1,2) xeq(1,2)],'Color','green')
hold on
line([0 tmax],[xeq(1,3) xeq(1,3)],'Color','cyan')
hold on
plot(t,Cep,'k')
ylim([0 1.2])
ylabel('x1')
grid on

subplot(2,2,2)
plot(t,x(:,2),'b')
hold on
line([0 tmax],[xeq(2,1) xeq(2,1)],'Color','red')
hold on
line([0 tmax],[xeq(2,2) xeq(2,2)],'Color','green')
hold on
line([0 tmax],[xeq(2,3) xeq(2,3)],'Color','cyan')
hold on
plot(t,Tep,'k')
ylim([1.5 2.4])
ylabel('x2')
grid on

subplot(2,2,3)
plot(t,u,'b')
hold on
line([0 tmax],[Tc Tc],'Color','k')
ylabel('u')
xlabel('t')
grid on

subplot(2,2,4)
plot(x(:,1),x(:,2),'b')
xlabel('C')
ylabel('T')
ylim([1.75 2.05])
xlim([0.4 1])
grid on

%---------- CONTROL ÓPTIMO

Ca=eye(3,3);
Da=zeros(3,1);

% Qa=eye(3);
% R=1;

% Qa=[0.3 0 0 ; 0 5 0 ; 0 0 5];
% R=0.75;

% Qa=[0.3 0 0 ; 0 5 0 ; 0 0 5];
% R=0.01;

 Qa=[0.3 0 0 ; 0 5 0 ; 0 0 5];
 R=0.03;

sys=ss(Aa,Ba,Ca,Da);

Kalqr=lqr(sys,Qa,R)
ke=Kalqr(1:2);
ki=Kalqr(3);
eabkilqr=eigs(Aa-Ba*Kalqr)

[var]=sim('rctanolineal_integrador',[0 tmax])
x=var.x;
% x1(:,1)=x(1,1,:);
% x2(:,1)=x(2,1,:);
% x = [x1,x2];
u=var.u;
t=var.t;
Cep=var.Cep;
Tep=var.Tep;

subplot(2,2,1)
hold on
plot(t,x(:,1),'m')
hold on
plot(t,Cep,'--k')
ylim([0 1.2])
xlabel('t')
ylabel('x1')

legend('Col Pol Int','PE1','PE2','PE3','Ce,Te','Óptimo Int')

subplot(2,2,2)
hold on
plot(t,x(:,2),'m')

hold on
plot(t,Tep,'k--')
ylim([1.5 2.4])
xlabel('t')
ylabel('x2')

subplot(2,2,3)
hold on
plot(t,u,'m')
ylabel('u')
xlabel('t')
 
subplot(2,2,4)
hold on
plot(x(:,1),x(:,2),'m')
xlabel('C')
ylabel('T')
