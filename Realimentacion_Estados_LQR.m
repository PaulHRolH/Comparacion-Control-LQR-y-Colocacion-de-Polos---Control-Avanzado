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
Qa=[0.5 0 0 0 0; 0 0.4 0 0 0 ; 0 0 0.3 0 0;0 0 0 0.3 0;0 0 0 0 0.2];
R=2;
 
sys=ss(Aa,Ba,Ca,Da);
 
Kalqr=lqr(sys,Qa,R)
K=Kalqr(:,1:4);
eabkilqr=eigs(Aa-Ba*Kalqr)
x0=[0.087 ,0.34 ,0.5, 0.5];
[var]=sim('rctalineal_simple_bicicleta',[0 tmax])
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