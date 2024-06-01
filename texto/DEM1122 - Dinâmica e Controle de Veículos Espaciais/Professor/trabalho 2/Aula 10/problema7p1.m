function problema7p1

% Dinamica e Controle de Veiculos Espaciais

% Exemplo 7.1 da referencia Wie, B. Space Vehicle Dynamics and Control. 2 . ed., AIAA
% Education Series, Reston, VA: AIAA, 2008.

% Adotar a sequencia 123 de angulos de Euler.

clc;close all;clear all;
%% Entrada de parametros
global Ixx Iyy Izz MX MY MZ T
% Momentos de inercia principais do VE quase axis simetrico
Ixx=4223; % kg.m^2
Iyy=4133; % kg.m^2
Izz = 768; % kg.m^2
% Momento da manobra de spin up
MX=0;MY=0;MZ=10; %Nm
% Tempo da manobra de spin up
T=120; %s
% Condicoes iniciais
wx0=0.0001; % rad/s
wy0=0; % rac/s
wz0=0; % rad/s
%% Simulacao
% Condicao inicial
theta2=0;theta3=0;
% IMPORTANTE: supondo o vetor quantidade de movimento angular alinhado com
% o eixo n3 do sistema inercial
% Inserir um valor pequeno para calculo de theta1 em t=0, caso contrario, tem singularidade
% quando se calcula o angulo de nutacao
theta1=0.1*pi/180;

x0=[wx0;wy0;wz0;theta1;theta2;theta3];
% Integracao da dinamica
% As tolerâncias são importantes, pois pequenos desvios prejudicam os resultados
opt=odeset('RelTol',1e-12,'AbsTol',1e-12,'MaxStep',1);
[t,x]=ode45(@dinamica_torque_VE_rigido,[0 1.5*T],x0,opt);
%% Calcula os angulos de nutacao theta, spin relativo psi e precessao phi

phi=x(:,4)+pi/2; % phi=theta1+pi/2

theta=x(:,5); % theta=theta2

psi=x(:,6)-pi/2; % psi=theta3-pi/2

%% Calcula a orientacao dos vetores quantidade de movimento angular e velocidade angular

N=length(t);

phih=zeros(N,1);thetah=phih;psih=phih;

phiw=zeros(N,1);thetaw=phiw;psiw=phiw;

for i=1:N
    C=angle2dcm(x(i,4),x(i,5),x(i,6),'XYZ'); % Matriz de atitude, transforma do fixo para o corpo
    wb=[x(i,1);x(i,2);x(i,3)]; % Velocidade angular no referencial do corpo
    Hb=[Ixx*x(i,1);Iyy*x(i,2);Izz*x(i,3)]; % Quantidade de movimento angular no referencial do corpo
    w0=transpose(C)*wb; % Transformacao do referencial do corpo para o fixo
    H0=transpose(C)*Hb; % Transformacao do referencial do corpo para o fixo
    psiw(i)=atan2(w0(2),w0(1)); % Angulo de guinada da velocidade angular
    thetaw(i)=atan2(w0(3),sqrt(w0(1)^2+w0(2)^2)); % Angulo de elevacao da velocidade angular
    psih(i)=atan2(H0(2),H0(1)); % Anguio de guinada da quantidade de movimento anguiar
    thetah(i)=atan2(H0(3),sqrt(H0(1)^2+H0(2)^2)); % Angulo de elevacao da quantidade de movimento angular
end
%% Estudo analítico
disp('Velocidade de spin final pelo solucao analítica - rad/s')
disp(MZ*T/Izz);
disp('Velocidade de spin final pela simulacao nao linear - rad/s')
disp(x(end,3));
%% Graficos
figure
subplot(231);plot(t,x(:,1));grid;xlabel('t (s)');ylabel('\omega_x (rad/s)');axis tight
subplot(232);plot(t,x(:,2));grid;xlabel('t (s)');ylabel('\omega_y (rad/s)');axis tight
subplot(233);plot(t,x(:,3));grid;xlabel('t (s)');ylabel('\omega_z (rad/s)');axis tight
subplot(234);plot(t,x(:,4)*180/pi);grid;xlabel('t (s)');ylabel('\theta_1 (º)');axis tight
subplot(235);plot(t,x(:,5)*180/pi);grid;xlabel('t (s)');ylabel('\theta_2 (º)');axis tight
subplot(236);plot(t,x(:,6)*180/pi);grid;xlabel( 't (s)');ylabel('\theta_3 (º)');axis tight
figure
subplot(221);plot(x(:,4)*180/pi,x(:,5)*180/pi);grid;xlabel('\theta 1 (º)');ylabel('\theta_2 (º)');axis tight
subplot(222);plot(t,phi*180/pi);grid;xlabel('t (s)');ylabel('precessao \phi (º)');axis tight
subplot(223);plot(t,theta*180/pi);grid;xlabel('t (s)');ylabel('nutacao \theta (º)');axis tight
subplot(224);plot(t,psi*180/pi);grid;xlabel('t (s)');ylabel('spin relativo \psi (º)');axis tight
%% Animacao 3D
%Animação desativada
%{ 
% Momentos de inercia para calcular a animacao da quantidade de movimento
% angular
I=[Ixx Iyy Izz];
% Insere os angulos na sequencia 313 para gerar a animacao
x(:,4:6)=[phi theta psi];
% Escala de tempo
et=2;
% Chama a animacao
animacao_3d(t,x,l,et);
%} 

end


