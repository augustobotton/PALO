function problema7p3
% Dinamica e Controle de Veiculos Espaciais
% Exemplo 7.3 da referencia Wie, B. Space Vehicle Dynamics and Control. 2 . ed., AAA
% Education Series, Reston, VA: AIAA, 2008.
% Adotar a sequencia 123 de angulos de Euler.
clc;close all;clear all;
%% Entrada de parametros
global Ixx Iyy Izz MX MY MZ T mu J
% Momentos de inercia principais do VE quase axis simetrico
Ixx=2000; % kg.m^2
Iyy=1500; % kg.m^2
Izz=1000; % kg.m^2
% Momento de inercia comum a todos os eixos da massa de propelente
J=18; % kg.m^2
% Coeficiente de amortecimento viscoso nas paredes do reservatorio de
% propelente
mu=30; % N.m.s
% Momento da manobra flap-spin (Movimento livre sem momentos)
MX=0;MY=0;MZ=0; % N.m
% Tempo da manobra de spin up
T=1000; %s
% Primeira condicao inicial
wx0a=0.1224; % rad/s
wy0a=0; % rad/s
wz0a=2.99; % rad/s
% Segunda condicao inicial
wx0b=0.125; % rad/s
wy0b=0; % rad/s
wz0b=2.99; % rad/s
%% Estuda a resposta para a primeira condicao inicial
% Condicao inicial
theta2=0;theta3=0;
% IMPORTANTE: supondo o vetor quantidade de movimento angular alinhado com
% o eixo n3 do sistema inercial
% Inserir um valor pequeno para calculo de theta1 em t=0, caso contrario, tem singulaidade —
% quando se calcula o angulo de nutacao
theta1=0.1*pi/180;
x0=[wx0a;wy0a;wz0a;theta1;theta2;theta3;0;0;0];
% Integracao da dinamica
% As tolerâncias são importantes, pois pequenos desvios prejudicam os resultados
opt=odeset('RelTol',1e-12,'AbsTol',1e-12,'MaxStep',1);
[t,x]=ode45(@dinamica_VE_semi_rigido,[0 T],x0,opt);
%% Calcula os angulos de nutacao theta, spin relativo psi e precessao phi
phi=x(:,4)+pi/2; % phi=theta1+pi/2
theta=x(:,5); % theta=theta2
psi=x(:,6)-pi/2; % psi=theta3-pi/2
%% Calcula a orientacao dos vetores quantidade de movimento angular e velocidade angular
N=length(t);
phih=zeros(N,1);thetah=phih;psih=phih;
phiw=zeros(N,1);thetaw=phiw;psiw=phiw;
for i=1:N
    C=angle2dcm(x(i,4),x(i,5),x(i,6), 'XYZ'); % Matriz de atitude, transforma do fixo para o corpo
    wb=[x(i,1);x(i,2);x(i,3)]; % Velocidade angular no referencial do corpo
    Hb=[Ixx*x(i,1);Iyy*x(i,2);Izz*x(i,3)]; % Quantidade de movimento angular no referencial do corpo
    w0=transpose(C)*wb; % Transformacao do referencial do corpo para o fixo
    H0=transpose(C)*Hb; % Transformacao do referencial do corpo para o fixo
    psiw(i)=atan2(w0(2),w0(1));  % Angulo de guinada da velocidade angular
    thetaw(i)=atan2(w0(3),sqrt(w0(1)^2+w0(2)^2)); % Angulo de elevacao da velocidade angular
    psih(i)=atan2(H0(2), H0(1)); % Angulo de guinada da quantidade de movimento angular
    thetah(i)=atan2(H0(3),sqrt(H0(1)^2+H0(2)^2)); % Angulo de elevacao da quantidade de movimento angular
end

%% Mostra o valor final das velocidades angulares comparando com o inicial
disp('Vetor velocidade angular inícial do VE - rad/s');
disp([wx0a;wy0a;wz0a]);
disp('Vetor velocidade angular final do VE - rad/s');
disp(x(end,1:3));
disp('Vetor velocidade angular relativa final da massa de propelente - rad/s');
disp(x(end,7:9));
%% Graficos
figure
subplot(231);plot(t,x(:,1));grid;xlabel('t (s)');ylabel('\omega_x (rad/s)');axis tight
subplot(232);plot(t,x(:,2));grid;xlabel('t (s)');ylabel('\omega_y (rad/s)');axis tight
subplot(233);plot(t,x(:,3));grid;xlabel('t (s)');ylabel('\omega_z (rad/s)');axis tight
subplot(234);plot(t,x(:,4)*180/pi);grid;xlabel('t (s)');ylabel('\theta_1 (º)');axis tight
subplot(235);plot(t,x(:,5)*180/pi);grid;xlabel('t (s)');ylabel('\theta_2 (º)');axis tight
subplot(236);plot(t,x(:,6)*180/pi);grid;xlabel('t (s)');ylabel('\theta_3 (º)');axis tight

figure
subplot(231);plot(t, phi*180/pi);grid;xlabel('t (s)');ylabel('precessao \phi (º)');axis tight
subplot(232);plot(t,theta*180/pi);grid;xlabel('t (s)');ylabel('nutacao \theta (º)');axis tight
subplot(233);plot(t, psi*180/pi);grid;xlabel('t (s)');ylabel('spin relativo \psi (º)');axis tight
subplot(234);plot(t,x(:,7));grid;xlabel('t (s)');ylabel('s_x (rad/s)');axis tight
subplot(235);plot(t,x(:,8));grid;xlabel('t (s)');ylabel('s_y (rad/s)');axis tight
subplot(236);plot(t,x(:,9));grid;xlabel('t (s)');ylabel('s_z (rad/s)');axis tight

figure
plot(x(:,1),x(:,3));grid;xlabel('\omega_x (rad/s)');ylabel('\omega_z (rad/s)');axis tight
%% Animacao 3D 
% Animação desativada
%{
% Momentos de inercia para calcular a animacao da quantidade de movimento
% angular
I=[Ixx Iyy Izz];
% Insere os angulos na sequencia 313 para gerar a animacao
X(:,4:6)=[phi theta psi];
% Escala de tempo
et=0.1;
% Chama a animacao
%animacao 3d(t,x,l,et);
%}

%% Estuda a resposta para a segunda condicao inicial
% Condicao inicial
theta2=0;theta3=0;
% IMPORTANTE: supondo o vetor quantidade de movimento angular alinhado com
% o eixo n3 do sistema inercial
% Inserir um valor pequeno para calculo de theta1 em t=0, caso contrario, tem singularidade
% quando se calcula o angulo de nutacao  
theta1=0.1*pi/180;
    x0=[wx0b;wy0b;wz0b;theta1;theta2;theta3;0;0;0];
% Integracao da dinamica
% As tolerâncias são importantes, pois pequenos desvios prejudicam os resultados
opt=odeset('RelTol',1e-12,'AbsTol',1e-12,'MaxStep',1);
[t,x] = ode45(@dinamica_VE_semi_rigido, [0 T],x0,opt);
%% Calcula os angulos de nutacao theta, spin relativo psi e precessao phi
phi=x(:,4)+pi/2; % phi=theta1+pi/2
theta=x(:,5); % theta=theta2
psi=x(:,6)-pi/2; % psi=theta3-pi/2
%% Calcula a orientacao dos vetores quantidade de movimento angular e velocidade angular
N=length(t);
phih=zeros(N,1);thetah=phih;psih=phih;
phiw=zeros(N, 1);thetaw=phiw;psiw=phiw;
for i=1:N
C=angle2dcm(x(i,4),x(i,5),x(i,6), 'XYZ'); % Matriz de atitude, transforma do fixo para o corpo
wb=[x(i,1);x(i,2);x(i,3)]; % Velocidade angular no referencial do corpo
Hb=[Ixx*x(i,1);Iyy*x(i,2);Izz*x(i,3)]; % Quantidade de movimento angular no referencial do corpo
w0=transpose(C)*wb; % Transformacao do referencial do corpo para o fixo
H0=transpose(C)*Hb; % Transformacao do referencial do corpo para o fixo
psiw(i)=atan2(w0(2),w0(1));  % Angulo de guinada da velocidade angular
thetaw(i)=atan2(w0(3),sqrt(w0(1)^2+w0(2)^2)); % Angulo de elevacao da velocidade angular
psih(i)=atan2(H0(2),H0(1));  % Angulo de guinada da quantidade de movimento angular
thetah(i)=atan2(H0(3),sqrt(H0(1)^2+H0(2)^2)); % Angulo de elevacao da quantidade de movimento angular
end
%% Mostra o valor final das velocidades angulares comparando com o inicial
disp('Vetor velocidade angular inicial do VE - rad/s');
disp([wx0b;wy0b;wz0b]);
disp('Vetor velocidade angular final do VE - rad/s');
disp(x(end, 1:3));
disp('Vetor velocidade angular relativa final da massa de propelente - rad/s');
disp(x(end,7:9));
%% Graficos
figure
subplot(231);plot(t,x(:,1));grid;xlabel('t (s)');ylabel('\omega_x (rad/s)');axis tight
subplot(232);plot(t,x(:,2));grid;xlabel('t (s)');ylabel('\omega_y (rad/s)');axis tight
subplot(233);plot(t,x(:,3));grid;xlabel('t (s)');ylabel('\omega_z (rad/s)');axis tight
subplot(234);plot(t,x(:,4)*180/pi);grid;xlabel('t (s)');ylabel('\theta_1 (º)');axis tight
subplot(235);plot(t,x(:,5)*180/pi);grid;xlabel('t (s)');ylabel('\theta_2 (º)');axis tight
subplot(236);plot(t,x(:,6)*180/pi);grid;xlabel('t (s)');ylabel('\theta_3 (º)');axis tight
figure
subplot(231);plot(t, phi*180/pi);grid;xlabel('t (s)');ylabel('precessao \phi (º)');axis tight
subplot(232);plot(t,theta*180/pi);grid;xlabel('t (s)');ylabel('nutacao \theta (º)');axis tight
subplot(233);plot(t,psi*180/pi);grid;xlabel('t (s)');ylabel('spin relativo \psi (º)');axis tight
subplot(234);plot(t,x(:,7));grid;xlabel('t (s)');ylabel('s_x (rad/s)');axis tight
subplot(235);plot(t,x(:,8));grid;xlabel('t (s)');ylabel('s_y (rad/s)');axis tight
subplot(236);plot(t,x(:,9));grid;xlabel('t (s)');ylabel('s_z (rad/s)');axis tight
figure
plot(x(:,1),x(:,3));grid;xlabel('\omega_x (rad/s)');ylabel('\omega_z (rad/s)');axis tight
%% Animacao 3D 
%Animação desativada
%{
% Momentos de inercia para calcular a animacao da quantidade de movimento
% angular
I=[Ixx Iyy Izz];
% Insere os angulos na sequencia 313 para gerar a animacao
x(:,4:6)=[phi theta psi];
% Escala de tempo
et=0 1;
% Chama a animacao
%animacao 3d(t,x,|,et);
%}
end

%% Dinamica de um veiculo espacial semi rigido com um reservatorio de propelente esferico
function xp=dinamica_VE_semi_rigido(t,x)
% Vetor de estado: x=[wx;wy;wz;theta1;theta2;theta3;sx;sy;sz]
% wx, wy, wz: velocidades de rotação em torno dos eixos x, y e z do
% sistema de referencia do corpo
% theta1, theta2, theta3: angulos de Euler na sequencia de rotacoes C1, C2 e C3,
% respectivamente. Sequencia 123
% sx, sy, sz: velocidades angulares relativas de uma massa interna de
% propelente
%% Passagem de parametros por variaveis globais
global Ixx Iyy Izz J mu
%% Desmembra o vetor de estado
wx=x(1);wy=x(2);wz=x(3);
theta1=x(4);theta2=x(5);theta3=x(6);
sx=x(7);sy=x(8);sz=x(9);
%% Vetor de controle
[Mx,My,Mz]=controle(t);
%% Equacoes de dinamica
wxp=((Iyy-Izz)*wz*wy+Mx+mu*sx)/(Ixx-J);
wyp=((Izz-Ixx)*wz*wx+My+mu*sy)/(Iyy-J);
wzp=((Ixx-Iyy)*wx*wy+Mz+mu*sz)/(Izz-J);

sxp=-wxp-(mu/J)*sx-wy*sz+wz*sy;
syp=-wyp-(mu/J) *sy-wz*sx+wx*sz;
szp=-wzp-(mu/J)*sz-wx*sy+wy*sx;
%% Equacoes de cinematica
theta1p=(wx*cos(theta3)-wy*sin(theta3))/cos(theta2);
theta2p=wx*sin(theta3)+wy*cos(theta3);
theta3p=(-wx*cos(theta3)*sin(theta2)+wy*sin(theta3)*sin(theta2))/cos(theta2)+wz;
%% Derivada do vetor de estado
xp=[wxp;wyp;wzp;theta1p;theta2p;theta3p;sxp;syp;szp];
end

%% Vetor de controle
function [Mx,My,Mz]=controle(t)
% Sem momento externo resultante
Mx=0;My=0;Mz=0;
end



