function problema7p14
% Dinamica e Controle de Veiculos Espaciais
% Problema 7.14 da referencia Wie, B. Space Vehicle Dynamics and Control. 2 . ed., AJAA
% Education Series, Reston, VA: AIAA, 2008.
clc;close all;clear all;
%% Entrada de parametros
global J tipo k c1 c2 c3 a b vu tu
% Matriz de momentos de inercia
J=[1200 100 -200
     100 2200 300
    -200 300 3100];  % kg.m'2
% Condicao inicial
Q0=[0.5;0.5;0.5;-0.5];
w0=[0;0;0];
% Intervalo de simulacao
T=500; %s
%% Projeto do controle
tipo=input('Que tipo de matriz de ganhos pretende usar? 1, 2, 3, ou 4?');
switch tipo
    case 1
        disp('K=k*I33, C=diag(c1,c2,c3) - Digite os valores: ');
        k=input('k: ');c1=input('c1: ');c2=input('c2: ');c3=input('c3: ');
    case 2
        disp('K=(k/g4^3)*I33, C=diag(c1,c2,c3) - Digite os valores: ');
        k=input('k: ');c1=input('c1: ');c2=input('c2: ');c3=input('c3: ');
    case 3
        disp('K=k*sgn(q4)*I33, C=diag(c1,c2,c3) - Digite os valores: ');
        k=input('k: '); c1=input('c1: ');c2=input('c2: ');c3=input('c3:');
    case 4
        disp('K=(a*diag(J)+b*I33)^-1, C=diag(c1,c2,c3) - Digite os valores: ');
        a=input('a>0: ');b=input('b>0: ');c1=input('c1>0: ');c2=input('c2>0: ');c3=input('c3>0: ');
    otherwise
        disp('Digite um valor de tipo valido.');
end
%% Integracao da dinamica
% Condicao inicial
x0=[w0;Q0];
% Vetores para guardar os controles e o tempo associado
vu=[0 0 0];
tu=0;
% As tolerâncias são importantes, pois pequenos desvios prejudicam os resultados
opt=odeset('RelTol', 1e-12,'AbsTol',1e-12,'MaxStep',1);
[t,x]=ode45(@dinamica_VE_quat,[0 T],x0,opt);
% Angulo principal
PHI=2*acos(x(:,7));
%% Graficos
figure
subplot(421);plot(t,x(:,1));grid;xlabel('t (s)');ylabel('\omega_x (rad/s)');axis tight
subplot(423);plot(t,x(:,2));grid;xlabel('t (s)');ylabel('\omega_y (rad/s)');axis tight
subplot(425);plot(t,x(:,3));grid;xlabel('t (s)');ylabel('\omega_z (rad/s)');axis tight
subplot(422);plot(t,x(:,4));grid;xlabel('t (s)');ylabel('q_1');axis tight

subplot(424);plot(t,x(:,5));grid;xlabel('t (s)');ylabel('q_2');axis tight
subplot(426);plot(t,x(:,6));grid;xlabel('t (s)');ylabel('q_3');axis tight
subplot(428);plot(t,x(:,7));grid;xlabel('t (s)');ylabel('q_4');axis tight
figure
subplot(311);plot(tu,vu(:,1));grid;xlabel('t (s)');ylabel('u_x - N.m');axis tight
subplot(312);plot(tu,vu(:,2));grid;xlabel('t (s)');ylabel('u_y - N.m');axis tight
subplot(313);plot(tu,vu(:,3));grid;xlabel('t (s)');ylabel('u_z - N.m');axis tight
figure
subplot(221);
plot(t, PHI*180/pi);grid;xlabel('t (s)');ylabel('angulo principal \Phi (graus)');axis tight
subplot(222);plot(x(:,4),x(:,5));grid;xlabel('q_1');ylabel('q_2');axis tight
subplot(223);plot(x(:,4),x(:,6));grid;xlabel('q_1');ylabel('q_3');axis tight
subplot(224);plot(x(:,5),x(:,6));grid;xlabel('q_2');ylabel('q_3');axis tight

%% Estudo de angulos de Euler
% Converte de quaternion para os angulos de Euler na sequencia 123
% Tomar cuidado, usamos o quaternion com a parte real no elemento 4,
% mas o MATLAB usa no elemento 1
[the1,the2,the3]=quat2angle([x(:,7),x(:,4),x(:,5) x(:,6)], 'XYZ');
% Calcula os angulos de nutacao: theta, spin relativo: psi e precessao: phi
phi=the1+pi/2; % phi=theta1+pi/2
theta=the2; % theta=theta2
psi=the3-pi/2; % psi=theta3-pi/2
% Graficos
figure
subplot(231);plot(t,the1*180/pi);grid;xlabel('t (s)');ylabel('\theta_1 (º)');axis tight
subplot(232);plot(t,the2*180/pi);grid;xlabel('t (s)');ylabel('\theta_2 (º)');axis tight
subplot(233);plot(t,the3*180/pi);grid;xlabel('t (s)');ylabel('\theta_3 (º)');axis tight
subplot(234);plot(t, phi*180/pi);grid;xlabel('t (s)');ylabel('\phi (º)');axis tight
subplot(235);plot(t,theta*180/pi);grid;xlabel('t (s)');ylabel('\theta (º)');axis tight
subplot(236);plot(t, psi*180/pi);grid;xlabel('t (s)');ylabel('\psi (º)');axis tight

%% Animacao 3D
%ANIMACAO 3D DESATIVADA
%{
% Insere os angulos na sequencia 313 para gerar a animacao

x(:,4:6)=[phi theta psi];

% Escala de tempo

et=10;

% Chama a animacao

animacao 3d(tx,J,et);
%}
end

%% Dinamica de um veiculo espacial rigido com cinematica pelos quaternions

function xp=dinamica_VE_quat(t,x)
% Vetor de estado: x=[wx;wy;wz;q1;q2;q3;q4]
% wx, wy, wz: velocidades de rotação em torno dos eixos x, y e z do SRC
% q1, q2, q3, q4: Componentes do quaternion de atitude do referencial do
% corpo com respeito ao referencial inercial

%% Passagem de parametros por variaveis globais

global J vu tu

%% Desmembra o vetor de estado
w=[x(1);x(2);x(3)];
Q=[x(4);x(5);x(6);x(7)];
%% Vetor de controle
u=controle_real_quat(x);
% Salva o historico do controle em uma variavel global
vu=[vu;u'];
tu=[tu;t];
%% Equacoes de dinamica
wp=J^(-1)*(u-skew(w)*J*w);
%% Equacoes de cinematica de quaternion
Qp=0.5*[-skew(w) w;-w' 0]*Q;
%% Derivada do vetor de estado
xp=[wp;Qp];
end
%% Vetor de controle
function u=controle_real_quat(x)
% O controle objetiva regular o satelite para o equilibrio, ou seja, para
% velocidade angular nula, com o referencial do corpo alinhado ao
% referencial inercial. Assim, a referencia do quaternion eh o quaternion
% identidade e o erro de rastreio eh o proprio quaternion medido.
% A estrutura do controle eh generica, consistindo de um analogo de
% controle proporcional derivativo: u=-K*q-C*w ("q" eh a parte vetorial do
% quaternion e "w" a velocidade angular com respeito ao referencial
% inercial).
% De acordo com o conteudo da secao 7.3.1, sao 4 possibilidades sugeridas
% para montar as matrizes de realimentacao K e C. Que sao identificadas
% pela variavel "tipo"
global tipo k c1 c2 c3 J a b
% Desmembra a entrada
w=[x(1);x(2);x(3)]; % Velocidade angular
q=[x(4);x(5);x(6)]; % Parte vetorial do quaternion
q4=x(7); % Parte escalar do quaternion
% Escolhe a forma de calculo dos ganhos
switch tipo
    case 1
        K=k*eye(3);
        C=diag([c1,c2,c3]);
    case 2
        K=(k/q4^3)*eye(3);
        C=diag([c1,c2,c3]);
    case 3
        K=k*sign(q4)*eye(3);
        C=diag([c1,c2,c3]);
    case 4
        % Usar somente os momentos de inercia no calculo dos controles
        K=(a*diag([J(1,1),J(2,2),J(3,3)])+b*eye(3))^(-1);
        C=diag([c1,c2,c3]);
end
% Lei de controle
u=-K*q-C*w;
end
%% Matriz anti simetrica do produto vetorial

function S = skew(w)
    S=[0    -w(3) w(2)
         w(3) 0     -w(1)
        -w(2) w(1) 0];
end