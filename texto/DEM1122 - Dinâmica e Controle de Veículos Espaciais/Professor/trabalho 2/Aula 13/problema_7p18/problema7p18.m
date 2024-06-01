% IMPORTANTE: para passar parametros ao Simulink, o programa nao pode ser
% uma funcao, mas um script
%% Problema 7.18 - projeto do controle de arfagem por alocacao de polos
% e sem filtro
%%
clc;clear all;close all;

% Dados
hmax=14e3; % lb.fts
hmax=hmax*0.3048*4.4482216; % N.m
n=0.0011; % rad/s

% Requisito - polos de malha fechada desejados
pr(1) = -1.5*n;
pr(2) = -1*n;
pr(3) = (-1.5+1i*1.5)*n;
pr(4) = (-1.5-11*1.5)*n;

% Chama a funcao da dinamica
[A,B,B2,ftarf,ftrg] = modelo_iss_cmg;
% Matrizes do modelo de espaco de estados de ordem reduzida de arfagem
%  1      2     3   4      5     6  7    8   9
% [wx Dwy wz phi theta psi hx hy hz] -> [thetap theta hy] (Dwy=thetap)
%  1    2   3
% [ux uy uz] -> [uy].
Aa=[A(2,2) A(2,5) A(2,8)
        A(5,2) A(5,5) A(5,8)
        A(8,2) A(8,5) A(8,8)];

Ba=[B(2,2)
        B(5,2)
        B(8,2)];

B2a=[B2(2,2)
          B2(5,2)
          B2(8,2)];
% Insere uma variavel de estado adicional referente a integral do momento
% angular hy
Aa=[Aa zeros(3,1)
        0 0 1 0];
Ba = [Ba;0];
B2a = [B2a;0];
% Matriz de ganhos de realimentacao de estado determinada por alocacao de
% polos
K=place(Aa,Ba,pr);
% Desmembra os ganhos do controlador - o diagrama usa o controle com sinal
% positivo ao inves de negativo, necessario trocar nos ganhos.
kyp=-K(2) % Realimentacao de theta
kyd=-K(1) % Realimentacao de thetap
kyh=-K(3) % Realimentacao de hy
kyi=-K(4) % Realimentacao da int(hy)
% Matriz A de malha fechada
Ac=Aa-Ba*K;
% Polos de malha fechada
damp(Ac)
% Modelo de espaco de estados de malha fechada - saidas theta e hy, entrada
%dy
C=[0 1 0 0;0 0 1 0]; 
D=[0;0];
mfarf=ss(Ac,B2a,C,D);

% Funcoes de transferencia de malha fechada
mfftarf=zpk(mfarf)

% Diagrama de bode da saída theta com respeito a entrada dy em malha
% fechada

figure();
bode(mfftarf(1,1));grid;

% Diagrama de bode da saída hy com respeito a entrada dy em malha
% fechada
figure();
bode(mfftarf(2,1));grid;
% Entradas do simulink

global T disturbio_y

T=5*2*pi/n; % 4 orbitas - intervalo da simulacao
t=0:30:T; % Vetor de tempos
t=t'; % Vetor coluna
N=length(t);
d=zeros(N,1); % Vetor de perturbacao dy em funcao do tempo

for j=1:N
    d(j)=4+2*sin(n*t(j))+0.5*sin(2*n*t(j)); % foot pound
    d(j)=d(j)*0.3048*4.4482216; % N.m
    
end

% Estrutura enviada ao Simulink
disturbio_y.time=t;
disturbio_y.signals.values=d;
% Condicoes iniciais

theta0=1*pi/180;
thetap0=0.001*pi/180;

% Simula

sim('diagrama_problema7p18');

% Resultados
t=theta.time;
N=length(t);
theta=theta.signals.values;
hy=hy.signals.values;
u=u.signals.values;
dy=dy.signals.values;
% Graficos
figure()
subplot(221);plot(t,theta*180/pi);grid;xlabel('(s)');ylabel('\theta (º)')
subplot(222);plot(t,hy,t,hmax*ones(1,N),t, -hmax*ones(1,N));grid;xlabel('t (s)');ylabel('h_y (N.m.s)');
subplot(223);plot(t,dy);grid;xlabel('t (s)');ylabel('d_y (N.m)');
subplot(224);plot(t,u);grid;xlabel('t (s)');ylabel('u_y (N.m)');

