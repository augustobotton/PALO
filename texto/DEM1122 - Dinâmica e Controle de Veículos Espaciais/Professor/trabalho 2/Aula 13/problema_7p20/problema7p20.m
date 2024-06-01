% IMPORTANTE: para passar parametros ao Simulink, o programa nao pode ser
% uma funcao, mas um script
clc;clear all;close all;
% Dados
hmax=14e3; % lb.ft.s
hmax=hmax*0.3048*4.4482216; % N.m
n=0.0011; % rad/s
% Requisito - polos de malha fechada desejados
% pr(l)=(-1.05+1i*0.68) *n;pr(2)=(-1.05-1i*0.68) *n;pr (3) =(-1.04+1i*0.72) *n;pr (4)=(-1.04-1i*0.72) *n;

% pr(5)=(-1.42+1i*1.38) *n;pr (6)=(—1.42-1i*1.38) *n;pr (7)=(-1.42+1i*1.38) *n;pr (8)=(-1.42-1i*1.38) *n;

pr(1)=(-1.05+1i*0.68)*n;
pr (2)=(-1.05-1i*0.68)*n;
pr (3) =(-1.04+1i*0.72)*n;
pr (4)=(-1.04-1i*0.72) *n;

pr(5)=(-1.42+1i*1.38) *n;
pr(6)=(-1.42-1i*1.38) *n;
pr(7)=(-1.3+1i*1.38) *n;
pr (8)=(-1.3-1i*1.38)*n;

% Chama a funcao da dinamica
[A,B,B2,ftarf,ftrg]=modelo_iss_cmg;
% Matrizes do modelo de espaco de estados de ordem reduzida de

% rolamento/guinada
%   1    2     3    4     5     6   7  8   9
% [wx Dwy wz phi theta psi hx hy hz)] -> [wx wz phi psi hx hz]

%   1   2   3
% [ux uy uz) -> [ux uz].

Arg=A;Arg(:,8)=[];Arg(8,:)=[];Arg(:,5)=[];Arg(5,:)=[];Arg(:,2)=[];Arg(2,:)=[];
Brg= B;Brg(:,2)=[];Brg(8,:)=[];Brg(5,:)=[];Brg(2,:)=[];
B2rg= B2;B2rg(:,2)=[];B2rg(8,:)=[];B2rg(5,:)=[];B2rg(2,:)=[];

% Insere 2 variaveis de estado adicionais referentes as integrais dos
% momentos angulares hx e hz
Arg=[Arg zeros(6,2)
         0 0 0 0 1 0 0 0
         0 0 0 0 0 1 0 0];
Brg=[Brg;0 0;0 0];
B2rg=[B2rg;0 0;0 0];
%% Verifica a controlabilidade
% Co=ctrb(Arg,Brg)
%% Numero de estados incontrolaveis
% inco=length(Arg)-rank(Co,le-20)
% Matriz de ganhos de realimentacao de estado determinada por alocacao de
% polos
K = place (Arg,Brg,pr);
%           1    2   3   4   5   6     7       8
% 1 ux: wx wz phi psi hx hz int.hx int.hz
% 2 uz: wx wz phi psi hx hz int.hx int.hz
% Desmembra os ganhos do controlador - o diagrama usa o controle com sinal
% positivo ao inves de negativo, necessario trocar nos ganhos.
K11=-[K(1,3) K(1,1) K(1,5) K(1,7)] % Realimentacao de phi, wx, hx, int.hx (nesta ordem) para ux
K13=-[K(1,4) K(1,2) K(1,6) K(1,8)] % Realimentacao de psi, wz, hz, int.hz (nesta ordem) para ux

K31=-[K(2,3) K(2,1) K(2,5) K(2,7)] % Realimentacao de phi, wx, hx, int.hx (nesta ordem) para uz

K33=-[K(2,4) K(2,2) K(2,6) K(2,8)] % Realimentacao de psi, wz, hz, int.hz (nesta ordem) para uz

% Matriz A de malha fechada

Ac=Arg-Brg*K;

% Polos de malha fechada
damp (Ac )
% Modelo de espaco de estados de malha fechada - saidas phi, hx, psi e hz,
% entradas dx e dz

C=[0 0 1 0 0 0 0 0;0 0 0 0 1 0 0 0;0 0 0 1 0 0 0 0;0 0 0 0 0 1 0 0];

D=[0 0;0 0;0 0;0 0];

mfrg=ss(Ac,B2rg,C,D);
% Funcoes de transferencia de malha fechada
mfftrg=zpk(mfrg)
% Diagrama de bode da saida phi com respeito a entrada dx em malha fechada
figure ();bode(mfftrg(1,1));grid;title('Entrada d_x saida \phi');
% Diagrama de bode da saida hx com respeito a entrada dx em malha fechada
figure();bode(mfftrg(2,1));grid;title('Entrada d_x saida h_x');

% Diagrama de bode da saida psi com respeito a entrada dz em malha fechada
figure();bode(mfftrg (3,2) );grid;title('Entrada d_z saida \psi');
% Diagrama de bode da saida hz com respeito a entrada dz em malha fechada
figure();bode(mfftrg(4,2) );grid;title('Entrada d_z saida h_z');

% Entradas do simulink
T=5*2*pi/n; % 4 orbitas - intervalo da simulacao
t=0:30:T; % Vetor de tempos
t=t'; % Vetor coluna
N=length(t);
d=zeros(N,2); % Vetor de perturbacao dy em funcao do tempo

for j=1:N
    d(j,1)=1+sin(n*t(j))+0.5*sin(2*n*t(j)); % foot.pound
    d(j,2)=1+sin(n*t(j))+0.5*sin(2*n*t(j));
    d(j,:)=d(j,:)*0.3048*4.4482216; % N.m
end

% Estrutura enviada ao Simulink

disturbios_xz.time=t;
disturbios_xz.signals.values=d;
% Condicoes iniciais
wx0=0.001*pi/180;
wz0=0.001*pi/180;
phi0=1*pi/180;
psi0=1*pi/180;

% Simula no Simulink
sim  = sim('diagrama_7p20');

% Resultados

t = sim.estado_rg.time;
N = length(t);
wx_wz_phi_psi = sim.estado_rg.signals.values;
hxz = sim.hxz.signals.values;
u = sim.u.signals.values;
dxz = sim.dxz.signals.values;
%% Graficos
figure()
subplot (231);plot (t,wx_wz_phi_psi(:,2)*180/pi);grid;xlabel('t (s)');ylabel('\omega_x (º/s)');
subplot (232);plot (t,wx_wz_phi_psi(:,1)*180/pi) ;grid;xlabel ('t (s)');ylabel('\phi (º)')
subplot (233);plot (t,hxz(:,1),t,hmax*ones(1,N),t,-hmax*ones(1,N));grid;xlabel ('t (s)');ylabel('h_x (N.m.s)')
subplot (234);plot (t,wx_wz_phi_psi(:,4)*180/pi) ;grid;xlabel('t (s)');ylabel('\omega_z (“/s)')
subplot (235);plot (t,wx_wz_phi_psi(:,3)*180/pi) ;grid;xlabel('t (s)');ylabel('\psi (º)')
subplot (236);plot (t,hxz (:,2),t,hmax*ones (1,N),t,-hmax*ones(1,N)) ;sgrid;xlabel ('t (s)');ylabel('h_z (N.m.s)')

figure()
subplot (221) ;plot (t,dxz (:,1)) ;grid;xlabel ('t (s)');ylabel('d_x (N.m)');
subplot (222) ;plot (t,dxz (:,2)) ;grid;xlabel ('t (s)');ylabel('d_z (N.m)');
subplot (223) ;plot (t,u(:,1)) ;grid;xlabel('t (s)');ylabel('u_x (N.m)');
subplot (224) ;plot (t,u(:,2));grid;xlabel ('t (s)');ylabel('u_z (N.m)');