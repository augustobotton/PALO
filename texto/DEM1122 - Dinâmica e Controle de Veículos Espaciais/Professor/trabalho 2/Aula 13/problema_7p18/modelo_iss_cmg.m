function [A,B,B2,ftarf,ftrg] = modelo_iss_cmg
% Dinamica e Controle de Veiculos Espaciais
% Secao 7.4 da referencia Wie, B. Space Vehicle Dynamics and Control. 2. ed., AJAA
% Education Series, Reston, VA: AIAA, 2008.
% Modelo da ISS com atuacao de atitude por CMG
% Saidas
% A,B,B2: matrizes da dinamica linearizada para a condicao de equilibrio
% nadir apontado em malha aberta (controle e perturbacoes nulos, matriz de inercia
% diagonal)
% farf, ftrg: matrizes de funcao de transferencia dos eixos de arfagem e
% rolamento/guinada, respectivamente.
clc;close all;clear al|;
%% Entrada de parametros
global II n 
% Matriz de momentos de inercia
J=[50.28, -0.39, 0.16
    -0.39, 10.80, 0.16
     0.16, 0.16, 58.57]*10^6; % slug.f2
J=J*14.5939*0.3048^2;
% Movimento medio da orbita
n=0.0011; %rad/s

%% Equilibrio e linearizacao
% Matriz de momentos de inercia passada como parametro para as funcoes
II = J;
per = input('Deseja calcular o equilibrio para a parcela constante da pertubacao atuante? sim (1), nao (0): ');
if per
    d=perturbacao_ISS(0);
    [Xe,Ue,hp]=eqOrbCirc(d);
    fprintf('Equilibrio para controle nulo e perturbacao d_x = %f N.m, d_y = %f N.m, d_z = %f N.m\n',...
    d(1),d(2),d(3));
    fprintf('Velocidades angulares: w_x=%f rad/s, w_y=%f rad/s, w_z=%f rad/s, W',Xe(1),Xe(2),Xe(3));
    fprintf('Angulos de atitude com respeito aa LVLH: phi=%f, theta=%f, psi=%f \n',...
    Xe(4)*180/pi,Xe(5)*180/pi,Xe(6)*180/pi);
    fprintf('Derivada do momento angular do CMG: hp_x=%f N.m/s, hp_y=%f N.m/s, hp_z=%f N.m/s \n',...
    hp(1),hp(2),hp(3));
end

% Determinacao do equilibrio da atitude para perturbacao e controle nulos
d=[0;0;0];
[Xe,Ue,hp] = eqOrbCirc(d);
disp('***** Equilbrio para perturbacao e controle nulos *****');
fprintf('Velocidades angulares: w_x=%f rad/s, w_y=%f rad/s, w_z=%f rad/s, \n',Xe(1),Xe(2),Xe(3));
fprintf('Angulos de atitude com respeito ao LVLH: phi=%f, theta=º%f, psi=%f, \n',...
Xe(4)*180/pi,Xe(5)*180/pi,Xe(6)*180/pi);
fprintf('Derivada do momento angular do CMG: hp_x=%f N.m/s, hp_y=9%f N.m/s, hp_z=%f N.m/s, \n',...
hp(1),hp(2),hp(3));

% Linearizacao
[A,B,B2]=linearizar(Xe,Ue,d);
fprintf('\n Matriz A\n ');
fprintf('\t%0.2e\t%0.2e\t%0.2e\t%0.2e\t%0.2e\t%0.2e\t%0.2e\t%0.2e\t%0.2e\n',A');
fprintf('\n Matriz B\n');
fprintf('\t%0.2e\t%0.2e\t%0.2e\n', B');
fprintf('\n Matriz B_2 \n');
fprintf('\t%0.2e\t%0.2e\t%0.2e\n', B2');

%% Aproximacao de matriz de inercia diagonal
% Matriz de momentos de inercia passada como parametro para as funcoes
II = diag( [J(1,1),J(2,2),J(3,3)] );
disp('**** Aproximacao de matriz de inercia diagonal *****');
per=input( 'Deseja calcular o equilibrio para a parcela constante da pertubacao atuante? sim (1), nao (0): ');
if per
    d=perturbacao_ISS(0);
    [Xe,Ue,hp]=eqOrbCirc(d);
    fprintf('Equilibrio para controle nulo e perturbacao d_x = %f N.m, d_y= %f N.m, d_z= %f N.m\n',...
    d(1),d(2),d(3));
    fprintf('Velocidades angulares: w_x=%f rad/s, w_y=%f rad/s, w_z=%f rad/s, \n',Xe(1),Xe(2),Xe(3));
    fprintf('Angulos de atitude com respeito ao LVLH: phi=%f, theta=%f, psi=%f, \n',...
    Xe(4)*180/pi,Xe(5)*180/pi,Xe(6)*180/pi);
    fprintf('Derivada do momento angular do CMG: hp_x=%f N.m/s, hp_y=%f N.m/s, hp_z=%f N.m/s, \n',...
    hp(1),hp(2),hp(3));
end

% Determinacao do equilibrio da atitude para perturbacao e controle nulos

d=[0;0;0];

[Xe,Ue,hp]=eqOrbCirc(d);

disp('***** Equilibrio para perturbacao e controle nulos ****');

fprintf('Velocidades angulares: w_x=%f rad/s, w_y=%frad/s, w_z=%f rad/s, \n',Xe(1),Xe(2),Xe(3));
fprintf('Angulos de atitude com respeito ao LVLH: phi=%f, theta=%f, psi=%f, \n',...
Xe(4)*180/pi,Xe(5)*180/pi,Xe(6)*180/pi);

fprintf('Derivada do momento angular do CMG: hp_x=%f N.m/s, hp_y=%f N.m/s, hp_z=%f N.m/s, \n',...
hp(1),hp(2),hp(3));

[A,B,B2]=linearizar(Xe,Ue,d);

fprintf('\n Matriz A\n ');
fprintf('\t%0.2e\t%0.2e\t%0.2e\t%0.2e\t%0.2e\t%0.2e\t%0.2e\t%0.2e\t%0.2e\n',A');

fprintf('in Matriz Bin');
fprintf('t%0.2e\t%0.2e\t%0.2e\n', B');

fprintf('\n Matriz B_2 \n');
fprintf('t%0.2e\t%0.2e\t%0.2e\n', B2');

%% Estabilidade
damp(A)
%% Modelos de espaco de estado e funcoes de transferencia

% Eixo de arfagem

b=B(:,2); % Entada u y

C=[0 0 0 0 1 0 0 0 0;0 0 0 0 0 0 0 1 0]; D=[0;0]; % Saidas theta e hy
ssarf=ss(A,b,C,D);

disp('*** Funcoes de transferencia do eixo de arfagem - saídas theta e hy, entrada uy ***');
zpk(ssarf)
ftarf=tf(ssarf);

% Eixos de rolamento e guinada

b=[B(:,1) B(:,3)]; % Entadas u_x e u_z

% Saidas phi, psi, hx e hz
C=[0 0 0 1 0 0 0 0 0;0 0 0 0 0 1 0 0 0;0 0 0 0 0 0 1 0 0;0 0 0 0 0 0 0 0 1];
D=[0 0;0 0;0 0;0 0];

ssrg = ss(A,b,C,D);

disp('*** Funcoes de transferencia dos eixos de rolamento/guinada - saidas phi, psi, hx e hz, entradas ux e uz ***');
zpk(ssrg)
ftrg=zpk(ssrg);
end

%% Dinamica de atitude de um veiculo espacial em orbita circular baixa
% com torque de gradiente gravitacional e atuador CMG idealizado (sem
% dinamica dos guimbais)
function xp=dinamica_VE_tgg_cmg(t,x,u,d)
% Vetor de estado: x=[wx;wy;wz;phi;theta;psi;hx;hy;hz]
% wx, wy, wz: velocidades de rotação inerciais em tomo dos eixos x, y e z do SRC
% phi, theta, psi [rad]: Angulos de Euler (rolamento, arfagem e guinada) do
% SRC com respeito ao LVLH de acordo com a sequencia 321
% hx, hy, hz [N.m.s]: quantidade de movimento angular do CMG nos eixos do SRC
% u [N.m]: vetor de controle, torques nos eixos x, y e z do SCR aplicados
% pelo atuador
% d [N.m]: vetor de disturbio, torques nos eixos x, y e z do SCR aplicados
% por alguma causa externa
%%
% Passagem de parametros por variaveis globais
global II n
% Desmembra o vetor de estado
w=[x(1);x(2);x(3)];
phi=x(4);theta=x(5);psi=x(6);
h=[x(7);x(8);x(9)]; % Vetor quantidade de movimento angular do atuador
% Matriz de rotacao do SRC com respeito ao LVLH
C1 = [1, 0, 0;0, cos(phi), sin(phi);0, -sin(phi), cos(phi)];
C2 = [cos(theta), 0, -sin(theta);0, 1, 0;sin(theta), 0, cos(theta)];
C3 = [cos(psi), sin(psi), 0;-sin(psi), cos(psi), 0;0, 0, 1];
C= C1*C2*C3;
% Velocidade relativa do SRC com respeito ao LVLH escrita no SRC
wlvlhc = C*[0;-n;0]; % Velocidade do LVLH escrita no SRC
wblvlh = w-wlvlhc;
% Equacoes de cinematica de angulos de Euler 321
Angp = (1/cos(theta))*[cos(theta) sin(phi)*sin(theta) cos(phi)*sin(theta)
                0 cos(phi)*cos(theta) -sin(phi)*cos(theta)
                0 sin(phi) cos(phi)]*wblvlh;
% Torque de gradiente gravitacional em orbita circular
Mg=TGG(C,n);
% Equacao de dinamica de rotacao do satelite. Equacao classica do corpo rigido,
% mais torque de gradiente gravitacional, acao de controle "u" aplicada pelo
% atuador e disturbio "d”
wp = II^(-1)*(-skew(w)*II*w+Mg-u+d);
% Dinamica do CMG ideal (ignorando o tempo de resposta do driver de controle
% e a resposta dos guimbais). Neste caso, se assemelha as rodas de reacao
hp = -skew(w)*h+u;
% Derivada do vetor de estado
xp=[wp;Angp;hp];
end

%% Perturbacao
function d=perturbacao_ISS(t)
global n
% Torque perturbativo atuante sobre a ISS, resultante da acao aerodinamica
    d= [1+sin(n*t)+0.5*sin(2*n*t)
    4+2*sin(n*t)+0.5*sin(2*n*t)
    1+sin(n*t)+0.5*sin(2*n*t)]; % foot.pound
    d = d*0.3048*4.4482216; % N.m
end

%% Torque de gradiente gravitacional
function Mg=TGG(C,n)
% Funcao para calculo de torque de gradiente gravitacional em orbita
% circular
% Entradas
% C: Matriz de rotacao do SRC com respeito ao LVLH
% n [rad/s]: movimento medio da orbita circular
%%
% Passagem de parametros por variaveis globais
global II
% Torque de gradiente gravitacional no SRC.
Mg=3*n^2*cross(C(:,3),II*C(:,3));
end

%% Matriz anti simetrica do produto vetorial
function S=skew(w)
    S=[0 -w(3) w(2)
         w(3) 0 -w(1)
        -w(2) w(1) 0];
end

%% Funcao de linearizacao do modelo
function [A,B,B2]=linearizar(Xe,Ue, De)
% Calcula as matrizes do modelo linear em torno de um ponto de equilibrio
% Entradas:
% Xe: vetor de estado no equilibrio
% Ue: vetor de controle no equilibrio
% De: vetor de disturbio no equilibrio
% Saidas:
% A: matriz de estado do modelo linear
% B: matriz de controle do modelo linear
% B2: matriz do disturbio do modelo linear
%%
% Incrementos para calcular diferenças finitas das variáveis de estado
dx=1e-3;inc=dx*eye(9,9);
% Inicializa a matriz À
A=zeros(9,9);
% Função de estado no equilíbrio
fe = dinamica_VE_tgg_cmg(0,Xe,Ue,De);
% Calcula a matriz "A” a partir de diferenças finitas das variáveis de estado
for j=1:9
    f = dinamica_VE_tgg_cmg(0,Xe+inc(:,j),Ue,De);
    A(:,j) = (f-fe)/dx;
end
% Incrementos para calcular as diferenças finitas das variáveis de controle
du=1e-2; inc=du*eye(3,3);
% Inicializa a matriz B
B=zeros(9,3);
% Calcula a matriz "B" a partir de diferenças finitas das variáveis de controle
for j=1:3
    f=dinamica_VE_tgg_cmg(0,Xe,Ue+inc(:,j),De);
    B(:,j)=(f-fe)/du;
end

% Incrementos para calcular as diferenças finitas das variáveis de disturbio
dd = 1e-2; inc=dd*eye(3,3);

% Inicializa matriz B2
B2=zeros(9,3);
% Calcula a matriz "B2" a partir de diferenças finitas das variáveis de disturbio

for j=1:3
    f=dinamica_VE_tgg_cmg(0,Xe,Ue,De+inc(:,j));
    B2(:,j)=(f-fe)/dd;
end
end

%% Funcao para calcular o equilibrio da atitude do satelite

function [Xe,Ue,hp]=eqOrbCirc(d)

% Calcula o equilibrio da atitude do satelite em orbita baixa submetido a torque de
% gradiente gravitacional e disturbio constante nos eixos do SRC.
% Entradas:
% d [N.m]: Disturbio constante nos eixos do SRC
% Saidas:
% Xe: vetor de estado no equilibrio
% Ue: vetor de controle no equilibrio
% hp: derivada da quantidade de movimento angular do CMG
% Condicoes de calculo de equilibrio
% O controle no equilibrio eh zero - sem uso/do atuador
% Sao consideradas as equacoes diferenciais do corpo somente, as do atuador
% sao ignoradas (variaveis ignoraveis)

%%
% Passagem e recebimento de parâmetros por variáveis globais
global D n

% Passagem de parâmetros para a função objetivo
D=d;
% Chute inicial para as incógnitas
% z=[wx wy wz phi theta psi]
z0=[0 -n 0 0 0 0 0]';% Solucao do caso sem disturbio

% Chamada da funcao objetivo
options=optimset('TolFun',1e-20,'MaxFunEvals',1e6, 'MaxIter',1e5,'TolX',1e-20);
[z, fval, flag]=fsolve(@obj_eq,z0,options);
disp('Valor final da função objetivo:');disp(fval);
disp('flag');disp('flag');

% Resultado
% X=[wx wy wz phi theta psi hx hy hz]
Xe=[z(1) z(2) z(3) z(4) z(5) z(6) 0 0 0]'; % Valores nulos das variaveis ignoraveis
% U=[ux uy uz]'
Ue = [0 0 0]'; % O calculo do equilibrio eh feito sob a condicao de controle nulo
% Derivada da quantidade de movimento angular do CMG no equilibrio
Xp=dinamica_VE_tgg_cmg(0,Xe,Ue,D);
hp=Xp(7:9);
end
%% Função objetivo para cálculo do equilibrio
function f=obj_eq(z)
% Recebimento de parâmetros por variáveis globais
global D
% Montagem dos vetores de estado e controle a partir das incognitas,
% parametros fornecidos e condicoes de contorno
% z=[wx wy wz phi theta psil'
% X=[wx wy wz phi theta psi hx hy hz]
Xe =[z(1) z(2) z(3) z(4) z(5) z(6) 0 0 0]';
% Us=[ux uy uz]'
Ue=[0 0 0]'; % O calculo do equilibrio eh feito sob a condicao de controle nulo
% Calcula a função de estado (derivadas das variaveis de estado)
Xp = dinamica_VE_tgg_cmg(0,Xe,Ue,D);
% Elimina as derivadas das variaveis ignoraveis
f=Xp(1:6);
end
