% Script para simular um voo de sondagem de veiculo de dois estagios com
% carga util.
% Eh utilizado o foguete de sondagem VSB30 com os dados abaixo, retirados
% de referencias indicadas em aula.
close all;clear all;clc;
global Re we mut J2 J3 J4 g lc dT Sr fc
global ms m0 mp  ti tq ts Isp h0 l_trilho
%% Parametros propulsivos
Isp(1)=260.6;     % s - Impulso específico do primeiro estágio (motor S31)
Isp(2)=261.1;     % s - Impulso específico do segundo estágio (motor S30)
mp(1)=677;     % kg - Massa de propelente do primeiro estagio (motor S31)
mp(2)=898;     % kg - Massa de propelente do segundo estagio (motor S30)
ti(1)=0; % s - Tempo da ignicao do primeiro estagio
ti(2)=15; % s - Tempo da ignicao do segundo estagio
tq(1)=13.5;    % s - Tempo do fim da queima do estágio 1 (motor S31)
tq(2)=44;       % s - Tempo do fim da queima do estágio 2 (motor S30)
ts(1)=13.5; % Tempo da separacao do primeiro estagio
ts(2)=59; % Tempo da separacao do segundo estagio
%% Parametros de massa estrutural e de carga util
ms(1)=284; % kg - Massa estrutural do primeiro estagio
ms(2)=320; % kg - Massa estrutural do segundo estagio
mL=400; % kg - Massa da carga util

%% Parametros aerodinamicos e ambientais
% Fator de correcao do arrasto a partir de dados de tunel de vento
fc=1.28;
% As areas de referencia sao das secoes transversais
S1=pi*(0.557/2)^2; % m^2 - Area aproximada da secao transversal do primeiro estagio
S2=pi*(0.557/2)^2; % m^2 - Area aproximada da secao longitudinal do segundo estagio
SL=pi*(0.46/2)^2; % m^2 -  Area aproximada da secao longitudinal da carga util
% Correcoes para levar em conta a area molhada, feitas em funcao do
% comprimento. Assume-se que a correcao soh ocorre em 50% da area de
% referencia, supondo uma influencia de 50% da area molhada no arrasto
% total
lt=12.6; % m - Comprimento total
l2=lt-3.214;    %  Comprimento sem o primeiro estagio
l3=l2-3.294; % Comprimento da carga util
f2=(l2/lt)*0.5+0.5; % Fator de correcao do segundo estagio
f3=(l3/lt)*0.5+0.5; % Fator de correcao da carga util
% Vetor de areas de referencia para calculo do arrasto
Sr=[S1
        S2*f2
        SL*f3];
lc=0.509; % Comprimento característico 
dT=10;  % K - Delta T em relação à atmosfera padrão (que é 15ºC no nível do mar)

%% Parametros da Terra - modelo axis simetrico (WGS-84)
Re=6378.1370e3; % m - Raio equatorial da Terra
we= 7.2921150e-5;  % (rad/s) - Velocidade de rotação da Terra com respeito ao espaco inercial
g=9.80665;   % m/s^2 - aceleracao da gravidade padrao ao nivel do mar
mut=3.986004418e14; % m3.s^-2
% Constantes de Jeffery
J2 = 0.00108263;
J3 = -0.00000254;
J4 = -0.00000161;

%% Condicoes iniciais - Centro espacial de Alcantara (do Google Maps)
h0=0;    % m - Altitude da base de lancamento
delta0 = -2.3267844*pi/180;   % rad - Latitude inicial
lon0 = -44.4111042*pi/180;  % rad - Longitude inicial
% Comprimento do trilho de lancamento
l_trilho=10;    % m

%% Dados fornecidos no teclado pelo usuario
T=input('Informe o tempo da simulacao (s): ');
V0=input('Informe o valor inicial da velocidade relativa (m/s): ');
A0=input('Informe o azimute inicial da velocidade relativa (graus): ');
A0=A0*pi/180;  % Azimute da velocidade relativa em t=0
phi0=input('Informe a elevacao inicial da velocidade relativa (graus): ');
phi0=phi0*pi/180; % Elevação da velocidade relativa em t=0

%% Parametros calculados
% Massa inicial do foguete
m0=sum(mp)+sum(ms)+mL;
% Distancia radial inicial
r0=Re+h0;

%% Estudo simplificado pela equacao de foguete
% Razão estrutural do primeiro e segundo estágios
sigma=ms./(ms+mp);
% Massa total na decolagem
m01=m0;
% Massa total na ignição do segundo estágio
m02=ms(2)+mp(2)+mL;
% Razão de carga útil do primeiro e segundo estágios
lamb(1)=m02/m01;
lamb(2)=mL/m02;
% Razão de carga útil total
lambL=prod(lamb);

% Velocidade de exaustão do primeiro e segundo estágios
ve=g*Isp;
% Delta v ideal da configuração original
Dv=-sum(ve.*log(sigma+(1-sigma).*lamb));

%% Mostra dados na tela
disp('Area de referencia do foguete com primeiro estagio (m^2):');disp(Sr(1));
disp('Area de referencia do foguete com segundo estagio (m^2):');disp(Sr(2));
disp('Area de referencia da carga util (m^2):');disp(Sr(3));
disp('Massa inicial antes da queima do primeiro estagio - kg');disp(m01);
disp('Massa inicial antes da queima do segundo estagio - kg');disp(m02);
disp('Massa da carga util - kg');disp(mL);
disp('Razoes estruturais');disp(sigma);
disp('Razoes de carga util');disp(lamb);
disp('Velocidades de exaustao - m/s');disp(ve);
disp('Razao de carga útil total');disp(lambL);
disp('Impulso de velocidade total ideal - m/s');disp(Dv);

%% Simulação
% Condição inicial
X0=[V0;A0;phi0;r0;delta0;lon0];
options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep ',0.5);
[t,X]=ode45('dinamica_foguete',[0 T],X0,options);

%% Pos processamento
N=length(t); % Numero de instantes de tempo
V=zeros(N,1);A=zeros(N,1);phi=zeros(N,1); % Velocidade - magnitude, azimute e elevacao
h=zeros(N,1);delta=zeros(N,1);lon=zeros(N,1); % Altitude, latitude e longitude PCPF
m=zeros(N,1);       % Massa
ft=zeros(N,1);      % Forca propulsiva
D=zeros(N,1);       % Forca de arrasto
q=zeros(N,1);   % Pressao dinamica
M=zeros(N,1);       % Numero de Mach
T=zeros(N,1);       % Temperatura
rho=zeros(N,1);       % Densidade
gc=zeros(N,1);       % Gravidade centripeta
gd=zeros(N,1);       % Gravidade latitudinal

for i=1:N
    V(i)=X(i,1);   % Velocidade relativa
    A(i)=X(i,2);   % Azimute da velocidade relativa
    phi(i)=X(i,3);   % Elevacao da velocidade relativa
    r=X(i,4);  % Distancia radial
    h(i)=r-Re; % Altitude
    delta(i)=X(i,5);    % Latitude
    lon(i)=X(i,6);   % Longitude no referencial fixo ao planeta
    [ft(i),m(i)]=propulsao_N_estagios(t(i));  % Forca propulsiva e massa
    [T(i),~,~,rho(i),~,M(i),~,~,Kn,~,~,R]=atm_padrao(h(i),V(i),lc,dT);  % Parametros atmosfericos
    [D(i),~,~]=aerodinamica_N_estagios(t(i),V(i),h(i),M(i),Kn,T(i),rho(i),R);    % Forcas aerodinamicas
    q(i)=0.5*rho(i)*V(i)^2; % Pressao dinamica
    [gc(i),gd(i)]=grav_axisimetrico(r,delta(i));  % Componentes da aceleracao da gravidade
end

%% Gráficos
figure(1)
subplot(231);plot(t,V,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('V (m/s)');
subplot(232);plot(t,A*180/pi,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('A (º)');
subplot(233);plot(t,phi*180/pi,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('\phi (º)');
subplot(234);plot(t,h/1e3,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('h (km)');
subplot(235);plot(t,delta*180/pi,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('\delta (º)');
subplot(236);plot(t,lon*180/pi,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('l (º)');
figure(2)
subplot(211);plot(t,ft,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('f_t (N)');
subplot(212);plot(t,m,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('m (kg)');
figure(3)
subplot(311);plot(t,D,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('D (N)');
subplot(323);plot(t,q,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('q (N/m^2)');
subplot(324);plot(t,M,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('M (-)');
subplot(325);plot(t,T-273.15,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('T (ºC)');
subplot(326);plot(t,rho,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('\rho (kg/m^3)');
figure(4)
subplot(121);plot(t,gc);grid;axis tight;xlabel('t (s)');ylabel('g_c (m/s^2)');
legend('grav. centripeta');
subplot(122);plot(t,gd);grid;axis tight;xlabel('t (s)');ylabel('g_\delta (m/s^2)');
legend('grav. latitudinal');