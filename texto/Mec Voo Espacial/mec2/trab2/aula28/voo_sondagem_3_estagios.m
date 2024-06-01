% Script para simular um voo de sondagem de veiculo de dois estagios com
% carga util.
% Eh utilizado o foguete de insercao orbital geosincrono conceitual
% (FIOGSC)
% de referencias indicadas em aula.
close all;clear all;clc;
global Re we mut J2 J3 J4 g lc dT Sr fc
global ms m0 mp  ti tq ts Isp h0 l_trilho tg
%% Parametros propulsivos
Isp(1)=251;     % s - Impulso específico do primeiro estágio (5xS50)
Isp(2)=271;     % s - Impulso específico do segundo estágio (1xS50)
Isp(3)=315;     % s - Impulso específico do segundo estágio (RD843)
mp(1)=5.5262e+04;     % kg - Massa de propelente do primeiro estagio (5xS50)
mp(2)=11058;     % kg - Massa de propelente do segundo estagio (1xS50)
mp(3)=  243.6;     % kg - Massa de propelente do segundo estagio (RD843)
Tq1=62;
Tq2=64.62;
Tq3=301;
ti(1)=0; % s - Tempo da ignicao do primeiro estagio
tq(1)=ti(1)+Tq1;    % s - Tempo do fim da queima do estágio 1
ts(1)=tq(1)+1; % Tempo da separacao do primeiro estagio
ti(2)=ts(1)+1; % s - Tempo da ignicao do segundo estagio
tq(2)=ti(2)+Tq2;       % s - Tempo do fim da queima do estágio 2
ts(2)=tq(2)+1; % Tempo da separacao do segundo estagio
ti(3)=ts(2)+1; % s - Tempo da ignicao do segundo estagio
tq(3)=ti(3)+Tq3;       % s - Tempo do fim da queima do estágio 2
ts(3)=tq(3)+1; % Tempo da separacao do segundo estagio
%% Parametros de massa estrutural e de carga util
ms(1)=7750; % kg - Massa estrutural do primeiro estagio
ms(2)=1367; % kg - Massa estrutural do segundo estagio
ms(3)=   64.7544; % kg - Massa estrutural do segundo estagio
mL=13; % kg - Massa da carga util

%% Parametros aerodinamicos e ambientais
% Fator de correcao do arrasto a partir de dados de tunel de vento
fc=1.28;
% As areas de referencia sao das secoes transversais
S1=4.6*5/3; % m^2 - Area aproximada da secao transversal do primeiro estagio
S2=1.5; % m^2 - Area aproximada da secao longitudinal do segundo estagio
S3=1.5; % m^2 - Area aproximada da secao longitudinal do terceiro estagio
SL=1.5; % m^2 -  Area aproximada da secao longitudinal da carga util
% Correcoes para levar em conta a area molhada, feitas em funcao do
% comprimento. Assume-se que a correcao soh ocorre em 50% da area de
% referencia, supondo uma influencia de 50% da area molhada no arrasto
% total
lt=7.33+7.1+6.28; % m - Comprimento total
l2=7.1+6.28;    %  Comprimento sem o primeiro estagio
l3=6.28;    %  Comprimento sem o segundo estagio
l4=1; % Comprimento da carga util
f2=(l2/lt)*0.5+0.5; % Fator de correcao do segundo estagio
f3=(l3/lt)*0.5+0.5; % Fator de correcao do terceiro estagio
f4=(l4/lt)*0.5+0.5; % Fator de correcao da carga util
% Vetor de areas de referencia para calculo do arrasto
Sr=[S1
       S2*f2
       S3*f3
       SL*f4];
lc=1.5; % Comprimento característico - diametro dos estagios 2 e superiores
dT=10;  % K - Delta T em relação à atmosfera padrão (que é 15ºC no nível do mar)

%% Parametros da Terra - modelo axis simetrico (WGS-84)
Re=6378.1370e3; % m - Raio equatorial da Terra
we= 7.2921150e-5;  % (rad/s) - Velocidade inercial de rotação da Terra
g=9.80665;   % m/s^2 - aceleracao da gravidade padrao ao nivel do mar
mut=3.986004418e14; % m3.s^-2
% Constantes de Jeffery
J2 = 0.00108263;
J3 = -0.00000254;
J4 = -0.00000161;
tg=0;   % s - Tempo em que o meridiano de referência tem longitude celeste nula

%% Condicoes iniciais - Centro espacial de Alcantara (do Google Maps)
h0=0;    % m - Altitude da base de lancamento
delta0 = -2.3267844*pi/180;   % rad - Latitude inicial
lon0 = -44.4111042*pi/180;  % rad - Longitude inicial
% Comprimento do trilho de lancamento
l_trilho=lt;    % m - igual ao comprimento total do foguete

%% Dados fornecidos no teclado pelo usuario
T=600; %input('Informe o tempo da simulacao (s): ');
V0=5;%input('Informe o valor inicial da velocidade relativa (m/s): ');
A0=85;%input('Informe o azimute inicial da velocidade relativa (graus): ');
A0=A0*pi/180;  % Azimute da velocidade relativa em t=0
phi0=85;%input('Informe a elevacao inicial da velocidade relativa (graus): ');
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
m02=ms(2)+mp(2)+ms(3)+mp(3)+mL;
% Massa total na ignição do terceiro estágio
m03=ms(3)+mp(3)+mL;
% Razão de carga útil do primeiro, segundo e terceiro estágios
lamb(1)=m02/m01;
lamb(2)=m03/m02;
lamb(3)=mL/m03;
% Razão de carga útil total
lambL=prod(lamb);

% Velocidade de exaustão
ve=g*Isp;
% Delta v ideal da configuração original
Dv=-sum(ve.*log(sigma+(1-sigma).*lamb));

%% Mostra dados na tela
disp('Area de referencia do foguete com primeiro estagio (m^2):');disp(Sr(1));
disp('Area de referencia do foguete com segundo estagio (m^2):');disp(Sr(2));
disp('Area de referencia do foguete com terceiro estagio (m^2):');disp(Sr(3));
disp('Area de referencia da carga util (m^2):');disp(Sr(4));
disp('Massa inicial antes da queima do primeiro estagio - kg');disp(m01);
disp('Massa inicial antes da queima do segundo estagio - kg');disp(m02);
disp('Massa inicial antes da queima do terceiro estagio - kg');disp(m03);
disp('Massa da carga util - kg');disp(mL);
disp('Razoes estruturais');disp(sigma);
disp('Razoes de carga util');disp(lamb);
disp('Velocidades de exaustao - m/s');disp(ve);
disp('Razao de carga útil total');disp(lambL);
disp('Impulso de velocidade total ideal - m/s');disp(Dv);

%% Simulação
% Condição inicial
X0=[V0;A0;phi0;r0;delta0;lon0];
options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.05);
[t,X]=ode45('dinamica_foguete',[0 T],X0,options);

%% Pos processamento
N=length(t); % Numero de instantes de tempo
% Magnitude, azimute e elevacao da velocidade relativa
V=zeros(N,1);A=zeros(N,1);phi=zeros(N,1); 
% Altitude, latitude e longitude no referencial fixo ao planeta
h=zeros(N,1);delta=zeros(N,1);lon=zeros(N,1); 
m=zeros(N,1);       % Massa
ft=zeros(N,1);      % Forca propulsiva
D=zeros(N,1);       % Forca de arrasto
q=zeros(N,1);   % Pressao dinamica
M=zeros(N,1);       % Numero de Mach
T=zeros(N,1);       % Temperatura
rho=zeros(N,1);       % Densidade
Vi=zeros(N,1);  % Magnitude da velocidade inercial
phii=zeros(N,1);    % Elevacao da velocidade inercial
Ai=zeros(N,1);  % Azimute da velocidade inercial
longc=zeros(N,1);   % Longitude celeste
ee=zeros(N,1);  % Energia especifica
a=zeros(N,1);   % Semi eixo maior da orbita
e=zeros(N,1);   % Excentricidade da orbita
tau=zeros(N,1);  % Tempo de periastro
OM=zeros(N,1);  % Ascencao reta do nodo ascendente
in=zeros(N,1);  % Inclinacao da orbita
om=zeros(N,1);  % Argumento de periastro
R0=zeros(N,3);  % Posicao no referencial ECI
for i=1:N
    % Magnitude, azimute e elevacao da velocidade relativa
    V(i)=X(i,1);A(i)=X(i,2);phi(i)=X(i,3);
    % Posicao no referencial PCPF
    h(i)=X(i,4)-Re;r=X(i,4);delta(i)=X(i,5);lon(i)=X(i,6);
    [ft(i),m(i)]=propulsao_N_estagios(t(i));  % Forca propulsiva e massa
    % Parametros atmosfericos
    [T(i),~,~,rho(i),~,M(i),~,~,Kn,~,~,R]=atm_padrao(h(i),V(i),lc,dT);  
    % Forcas aerodinamicas
    [D(i),~,~]=aerodinamica_N_estagios(t(i),V(i),h(i),M(i),Kn,T(i),rho(i),R);    
    q(i)=0.5*rho(i)*V(i)^2; % Pressao dinamica
    % Coordenadas da velocidade inercial no referencial LVLH
    [Vi(i),phii(i),Ai(i)]=Vrel2Vine(V(i),phi(i),A(i),we,r,delta(i));
    % Longitude celeste
     longc(i)=long_ECEF2ECI(t(i),lon(i),we,tg);
    % Energia especifica da orbita
    ee(i)=Vi(i)^2/2-mut/r;
    % Posicao e velocidade inercial no referencial ICP
    [rc0,vc0]=RvelPolar2RvelRet(Vi(i),Ai(i),phii(i),r,delta(i),longc(i));R0(i,:)=rc0';
    % Elemetnos orbitais
    par_orb=det_orbita(t(i),rc0,vc0,mut);
    a(i)=par_orb(1);e(i)=par_orb(2);tau(i)=par_orb(3);
    OM(i)=par_orb(4);in(i)=par_orb(5);om(i)=par_orb(6);
end
        %% Analise de orbita   
        n_gso=we; % Velocidade angular media da orbita geosincrona
        a_gso=(mut/n_gso^2)^(1/3); % Semi eixo maior da orbita geosincrona
        viref=sqrt(mut/a_gso); % Velocidade de uma orbita de referencia circular GSO
        % Geracao de vetores para tracar grafico       
        ar=a_gso*ones(1,N);  % Semi eixo maior de uma orbita circular de referencia GSO
        Vir=viref*ones(1,N);   % Velocidade de uma orbita de referencia GSO
        eer=-mut./(2*ar); % Energia especifica de uma orbita de referencia GSO
        % Altitude e velocidade inercial no fim da queima do terceiro estagio
        for i=1:N
            if t(i)>tq(3)
                break
            end
        end
        ifq=i-1;
        tfq=t(ifq); % Tempo do fim da queima do terceiro estagio
        Vfq=Vi(ifq)*ones(1,N);  % Velocidade inercial no fim da queima do terceiro estagio
        hfq=h(ifq)*ones(1,N);   % Altitude no fim da queima do terceiro estagio
        P=2*pi*sqrt((Re+hfq(1))^3/mut);    % Periodo da orbita obtida
        disp('Periodo da orbita obtida (min): ');
        disp(P/60);
        rp=a(ifq)*(1-e(ifq));   % Raio do periastro
        ra=a(ifq)*(1+e(ifq));   % Raio do apoastro
        disp('Raio do periastro (km): ');disp(rp/1e3);
        disp('Raio do apoastro (km): ');disp(ra/1e3);
        disp('Altitude do periastro (km): ');disp((rp-Re)/1e3);
        disp('Altitude do apoastro (km): ');disp((ra-Re)/1e3);

        %% Gráficos
        figure(1)
        subplot(231);plot(t,V,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('V (m/s)');
        subplot(232);plot(t,A*180/pi,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('A (º)');
        subplot(233);plot(t,phi*180/pi,'LineWidth',2);grid;axis tight;
        xlabel('t (s)');ylabel('\phi (º)');
        subplot(234);plot(t,h/1e3,'LineWidth',2);hold;plot(t,hfq/1e3,'--',tfq,hfq(1)/1e3,'*');
        grid;axis tight;xlabel('t (s)');ylabel('h (km)');
        legend('altitude','altitude no fim da queima do 3º estagio');
        subplot(235);plot(t,delta*180/pi,'LineWidth',2);grid;axis tight;
        xlabel('t (s)');ylabel('\delta (º)');
        subplot(236);plot(t,lon*180/pi,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('l(º)');

        figure(2)
        subplot(221);plot(t,Vi,'LineWidth',2);hold;plot(t,Vir,'--',t,Vfq,'-.',tfq,Vfq(1),'*');
        grid;xlabel('t (s)');ylabel('V_i (m/s)');
        legend('Velocidade inercial','Velocidade de orbita circular de referencia GSO',...
            'Velocidade no fim da queima do terceiro estagio');
        subplot(222);plot(t,Ai*180/pi,'LineWidth',2);grid;axis tight;
        xlabel('t (s)');ylabel('A_i (º)');
        subplot(223);plot(t,phii*180/pi,'LineWidth',2);grid;axis tight;
        xlabel('t (s)');ylabel('\phi_i (º)');
        subplot(224);plot(t,longc*180/pi,'LineWidth',2);grid;axis tight;
        xlabel('t (s)');ylabel('\lambda (º)');

        figure(3)
        subplot(211);plot(t,ft,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('f_t (N)');
        subplot(212);plot(t,m,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('m (kg)');

        figure(4)
        subplot(311);plot(t,D,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('D (N)');
        subplot(323);plot(t,q,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('q (N/m^2)');
        subplot(324);plot(t,M,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('M (-)');
        subplot(325);plot(t,T-273.15,'LineWidth',2);grid;axis tight;
        xlabel('t (s)');ylabel('T (ºC)');
        subplot(326);plot(t,rho,'LineWidth',2);grid;axis tight; 
        xlabel('t (s)');ylabel('\rho (kg/m^3)');

        figure(5)
        subplot(311);plot(t,ee,t,eer,'--','LineWidth',2);grid;
        xlabel('t (s)');ylabel('\epsilon (J/kg)');
        legend('Energia especifica','Energia especifica de orbita GSO');
        subplot(334);plot(t,a/1e3,'LineWidth',2);hold;legend('Semi eixo maior');
        plot(t,ar/1e3,'--',t,Re*ones(1,N)/1e3,'-.');grid;xlabel('t (s)');ylabel('a (km)');
        legend('Semi eixo maior','Semi eixo maior de orbita GSO','Raio da Terra');
        subplot(335);plot(t,e,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('e (-)');
        subplot(336);plot(t,tau,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('\tau (s)');
        subplot(337);plot(t,OM*180/pi,'LineWidth',2);grid;axis tight;
        xlabel('t (s)');ylabel('\Omega (º)');
        subplot(338);plot(t,in*180/pi,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('i (º)');
        subplot(339);plot(t,om*180/pi,'LineWidth',2);grid;axis tight;
        xlabel('t (s)');ylabel('\omega (º)');