function voo_insercao_orbita
% Script para simular um voo de insercao em orbita
% Eh utilizado o foguete de sondagem VSB30 com os dados abaixo, retirados
% de referencias indicadas em aula. Ele eh
close all;clear all;clc;
global Re we mut tg J2 J3 J4 g lc dT Sr fc
global ms m0 mp  ti tq ts Isp h0 l_trilho

%% Parametros de massa estrutural e de carga util
ms(1)=284; % kg - Massa estrutural do primeiro estagio
ms(2)=320; % kg - Massa estrutural do segundo estagio
% Ajuste dos dois primeiros estagios
ds(1)=0.1*ms(1);   % Reducao de massa estrutural do primeiro estagio, convertida em massa de propelente
ds(2)=0.1*ms(2);   % Reducao de massa estrutural do segundo estagio, convertida em massa de propelente
ms(1)=ms(1)-ds(1);
ms(2)=ms(2)-ds(2);
% Massa estrutural no terceiro estagio
mo=400; % kg - massa original
ms(3)=55;   % Supor 30kg - valor muito baixo
% Supor 5 kg de carga útil (CubeSat 3U mais PPOD)
mL=5;

%% Parametros propulsivos
Isp(1)=261;     % s - Impulso específico do primeiro estágio (motor S31)
Isp(2)=261;     % s - Impulso específico do segundo estágio (motor S30)
mp(1)=677+ds(1);     % kg - Massa de propelente do primeiro estagio (motor S31). Original mais ajuste
mp(2)=898+ds(2);     % kg - Massa de propelente do segundo estagio (motor S30). Original mais ajuste
ti(1)=0; % s - Tempo da ignicao do primeiro estagio
ti(2)=15; % s - Tempo da ignicao do segundo estagio
tq(1)=13.5;    % s - Tempo do fim da queima do estágio 1 (motor S31)
tq(2)=44;       % s - Tempo do fim da queima do estágio 2 (motor S30)
ts(1)=13.5; % Tempo da separacao do primeiro estagio
ts(2)=59; % Tempo da separacao do segundo estagio
% Massa de propelente do novo estagio
mp(3)=mo-ms(3)-mL;    % Propelente eh tudo que sobra depois da massa de estrutura e de carga util
% propelente solido com alto impulso especifico
Isp(3)=350; % s
% Tempo da ignicao do terceiro estagio. Eh um dos parametros mais
% importantes para se obter a orbita desejada
ti(3)=ts(2)+185; % s - Determinado em relacao ao tempo de separacao do segundo estagio
% Tempo do fim da queima do terceiro estagio.
tq(3)=ti(3)+30; % s - Sem estudo profundo, soh tentativa e erro para gerar a resposta desejada
% Tempo de separacao do terceiro estagio. 
ts(3)=tq(2)+5;  % s - Sem estudo profundo. Um chute.
% Na mecanica de voo, este parametro soh tem importancia sobre o arrasto, o qual nao 
% eh importante na altitude orbital durante o lancamento. Na pratica, esta relacionado a questoes operacionais,
% tal como a comunicao e monitoramento do satelite, visto que o sistema de aquisicao de dados e
% telemetria do terceiro estagio eh usado para rastrear o satelite nos primeiros
% instantes da orbita

%% Parametros aerodinamicos e ambientais
% Fator de correcao do arrasto a partir de dados de tunel de vento
fc=1.28;
% As areas de referencia sao das secoes transversais
S1=pi*(0.557/2)^2; % m^2 - Area aproximada da secao transversal do primeiro estagio
S2=pi*(0.557/2)^2; % m^2 - Area aproximada da secao longitudinal do segundo estagio
S3=pi*(0.46/2)^2; % m^2 -  Area aproximada da secao longitudinal do terceiro estagio
Sp=pi*(0.46/2)^2; % m^2 -  Area aproximada da secao longitudinal da carga util
% Correcao para levar em conta a area molhada, feitas em funcao do
% comprimento. Assume-se que a correcao soh ocorre em 50% da area de
% referencia, supondo uma influencia de 50% da area molhada no arrasto
% total
lt=12.6; % m - Comprimento total
l2=lt-3.214;    %  Comprimento sem o primeiro estagio
% Comprimento da carga util
l4=0.5;    % m - um chute para um cubesat 3U, PPOD e acoplamentos
% Comprimento do terceiro estagio
l3=l2-3.294-l4; 
% Fatores de correcao da area molhada (para ajuste da area de referencia do
% arrasto)
f2=(l2/lt)*0.5+0.5; % Fator de correcao do segundo estagio
f3=(l3/lt)*0.5+0.5; % Fator de correcao do terceiro estagio
f4=(l4/lt)*0.5+0.5; % Fator de correcao da carga util
% Vetor de areas de referencia para calculo do arrasto
Sr=[S1
       S2*f2
       S3*f3
       Sp*f4];
lc=0.5; % Comprimento característico 
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
tg=0;   % s - Tempo em que o meridiano de referência tem longitude celeste nula

%% Condicoes iniciais - Centro espacial de Alcantara (do Google Maps)
h0=0;    % m - Altitude da base de lancamento
delta0 = -2.3267844*pi/180;   % rad - Latitude inicial
lon0 = -44.4111042*pi/180;  % rad - Longitude inicial
% Pequena correcao para compatibilizar com o tracado de mapas do MATLAB
lon0 = lon0-0.03*pi/180;
% Comprimento do trilho de lancamento
l_trilho=10;    % m

%% Parametros calculados
% Massa inicial do foguete
m0=sum(mp)+sum(ms)+mL;
% Distancia radial inicial
r0=Re+h0;

%% Estudo simplificado pela equacao de foguete
% Razoes estruturais
sigma=ms./(ms+mp);
% Massa total na decolagem
m01=m0;
% Massa total na ignição do segundo estágio
m02=m01-mp(1)-ms(1);
% Massa inicial antes da queima do terceiro estagio
m03=m02-mp(2)-ms(2);  % kg
% Razão de carga útil do primeiro e segundo estágios
lamb(1)=m02/m01;
lamb(2)=m03/m02;
% Razão de carga útil do terceiro estagio
lamb(3)=mL/m03;
% Razão de carga útil total
lambL=prod(lamb);
% Velocidades de exaustão
ve=g*Isp;
% Delta v total ideal
Dv=-sum(ve.*log(sigma+(1-sigma).*lamb));
% Delta v ideal impresso pelo terceiro estagio
Dv3=-ve(3)*log(sigma(3)+(1-sigma(3))*lamb(3));

%% Mostra dados na tela
disp('Area de referencia do foguete com primeiro estagio (m^2):');
disp(Sr(1));
disp('Area de referencia do foguete com segundo estagio (m^2):');
disp(Sr(2));
disp('Area de referencia do foguete com terceiro estagio (m^2):');
disp(Sr(3));
disp('Area de referencia da carga util (m^2):');
disp(Sr(4));
disp('Massa inicial antes da queima de cada estagio - kg');
disp(m0);
disp('Razoes estruturais');
disp(sigma);
disp('Razoes de carga util');
disp(lamb);
disp('Razao de carga útil total');
disp(lambL);
disp('Velocidades de exaustao - m/s');
disp(ve);
disp('Impulso de velocidade total ideal - m/s');
disp(Dv);
disp('Impulso de velocidade ideal impresso pelo terceiro estagio - m/s');
disp(Dv3);
%% Verifica se o usuario deseja continuar
cont=input('Deseja continuar ? (1 - sim) (0 - nao): ');

if cont==0
    return
end

%% Simulação
simula=1;
    while simula==1
    %% Dados fornecidos no teclado pelo usuario
    TF=input('Informe o tempo da simulacao (s): ');
    v0=input('Informe o valor inicial da velocidade relativa (m/s): ');
    phi0=input('Informe a condicao inicial do angulo de elevacao (graus): ');
    phi0=phi0*pi/180;
    in=input('Informe a inclinacao da orbita desejada (graus): ');
    in=in*pi/180;
    % Faz um teste de factibilidade usando a inclinacao da orbita e a latitude
    % inicial
    y=cos(in)/cos(delta0);
    if abs(y)>1
        disp('Nao eh possivel atingir a inclinacao a partir da latitude inicial. Calculando a menor possivel');
        y=sign(y);
    end
    Ai_f=asin(y); % Condicao final do azimute de velocidade inercial
    % Estimativa da condicao inicial de azimute de velocidade relativa
    % Raio de uma orbita circular de referencia com 100 km de altitude
    rref=Re+100e3;  
    % Velocidade de uma orbita de referencia com 100km de altitude
    viref=sqrt(mut/rref);
    A0=atan(tan(Ai_f)-(rref*we*cos(delta0))/(viref*cos(Ai_f)));
    disp('Condicao final de azimute de velocidade inercial (º): ');
    disp(Ai_f*180/pi);
    disp('Condicao inicial de azimute de velocidade relativa (º): ');
    disp(A0*180/pi);    

    % Condição inicial
    X0=[v0;A0;phi0;r0;delta0;lon0];
    %options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep ',1);
    options = odeset('RelTol',1e-8,'AbsTol',1e-12);
    %[t,X]=ode45('dinamica_foguete',[0 TF],X0,options);
    [t,X]=ode15s('dinamica_foguete',[0 TF],X0,options);
    
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
        V(i)=X(i,1);   % Velocidade relativa
        A(i)=X(i,2);   % Azimute da velocidade relativa
        phi(i)=X(i,3);   % Elevacao da velocidade relativa
        h(i)=X(i,4)-Re; % Altitude
        r=X(i,4);  % Distancia radial
        delta(i)=X(i,5);    % Latitude
        lon(i)=X(i,6);   % Longitude no referencial fixo ao planeta
        [ft(i),m(i)]=propulsao_N_estagios(t(i));  % Forca propulsiva e massa
        [T(i),~,rho(i),~,M(i),~,~,Kn,~,~]=atm_padrao(h(i),V(i),lc,dT);  % Parametros atmosfericos
        [D(i),~,~]=aerodinamica_N_estagios(t(i),V(i),h(i),M(i),Kn,T(i),rho(i));    % Forcas aerodinamicas
        q(i)=0.5*rho(i)*V(i)^2; % Pressao dinamica
        % Coordenadas da velocidade inercial no referencial horizontal local
        [Vi(i),phii(i),Ai(i)]=Vrel2Vine(V(i),phi(i),A(i),we,r,delta(i));
        % Longitude celeste
         longc(i)=long_ECEF2ECI(t(i),lon(i),we,tg);
        % Energia especifica da orbita
        ee(i)=Vi(i)^2/2-mut/r;
        % Posicao e velocidade inercial no referencial celeste ECI
        [rc0,vc0]=RvelPolar2RvelRet(Vi(i),Ai(i),phii(i),r,delta(i),longc(i));
        R0(i,:)=rc0';
        % Parametros orbitais
        par_orb=det_orbita(t(i),rc0,vc0,mut);
        a(i)=par_orb(1);e(i)=par_orb(2);tau(i)=par_orb(3);OM(i)=par_orb(4);in(i)=par_orb(5);om(i)=par_orb(6);
    end
        %% Analise de orbita
        % Geracao de vetores para tracar grafico
        % Semi eixo maior de uma orbita circular de referencia com 100 km de altitude
        ar=rref*ones(1,N);  
        % Velocidade de uma orbita de referencia com 100km de altitude
        Vir=viref*ones(1,N);  
        % Energia especifica de uma orbita de referencia com a=100km
        eer=-mut./(2*ar); 
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
        disp('Raio do periastro (km): ');
        disp(rp/1e3);
        disp('Raio do apoastro (km): ');
        disp(ra/1e3);
        disp('Altitude do periastro (km): ');
        disp((rp-Re)/1e3);
        disp('Altitude do apoastro (km): ');
        disp((ra-Re)/1e3);

        %% Gráficos
        figure(1)
        subplot(231)
        plot(t,V,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('V (m/s)');
        subplot(232)
        plot(t,A*180/pi,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('A (º)');
        subplot(233)
        plot(t,phi*180/pi,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('\phi (º)');
        subplot(234)
        plot(t,h/1e3,'LineWidth',2);hold
        plot(t,hfq/1e3,'--',tfq,hfq(1)/1e3,'*');
        grid;axis tight;xlabel('t (s)');ylabel('h (km)');
        legend('altitude','altitude no fim da queima do 3º estagio');
        subplot(235)
        plot(t,delta*180/pi,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('\delta (º)');
        subplot(236)
        plot(t,lon*180/pi,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('l(º)');

        figure(2)
        subplot(221);
        plot(t,Vi,'LineWidth',2);hold;
        plot(t,Vir,'--',t,Vfq,'-.',tfq,Vfq(1),'*');grid;xlabel('t (s)');ylabel('V_i (m/s)');
        legend('Velocidade inercial','Velocidade de orbita circular de referencia de 100km de altitude',...
            'Velocidade no fim da queima do terceiro estagio');
        subplot(222);
        plot(t,Ai*180/pi,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('A_i (º)');
        subplot(223);
        plot(t,phii*180/pi,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('\phi_i (º)');
        subplot(224);
        plot(t,longc*180/pi,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('\lambda (º)');

        figure(3)
        subplot(211)
        plot(t,ft,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('f_t (N)');
        subplot(212)
        plot(t,m,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('m (kg)');

        figure(4)
        subplot(311)
        plot(t,D,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('D (N)');
        subplot(323)
        plot(t,q,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('q (N/m^2)');
        subplot(324)
        plot(t,M,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('M (-)');
        subplot(325)
        plot(t,T-273.15,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('T (ºC)');
        subplot(326)
        plot(t,rho,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('\rho (kg/m^3)');

        figure(5)
        subplot(331);
        plot(t,a/1e3,'LineWidth',2);hold;
        legend('Semi eixo maior');
        plot(t,ar/1e3,'--',t,Re*ones(1,N)/1e3,'-.');
        grid;xlabel('t (s)');ylabel('a (km)');
        legend('Semi eixo maior de orbita com 100km','Raio da Terra');
        subplot(332);
        plot(t,e,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('e (-)');
        subplot(333);
        plot(t,tau,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('\tau (s)');
        subplot(334);
        plot(t,OM*180/pi,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('\Omega (º)');
        subplot(335);
        plot(t,in*180/pi,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('i (º)');
        subplot(336);
        plot(t,om*180/pi,'LineWidth',2);grid;axis tight;xlabel('t (s)');ylabel('\omega (º)');
        subplot(325);
        plot(t,ee,t,eer,'--','LineWidth',2);grid;xlabel('t (s)');ylabel('\epsilon (J/kg)');
        legend('Energia especifica','Energia especifica de orbita com 100km');

        figure(6)
        traj=[delta lon]*180/pi;
        desenha_mapa_trajetoria([delta0*180/pi,lon0*180/pi,h0],traj)

        figure(7)
        plot3(R0(:,1)/1e3,R0(:,2)/1e3,R0(:,3)/1e3,'LineWidth',2);
        xlabel('X (km)');ylabel('Y (km)');zlabel('Z (km)');axis tight
        hold on
        h1=gca;
        earth_sphere(h1,'km') 

        simula=input('Deseja simular novamente ? (1) sim (0) nao: ');
    end
end