function [T,p,rho,ainf,M,mu,Pr,Kn,d,Re]=atm_padrao(h,v,lc,dT)
%% Funcao para calculo da atmosfera padrao para a altitude geometrica de zero ate 2000 km
% Eh utilizado o modelo apresentado na referencia 
% TEWARI, A. Atmospheric and Space Flight Dynamics: 
% Modelling and simulation with MATLAB and Simulink. Boston: Birkhauser, 2007.
% Capitulo 9 - Planetary atmosphere
% O modelo possui 21 camadas, sendo uma combinacao dos modelos de atmosfera
% padrao norte americanos de 1962 e 1976.
% Todas as camadas consideradas possuem variacao linear de temperatura ou
% sao isotermicas. Na faixa de altitude de zero ate 86km, eh utilizado o
% perfil de temperaturas da atmosfera padrao norte americana do ano de
% 1976. Para a faixa de altitude de 86km a 2000km, eh utilizado o modelo de
% 1962.
% Apesar de a atmosfera nao apresentar equilibrio termico e quimico acima
% de 86km, sendo os perfis de temperatura, nessas regioes, nao lineares,
% esta funcao adota do modelo linear de temperaturas de 1962, nessa faixa,
% como uma aproximacao.
% Entradas:
% h - altitude geometrica [m]
% v - velocidade do veiculo em relacao ao escoamento nao perturbado, em
% metros por segundo [m/s]
% lc - comprimento caracteristico do veiculo, em metros [m]
% dT - Variacao de temperatura em relacao a temperatura padrao [K ou °C]
% Saidas
% T - Temperatura, em Kelvin [K]
% p - Pressao, em Pascal [N/m^2]
% rho - Densidade [kg/m^2]
% ainf - Velocidade do som [m/s]
% M - Numero de Mach [adm]
% mu - Coeficiente de viscosidade dinamica [kg/m*s] 
% Pr - Numero de Prandtl [adm]
% Kn - Numero de Knudsen [adm]
% d - Parametro de regime de escoamento [adm]
% Re - Numero de Reynolds [adm]
%% Entrada de dados
% Vetores com os seguintes dados, cada linha eh uma camada do modelo de 
% atmosfera padrao: altitude no inicio da camada (m), temperatura no inicio
% da camada (K), constante de gas ideal do ar na camada (J/kg.K), taxa de
% lapso termico (K/m)
hi=1e3*[0
            11.0191
            20.0631
            32.1619
            47.3501
            51.4125
            71.8020
            86
            100
            110
            120
            150
            160
            170
            190
            230
            300
            400
            500
            600
            700];

Ti=[288.15
    216.65
    216.65
    228.65
    270.65
    270.65
    214.65
    186.946
    210.02
    257.0
    349.49
    892.79
    1022.2
    1103.4
    1205.4
    1322.3
    1432.1
    1487.4
    1506.1
    1506.1
    1507.6];

R=[287.0
    287.0
    287.0
    287.0
    287.0
    287.0
    287.02
    287.02
    287.84
    291.06
    308.79
    311.80
    313.69
    321.57
    336.68
    366.84
    416.88
    463.36
    493.63
    514.08
    514.08];

a=1e-3*[-6.5
            0.0
            1.0
            2.8
            0.0
            -2.8
            -2.0
            1.693
            5.0
            10.0
            20.0
            15.0
            10.0
            7.0
            5.0
            4.0
            3.3
            2.6
            1.7
            1.1
            0.0];
% Constantes
g0 = 9.80665; % Valor ao nivel do mar da aceleracao da gravidade (m/s^2)
Na = 6.0220978e23; % Numero de Avogadro
sigma = 3.65e-10; % Diametro de colisao para o ar (m)
m0 = 28.964e-3; % Massa molar do ar ao nivel do mar (kg/Mol)
P0 = 1.01325e5; % Pressao padrao ao nivel do mar (N/m^2)
Re = 6378.14e3; % Raio medio da Terra (m)
gamma = 1.405; % Razao de calores especificos ao nivel do mar
%% Constantes calculadas
% Numero beta associado a distancia radial media do nivel do mar
beta = 2/Re;
%% Identifica a camada a qual a altitude pertence
% Contador 
i=1;
if h<0
    %disp('Cuidado! Altitude negativa.');
    %disp('Os resultados apresentados dizem respeito a h=0.');
    i=1;h=0;
elseif h>2000e3
    %disp('A altitude fornecida esta acima do limite superior de 2.000 km.');
    %disp('Os resultados apresentados dizem respeito a h=2.000 km');
    i=21;h=2000e3;
else
    for i=1:21
        if (i==21)
            break
        elseif ((h>=hi(i))&&(h<hi(i+1)))
            break
        end
    end
end
%% Realiza os calculos
% Pressao no inicio da camada i
Pi=P0; % Pressao inicial da camada - inicializa com o valor ao nivel do mar
px=P0;  % Variavel auxiliar para guardar o valor de pressao inicial na camada anterior
for j=2:i
    if a(j-1)~=0 % Verifica se a camada nao eh isotermica
        % Calcula a pressao inicial da camada i a partir do modelo da
        % camada i-1 
        A=1+(a(j-1)*(hi(j)-hi(j-1)))/Ti(j-1);
        B=-(g0/(R(j-1)*a(j-1)))*(1+beta*(Ti(j-1)/a(j-1)-hi(j-1)));
        C=(g0*beta/(R(j-1)*a(j-1)))*(hi(j)-hi(j-1));
        Pi=px*(A^B)*exp(C);
        px=Pi;  % O valor atual sera o valor anterior na proxima iteracao
    else
        % Calcula a pressao inicial da camada i pelo modelo isotermico
        Pi=px*exp(-(g0/(R(j-1)*Ti(j-1)))*(hi(j)-hi(j-1))*(1-(beta/2)*(hi(j)-hi(j-1))));
        px=Pi;  % O valor atual sera o valor anterior na proxima iteracao
    end
end
% Temperatura padrao
T=Ti(i)+a(i)*(h-hi(i));
% Corrige pelo valor de delta T adotado
T=T+dT;
% Pressao
if a(i)~=0 % Verifica se a camada nao eh isotermica
        % Calcula a pressao
        A=1+(a(i)*(h-hi(i)))/Ti(i);
        B=-(g0/(R(i)*a(i)))*(1+beta*(Ti(i)/a(i)-hi(i)));
        C=(g0*beta/(R(i)*a(i)))*(h-hi(i));
        p=Pi*(A^B)*exp(C);
else
        % Calcula a pressao pelo modelo isotermico
        p=Pi*exp(-(g0/(R(i)*Ti(i)))*(h-hi(i))*(1-(beta/2)*(h-hi(i))));
end
% Calcula a densidade
rho=p/(R(i)*T);
% Velocidade do som
ainf=sqrt(gamma*R(i)*T);
% Numero de Mach
M=v/ainf;
% Coeficiente de viscosidade dinamica
mu=1.458e-6*(T^(3/2))/(T+110.4);
% Numero de Prandtl
cp=R(i)*gamma/(gamma-1);
kT=(2.64638e-3*(T^(3/2)))/(T+245.4*(10^(-12/T)));
Pr=mu*cp/kT;
% Numero de Knudsen
lam=m0/(sqrt(2)*pi*sigma^2*rho*Na);
Kn=lam/lc;
% Parametro de regime de escoamento
if Kn>=10
    d=1;
elseif Kn<=0.01
    d=2;
else
    d=3;
end
% Re - Numero de Reynolds [adm]
Re=rho*v*lc/mu;
end