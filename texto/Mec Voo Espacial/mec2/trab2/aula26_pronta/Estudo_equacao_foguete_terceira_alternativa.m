%% Aplicação prática de cálculo de impulso de velocidade total
% Dados apresentados na referência
% Pedro L. K. da Cás, Carlos A. G. Veras, Olexiy Shynkarenko, and Rodrigo Leonardi. A
% brazilian space launch system for the small satellite market. Aerospace, 6(11), 2019
% Veiculo  C-2  da referencia alterado no terceiro estagio. Chutar
% massas de propelente. Configuracao com 3, 4 e 5 motores de primeiro estagio.

%       1° estágio| 2º estágio | 3° estágio
% C-2-c:  3xS50     | 1xS50      | 1xRD843
% C-2-d:  4xS50     | 1xS50      | 1xRD843
% C-2-e:  5xS50     | 1xS50      | 1xRD843

clc;
clear all
%% Constantes
g=9.81; % m/s^2
%% Dados
% Massa de carga útil - Supor um CubeSat 8U (8*1,33+2,36 kg)
mL=13; % kg
% Variacoes da massa de propelente
fp=0.1:0.01:2;  % De metade ate o dobro
N=length(fp);
Dv_c2c=zeros(1,N);Dv_c2d=zeros(1,N);Dv_c2e=zeros(1,N);    % Impulsos resultantes
ms3=zeros(1,N); % Massa estrutural do 3° estagio
for i=1:N
% Massa de propelente
mp_c2c=[33157 11058 fp(i)*609];   % kg - C-2-c
mp_c2d=[33157*4/3 11058 fp(i)*609];   % kg - C-2-d
mp_c2e=[33157*5/3 11058 fp(i)*609];   % kg - C-2-e
% Massa estutural
ms3(i)=(0.21/(1-0.21))*mp_c2c(3);
ms_c2c=[4650 1367 ms3(i)];   % kg - C-2-c
ms_c2d=[4650*4/3 1367 ms3(i)];   % kg - C-2-d
ms_c2e=[4650*5/3 1367 ms3(i)];   % kg - C-2-e
% Impulso específico ideal do primeiro e segundo estágio
Isp_c2=[251 271 315]; % s - C-2
%% Cálculos
% Razões estruturais
sigma_c2c=ms_c2c./(ms_c2c+mp_c2c);
sigma_c2d=ms_c2d./(ms_c2d+mp_c2d);
sigma_c2e=ms_c2e./(ms_c2e+mp_c2e);
% Massa total no inicio da queima de cada estagio
% C-2-c
m0_c2c(1)=sum(ms_c2c)+sum(mp_c2c)+mL;
m0_c2c(2)=m0_c2c(1)-ms_c2c(1)-mp_c2c(1);
m0_c2c(3)=m0_c2c(2)-ms_c2c(2)-mp_c2c(2);
% C-2-d
m0_c2d(1)=sum(ms_c2d)+sum(mp_c2d)+mL;
m0_c2d(2)=m0_c2d(1)-ms_c2d(1)-mp_c2d(1);
m0_c2d(3)=m0_c2d(2)-ms_c2d(2)-mp_c2d(2);
% C-2-e
m0_c2e(1)=sum(ms_c2e)+sum(mp_c2e)+mL;
m0_c2e(2)=m0_c2e(1)-ms_c2e(1)-mp_c2e(1);
m0_c2e(3)=m0_c2e(2)-ms_c2e(2)-mp_c2e(2);

% Razões de carga útil
% C-2-c
lamb_c2c(1)=m0_c2c(2)/m0_c2c(1);
lamb_c2c(2)=m0_c2c(3)/m0_c2c(2);
lamb_c2c(3)=mL/m0_c2c(3);
% C-2-d
lamb_c2d(1)=m0_c2d(2)/m0_c2d(1);
lamb_c2d(2)=m0_c2d(3)/m0_c2d(2);
lamb_c2d(3)=mL/m0_c2d(3);
% C-2-e
lamb_c2e(1)=m0_c2e(2)/m0_c2e(1);
lamb_c2e(2)=m0_c2e(3)/m0_c2e(2);
lamb_c2e(3)=mL/m0_c2e(3);

% Razão de carga útil total
% C-2-c
lambL_c2c=prod(lamb_c2c);
% C-2-d
lambL_c2d=prod(lamb_c2d);
% C-2-e
lambL_c2e=prod(lamb_c2e);

% Impulso de velocidade
% C-2-c
Dv_c2c(i)=-sum(g*Isp_c2.*log(sigma_c2c+(1-sigma_c2c).*lamb_c2c));
% C-2-d
Dv_c2d(i)=-sum(g*Isp_c2.*log(sigma_c2d+(1-sigma_c2d).*lamb_c2d));
% C-2-e
Dv_c2e(i)=-sum(g*Isp_c2.*log(sigma_c2e+(1-sigma_c2e).*lamb_c2e));
end
figure;
subplot(141);plot(fp*100,Dv_c2c/1e3);xlabel('Fração de m_p original - %');ylabel('\Delta v - km/s');
grid;legend('Primeiro estágio com 3 motores S50');
subplot(142);plot(fp*100,Dv_c2d/1e3);xlabel('Fração de m_p original - %');ylabel('\Delta v - km/s');
grid;legend('Primeiro estágio com 4 motores S50');
subplot(143);plot(fp*100,Dv_c2e/1e3);xlabel('Fração de m_p original - %');ylabel('\Delta v - km/s');
grid;legend('Primeiro estágio com 5 motores S50');
subplot(144);plot(fp*100,ms3);xlabel('Fração de m_p original - %');ylabel('m_{s_3} - kg');grid;