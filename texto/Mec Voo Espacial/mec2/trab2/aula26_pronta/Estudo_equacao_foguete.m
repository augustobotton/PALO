%% Aplicação prática de cálculo de impulso de velocidade total
% Dados das 4 configurações de foguetes apresentados na referência
% Pedro L. K. da Cás, Carlos A. G. Veras, Olexiy Shynkarenko, and Rodrigo Leonardi. A
% brazilian space launch system for the small satellite market. Aerospace, 6(11), 2019
% VLM-1 e 3 veículos feitos com combinações de seus motores
%       1° estágio| 2º estágio | 3° estágio
% VLM-1 1xS50     | 1xS50      | 1xS44
% C-1:  3xS50     | 1xS50      | 1xS44
% C-2:  3xS50     | 1xS50      | 1xRD843
% C-3:  3xS50     | 1xS50      | 4xRD843
clc;
clear all
%% Constantes
g=9.81; % m/s^2
%% Dados
% Massa de carga útil - Supor um CubeSat 8U (8*1,33+2,36 kg)
mL=13; % kg
% Massa estutural
ms_vlm=[1367 1367 166.5];   % kg - VLM-1
ms_c1=[4650*6/3 1367 166.5];   % kg - C-1
ms_c2=[4650*6/3 1367 161.9];   % kg - C-2
ms_c3=[4650*6/3 1367 228.7];   % kg - C-3
% Massa de propelente
mp_vlm=[11058 11058 813];   % kg - VLM-1
mp_c1=[33157*6/3 11058 813];   % kg - C-1
mp_c2=[33157*6/3 11058 609];   % kg - C-2
mp_c3=[33157*6/3 11058 811];   % kg - C-3
% Impulso específico ideal do primeiro e segundo estágio
Isp_vlm=[271 271 270]; % s - VLM-1
Isp_c1=[251 271 270]; % s - C-1
Isp_c2=[251 271 315]; % s - C-2
Isp_c3=[251 271 315]; % s - C-3
%% Cálculos
% Razões estruturais
sigma_vlm=ms_vlm./(ms_vlm+mp_vlm);
sigma_c1=ms_c1./(ms_c1+mp_c1);
sigma_c2=ms_c2./(ms_c2+mp_c2);
sigma_c3=ms_c3./(ms_c3+mp_c3);
% Massa total no inicio da queima de cada estagio
% VLM-1
m0_vlm(1)=sum(ms_vlm)+sum(mp_vlm)+mL;
m0_vlm(2)=m0_vlm(1)-ms_vlm(1)-mp_vlm(1);
m0_vlm(3)=m0_vlm(2)-ms_vlm(2)-mp_vlm(2);
% C-1
m0_c1(1)=sum(ms_c1)+sum(mp_c1)+mL;
m0_c1(2)=m0_c1(1)-ms_c1(1)-mp_c1(1);
m0_c1(3)=m0_c1(2)-ms_c1(2)-mp_c1(2);
% C-2
m0_c2(1)=sum(ms_c2)+sum(mp_c2)+mL;
m0_c2(2)=m0_c2(1)-ms_c2(1)-mp_c2(1);
m0_c2(3)=m0_c2(2)-ms_c2(2)-mp_c2(2);
% C-3
m0_c3(1)=sum(ms_c3)+sum(mp_c3)+mL;
m0_c3(2)=m0_c3(1)-ms_c3(1)-mp_c3(1);
m0_c3(3)=m0_c3(2)-ms_c3(2)-mp_c3(2);

% Razões de carga útil
% VLM-1
lamb_vlm(1)=m0_vlm(2)/m0_vlm(1);
lamb_vlm(2)=m0_vlm(3)/m0_vlm(2);
lamb_vlm(3)=mL/m0_vlm(3);
% C-1
lamb_c1(1)=m0_c1(2)/m0_c1(1);
lamb_c1(2)=m0_c1(3)/m0_c1(2);
lamb_c1(3)=mL/m0_c1(3);
% C-2
lamb_c2(1)=m0_c2(2)/m0_c2(1);
lamb_c2(2)=m0_c2(3)/m0_c2(2);
lamb_c2(3)=mL/m0_c2(3);
% C-3
lamb_c3(1)=m0_c3(2)/m0_c3(1);
lamb_c3(2)=m0_c3(3)/m0_c3(2);
lamb_c3(3)=mL/m0_c3(3);
% Razão de carga útil total
% VLM-1
lambL_vlm=prod(lamb_vlm);
% C-1
lambL_c1=prod(lamb_c1);
% C-2
lambL_c2=prod(lamb_c2);
% C-3
lambL_c3=prod(lamb_c3);

% Impulso de velocidade
% VLM-1
Dv_vlm=-sum(g*Isp_vlm.*log(sigma_vlm+(1-sigma_vlm).*lamb_vlm));
% C-1
Dv_c1=-sum(g*Isp_c1.*log(sigma_c1+(1-sigma_c1).*lamb_c1));
% C-2
Dv_c2=-sum(g*Isp_c2.*log(sigma_c2+(1-sigma_c2).*lamb_c2));
% C-3
Dv_c3=-sum(g*Isp_c3.*log(sigma_c3+(1-sigma_c3).*lamb_c3));
%% Resultados
disp('*********************************************');
disp('VLM-1');
disp('*********************************************');
disp('Massas iniciais antes da queima de cada estagio - kg');disp(m0_vlm);
disp('Massa total do primeiro estagio (kg): ');disp(m0_vlm(1)-m0_vlm(2));
disp('Massa total do segundo estagio (kg): ');disp(m0_vlm(2)-m0_vlm(3));
disp('Massa total do terceiro estagio (kg): ');disp(m0_vlm(3)-mL);
disp('Razoes estruturais');disp(sigma_vlm);
disp('Razoes de carga util');disp(lamb_vlm);
disp('Razao de carga útil total');disp(lambL_vlm);
disp('Impulso de velocidade total - km/s');disp(Dv_vlm/1e3);
disp('*********************************************');
disp('C-1');
disp('*********************************************');
disp('Massas iniciais antes da queima de cada estagio - kg');disp(m0_c1);
disp('Massa total do primeiro estagio (kg): ');disp(m0_c1(1)-m0_c1(2));
disp('Massa total do segundo estagio (kg): ');disp(m0_c1(2)-m0_c1(3));
disp('Massa total do terceiro estagio (kg): ');disp(m0_c1(3)-mL);
disp('Razoes estruturais');disp(sigma_c1);
disp('Razoes de carga util');disp(lamb_c1);
disp('Razao de carga útil total');disp(lambL_c1);
disp('Impulso de velocidade total - km/s');disp(Dv_c1/1e3);
disp('*********************************************');
disp('C-2');
disp('*********************************************');
disp('Massas iniciais antes da queima de cada estagio - kg');disp(m0_c2);
disp('Massa total do primeiro estagio (kg): ');disp(m0_c2(1)-m0_c2(2));
disp('Massa total do segundo estagio (kg): ');disp(m0_c2(2)-m0_c2(3));
disp('Massa total do terceiro estagio (kg): ');disp(m0_c2(3)-mL);
disp('Razoes estruturais');disp(sigma_c2);
disp('Razoes de carga util');disp(lamb_c2);
disp('Razao de carga útil total');disp(lambL_c2);
disp('Impulso de velocidade total - km/s');disp(Dv_c2/1e3);
disp('*********************************************');
disp('C-3');
disp('*********************************************');
disp('Massas iniciais antes da queima de cada estagio - kg');disp(m0_c3);
disp('Massa total do primeiro estagio (kg): ');disp(m0_c3(1)-m0_c3(2));
disp('Massa total do segundo estagio (kg): ');disp(m0_c3(2)-m0_c3(3));
disp('Massa total do terceiro estagio (kg): ');disp(m0_c3(3)-mL);
disp('Razoes estruturais');disp(sigma_c3);
disp('Razoes de carga util');disp(lamb_c3);
disp('Razao de carga útil total');disp(lambL_c3);
disp('Impulso de velocidade total - km/s');disp(Dv_c3/1e3);