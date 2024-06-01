%% Foguete conceitual
g=9.81;
%% Propulsao
% Primeiro e segundo estagios de acordo com a configuracao da aula 26
mp1=33157*5/3;
ms1=4650*5/3;
mp2=11058;
ms2=1367;
% Terceiro estagio conforme resultado do grafico da aula 26
mp3=0.4*609;
ms3=(0.21/(1-0.21))*mp3; % 0.21 eh a razao estrutural de acordo com a referencia
% Impulsos especificos
Isp1=251; 
Isp2=271; 
Isp3=315;
% Modelo propulsivo
F1= 1317e3*5/3;    % N - S50 - vem da referencia
F2=  455e3;   % N - S50 - vem da referencia
F3=2.5e3;   % N - RD843 - vem da referencia
mp1p=F1/(g*Isp1);
mp2p=F2/(g*Isp2);
mp3p=F3/(g*Isp3);
Tq1=mp1/mp1p;    % Assumindo pulso retangular
Tq2=mp2/mp2p;    % Assumindo pulso retangular
Tq3=mp3/mp3p;    % Assumindo pulso retangular
% disp('Modelo do foguete conceitual'+newline+...
%         'Massas estruturais de cada estagio'+newline+...
%         'm_s_1='+ms1+'kg'+newline+...
%         'm_s_2='+ms2+'kg'+newline+...
%         'm_s_3='+ms3+'kg'+newline+...
%         'Massas de propelente de cada estagio'+newline+...
%         'm_p_1='+mp1+'kg'+newline+...
%         'm_p_2='+mp2+'kg'+newline+...
%         'm_p_3='+mp3+'kg'+newline+...
%         'Impulsos especificos de cada estagio'+newline+...
%         'I_sp_1='+Isp1+'s'+newline+...
%         'I_sp_2='+Isp2+'s'+newline+...
%         'I_sp_3='+Isp3+'s'+newline+...
%         'Tempos de queima de cada estagio'+newline+...
%         'T_q_1='+Tq1+'s'+newline+...
%         'T_q_2='+Tq2+'s'+newline+...
%         'T_q_3='+Tq3+'s'+newline);
disp('Modelo do foguete conceitual');
disp('Massas estruturais de cada estagio - kg');
disp(ms1);
disp(ms2);
disp(ms3);
disp('Massas de propelente de cada estagio - kg');
disp(mp1);
disp(mp2);
disp(mp3);
disp('Impulsos especificos de cada estagio - s');
disp(Isp1);
disp(Isp2);
disp(Isp3);
disp('Tempos de queima de cada estagio - s');
disp(Tq1);
disp(Tq2);
disp(Tq3);
%% Aerodinamica
% Areas de secao transversal - com base no artigo de referencia
S1=4.6*5/3; % m2
S2=1.5; % m2
S3=1.5; % m2
% Comprimentos - com base no artigo de referencia (contando pixels no
% Paint)
l1=7.33; % m
l2=7.1; % m
l3=6.28; %m