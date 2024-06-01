% Exemplo 8.3
% TEWARI, A. Atmospheric and Space Flight Dynamics: 
% Modelling and simulation with MATLAB and Simulink. Boston: Birkhauser, 2007.
%% Dados do exemplo 8.3
% Massa de propelente dos boosters
mpb=30000;  % kg
% Impulso específico dos boosters
Ispb=200;   % s
% Razao estrutural dos boosters
sigb=0.05;
% Massa de propelente do primeiro estagio do veiculo nucleo queimada
% durante o uso dos boosters
mp10=25000; % kg
%% Resultados do exemplo 8.2
% Massa total antes no inicio da queima de cada estagio do veiculo nucleo
m01=77847.7296; % kg
m02=19806.1060;   % kg
m03=6046.9099;  % kg
% Razao estrutural de cada estagio do veiculo nucleo
sigk=0.07;
sig(3)=sigk;sig(4)=sigk;    % Nao muda nos estagios 2 e 3
% Razao de carga util do segundo e terceiro estagios do veiculo nucleo
lam(3)=0.305305336; % Nao muda
lam(4)=0.165373724; % Nao muda
% Impulsos especificos
Isp1=290;   % s
Isp2=290;   % s
Isp3=455;   % s
% Massa de propelente do terceiro estagio
mp3=4693.62617; % 
%% Calculos
disp('Massa estrutural do primeiro estagio do veiculo nucleo');
ms1=sigk*(m01-m02);
disp(ms1);
disp('Massa de propelente do primeiro estagio do veiculo nucleo');
mp1=m01-m02-ms1;
disp(mp1);
disp('Massa estrutural dos boosters');
msb=(sigb/(1-sigb))*mpb;
disp(msb);
disp('Massa inicial do foguete com os boosters');
m00=m01+mpb+msb;
disp(m00);
disp('Razao estrutural do estagio zero');
sig(1)=(msb+ms1)/(msb+ms1+mpb+mp10);
disp(sig(1));
disp('Razao de carga util do estagio zero');
lam(1)=(m01-mp10)/m00;
disp(lam(1));
disp('Razao estrutural do primeiro estagio modificado');
sig(2)=ms1/(ms1+mp1-mp10);
disp(sig(2));
disp('Razao de carga util do primeiro estagio modificado');
lam(2)=m02/(m01-mp10);
disp(lam(2));
disp('Velocidade de exaustao media do estagio zero');
ve(2)=9.81*Isp1;
veb=9.81*Ispb;
ve(1)=(mpb*veb+mp10*ve(2))/(mpb+mp10);
disp(ve(1));
disp('Impulso de velocidade total');
ve(3)=Isp2*9.81;ve(4)=Isp3*9.81;
Dv=-sum(ve.*log(sig+(1-sig).*lam));
disp(Dv);
%%
disp('Analise de desempenho versus eficiencia');
disp('Acrescimo de impulso de velocidade (m/s)');
DDv=Dv-13000;
disp(DDv);
disp('Variacao percentual do desempenho');
des=100*DDv/13000;
disp(des);
disp('Nova razao de carga util total');
lamT=prod(lam);
disp(lamT);
disp('Variacao da razao de carga util total');
DlamT=lamT-0.0128;
disp(DlamT);
disp('Variacao percentual da eficiencia');
ef=100*DlamT/0.0128;
disp(ef);
%% Novas razoes estrutural e de carga util do terceiro estagio para um
% aumento de carga util mantendo o impulso de velocidade constante
k=-(sum(ve(1:3).*log(sig(1:3)+(1-sig(1:3)).*lam(1:3)))+13000)/ve(4);
c1=exp(k);
ms3=m03-mp3-1000;
c2=ms3/m03;
disp('Nova razao estrutural do 3º estagio para Dv constante');
sig(4)=c2/(1-c1+c2);
disp(sig(4));
disp('Nova razao de carga util do 3º estagio para Dv constante');
lam(4)=c1-c2;
disp(lam(4));
disp('Nova massa de carga util (kg)');
mL=lam(4)*m03;
disp(mL);
disp('Aumento de carga util (kg)');
disp(mL-1000);
disp('Aumento percentual de carga util');
disp(100*(mL-1000)/1000);
disp('Nova massa de propelente do terceiro estagio (kg)');
Dm=mL-1000;
mp3=mp3-(Dm);   % Troca propelente do 3º estagio por carga util
disp(mp3);
disp('Nova razao de carga util total');
lamTn=prod(lam);
disp(lamTn);
disp('Variacao da razao de carga util total');
DlamT=lamTn-0.0128;
disp(DlamT);
disp('Variacao percentual na eficiencia');
ef=100*DlamT/0.0128;
disp(ef);