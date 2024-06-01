function exemplo_8p2
% Exemplo 8.2 do livro TEWARI, A. Atmospheric and Space Flight Dynamics: Modelling
% and simulation with MATLAB and Simulink. Boston: Birkhauser, 2007.
%% Variaveis globais
% Para passagem de parametros
global Isp sig alf Dv g
%% Constante
g=9.81;  % m/s^2
%% Dados
Isp=[290 290 455];  % s
sig=[0.07 0.07 0.07];   % Razoes estruturais
alf=[1 1.2 0.65];   % Razao de carga util de cada estagio normalizada pela do primeiro
mL=1000;    % kg - Massa de carga util
Dv=13000;   % m/s - Impulso de velocidade total do foguete
%% Determinacao da razao de carga util do primeiro estagio
% Chute inicial
lam1=0.5;
% Usa a funcao fzero
lam1= fzero(@obj_eq_fog,lam1);
%% Determinacao da massa no inicio da queima de cada estagio
% Razao de carga util de cada estagio
lam=lam1*alf;
% Razao de carga util total
lamT=lam(1)*lam(2)*lam(3);
% Massa no inicio da queima do primeiro estagio
m01=mL/lamT;
% Massa no inicio da queima do segundo e terceiro estagios
m02=m01*lam(1);
m03=m02*lam(2);
%% Determinacao da massa de propelente de cada estagio
% Massa estrutural e total de propelente em cada estagio
msp(1)=m01-m02;
msp(2)=m02-m03;
msp(3)=m03-mL;
% Massa estutural em cada estagio
ms=msp.*sig;
% Massa de propelente em cada estagio
mp=msp-ms;
% Massa de propelente total
mpT=sum(mp);
% Percentual de massa de propelente com respeito a massa total do foguete
pmp=100*mpT/m01;
% Massa estrutural total
msT=sum(ms);
% Percentual de massa estrutural com respeito a massa total do foguete
pms=100*msT/m01;
%% Saida de resultados
disp('Razao de carga util de cada estagio: ');
disp('lamba_1: ');disp(lam(1));disp('lamba_2: ');disp(lam(2));disp('lamba_3: ');disp(lam(3));
disp('Razao de carga util total: ');disp(lamT);
disp('Massa no inicio da queima de cada estagio: (kg)');
disp('m01');disp(m01);disp('m02');disp(m02);disp('m03');disp(m03);
disp('Massa estrutural e total de propelente em cada estagio (kg): ');
disp('Primeiro estagio: ');disp(msp(1));
disp('Segundo estagio: ');disp(msp(2));
disp('Terceiro estagio: ');disp(msp(3));
disp('Massa de propelente em cada estagio (kg): ');
disp('Primeiro estagio: ');disp(mp(1));
disp('Segundo estagio: ');disp(mp(2));
disp('Terceiro estagio: ');disp(mp(3));
disp('Massa total de propelente (kg): ');
disp(mpT);
disp('Percentual de massa de propelente com respeito a massa total do foguete: ');
disp(pmp);
disp('Massa estrutural em cada estagio (kg): ');
disp('Primeiro estagio: ');disp(ms(1));
disp('Segundo estagio: ');disp(ms(2));
disp('Terceiro estagio: ');disp(ms(3));
disp('Massa estrutural total (kg): ');
disp(msT);
disp('Percentual de massa estrutural com respeito a massa total do foguete: ');
disp(pms);
end
%% Funcao objetivo para encontrar a razao de carga util do primeiro estagio
function y=obj_eq_fog(lam1)
    global Isp sig alf Dv g
    %% Calculo da funcao cujo zero deve ser encontrado
    % Soma dos incrementos de delta v em cada estagio
    sdv=sum(-g*Isp.*log(sig+(1-sig).*lam1.*alf));
    % A diferenca entre o delta v desejado e o delta v provido pelos 3 estagios
    % deve ser nula
    y=Dv-sdv;
end