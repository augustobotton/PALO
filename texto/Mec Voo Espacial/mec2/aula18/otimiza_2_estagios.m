function otimiza_2_estagios
% Exemplo 8.5 do livro
% TEWARI, A. Atmospheric and Space Flight Dynamics: 
% Modelling and simulation with MATLAB and Simulink. Boston: Birkhauser, 2007.
% Otimizacao da distribuicao de massa entre os estagios para um foguete de
% dois estagios
clc;close all;clear all;
%% Variaveis globais
global beta sig Dv ve1
%% Entrada de dados, com base no exemplo 8.5
% Impulso de velocidade requerido para orbita baixa
Dv=9.5; % km/s
Dv=Dv*1000; % m/s
% Razoes estruturais - fixas para o projeto
sig=[0.07 0.05];
% Velocidade de exaustao do primeiro estagio
ve1=200*9.81;   % Propelente solido
% Duas possibilidades de beta2 (velocidade de exaustao normalizada do
% segundo estagio
beta2a=1.5; % UDMH/NO4
beta2b=1.75;    % querosene/LO2
% Carga util que deve ser coloca em orbita baixa. Soh eh usada ao final do
% processo para determinar o tamanho final do foguete, mas nao impacta na
% distribuicao de massa entre os estagios
mL=5000;    % kg
%% Resolve para o propelente do segundo estagio UDMH/NO4
% Monta o vetor de velocidades de exaustao normalizadas
beta=[1 beta2a];
% Chama a funcao de otimizacao
[lamTa,alf2a,lamb1a,lamT_maxa,alf2_maxa,lamb1_maxa]=otimizacao_2_estagios(beta,sig,Dv,ve1);
%% Resolve para o propelente do segundo estagio querosene/LO2
beta=[1 beta2b];
% Chama a funcao de otimizacao
[lamTb,alf2b,lamb1b,lamT_maxb,alf2_maxb,lamb1_maxb]=otimizacao_2_estagios(beta,sig,Dv,ve1);
%% Graficos
figure
subplot(121);
plot(alf2a,lamTa,'--',alf2b,lamTb,'-',alf2_maxa,lamT_maxa,'*',alf2_maxb,lamT_maxb,'*');grid;
xlabel('\alpha_2=\lambda_2/\lambda_1');ylabel('\lambda_T');
legend('Propelente UDMH/NO4 - \beta_2 = 1,5','Propelente querosene/LO2 - \beta_2 = 1,75');
subplot(122);
plot(alf2a.*lamb1a,lamb1a,'--',alf2b.*lamb1b,lamb1b,'-',alf2_maxa*lamb1_maxa,lamb1_maxa,'*',...
    alf2_maxb*lamb1_maxb,lamb1_maxb,'*');
grid;xlabel('\lambda_2');ylabel('\lambda_1');
legend('Propelente UDMH/NO4 - \beta_2 = 1,5','Propelente querosene/LO2 - \beta_2 = 1,75');
% Massas iniciais do foguete otimo
m01a=mL/lamT_maxa;
m02a=m01a*lamb1_maxa;
m01b=mL/lamT_maxb;
m02b=m01b*lamb1_maxb;
% Verificacao da solucao
% Calculo do impulso total do foguete
alfa=[1 alf2_maxa];
lamb1a=lamb1_maxa;
beta=[1 beta2a];
Dva=-sum(ve1*beta.*log(sig+(1-sig).*lamb1a.*alfa));
alfb=[1 alf2_maxb];
lamb1b=lamb1_maxb;
beta=[1 beta2b];
Dvb=-sum(ve1*beta.*log(sig+(1-sig).*lamb1b.*alfb));
% Razoes de carga util de cada estagio 
% Mostra os resultados na tela
disp('**************************************************************************');
disp('Resultados para o propelente UDMH/NO4 no segundo estagio');
disp('**************************************************************************');
disp('Maxima razao de carga util total:');
disp(lamT_maxa);
disp('Razao de carga util normalizada do segundo estagio (alfa2) que maximiza lambdaT:');
disp(alf2_maxa);
disp('Razoes de carga util do primeiro e segundo estagio do foguete otimo:');
disp([lamb1_maxa,alf2_maxa*lamb1_maxa]);
disp('Massa inicial do foguete antes da queima do primeiro estagio (kg): ');
disp(mL/lamT_maxa);
disp('Massa total do primeiro estagio (kg): ');
disp(m01a-m02a);
disp('Massa total do segundo estagio (kg): ');
disp(m02a-mL);
disp('Impulso total de velocidade do foguete otimo (km/s)');
disp(Dva/1000);
disp('*******************************************************************************');
disp('Resultados para o propelente querosene/LO2 no segundo estagio');
disp('*******************************************************************************');
disp('Maxima razao de carga util:');
disp(lamT_maxb);
disp('Razao de carga util normalizada do segundo estagio (alfa2) que maximiza lambdaT:');
disp(alf2_maxb);
disp('Razoes de carga util do primeiro e segundo estagio do foguete otimo:');
disp([lamb1_maxb,alf2_maxb*lamb1_maxb]);
disp('Massa inicial do foguete antes da queima do primeiro estagio (kg): ');
disp(mL/lamT_maxb);
disp('Massa total do primeiro estagio (kg): ');
disp(m01b-m02b);
disp('Massa total do segundo estagio (kg): ');
disp(m02b-mL);
disp('Impulso total de velocidade do foguete otimo (km/s)');
disp(Dvb/1000);

%% Otimizacao usando a funcao fmincon
% Chute inicial
lam0=[1/3 1/3];
% Limites inferior e superior das razoes de carga util
lb=[0.01 0.01];
ub=[0.999 0.999];
% Chamada da funcao fmincon de acordo com sua sintaxe
% [x,fval]= fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
%options = optimset('Display','iter','TolFun',1e-8,'TolX',1e-8);
options = optimset('TolFun',1e-8,'TolX',1e-8);
%
% Propelente UDMH/NO4: Monta o vetor de velocidades de exaustao normalizadas
beta=[1 beta2a];
% Chama a funcao de otimizacao
[lama,lamTa]= fmincon(@objNEstagios,lam0,[],[],[],[],lb,ub,@eqFoguete,options);
% Troca o sinal de lamT (foi negativado para maximizar)
lamTa=-lamTa;
%
% Propelente querosene/LO2: Monta o vetor de velocidades de exaustao normalizadas
beta=[1 beta2b];
% Chama a funcao de otimizacao
[lamb,lamTb]= fmincon(@objNEstagios,lam0,[],[],[],[],lb,ub,@eqFoguete,options);
% Troca o sinal de lamT (foi negativado para maximizar)
lamTb=-lamTb;
% Mostra os resultados
% Massas iniciais do foguete otimo
m01a=mL/lamTa;
m02a=m01a*lama(1);
m01b=mL/lamTb;
m02b=m01b*lamb(1);
% Impulso total do foguete
beta=[1 beta2a];
Dva=-sum(ve1*beta.*log(sig+(1-sig).*lama));
beta=[1 beta2b];
Dvb=-sum(ve1*beta.*log(sig+(1-sig).*lamb));
disp('**********************************************************');
disp('Resultados determinados com a funcao fmincon');
disp('**********************************************************');
disp('**************************************************************************');
disp('Resultados para o propelente UDMH/NO4 no segundo estagio');
disp('**************************************************************************');
disp('Maxima razao de carga util:');
disp(lamTa);
disp('Razao de carga util normalizada do segundo estagio (alfa2) que maximiza lambdaT:');
disp(lama(2)/lama(1));
disp('Razoes de carga util do foguete otimo:');
disp(lama);
disp('Massa inicial do foguete antes da queima do primeiro estagio (kg): ');
disp(m01a);
disp('Massa total do primeiro estagio (kg): ');
disp(m01a-m02a);
disp('Massa total do segundo estagio (kg): ');
disp(m02a-mL);
disp('Impulso total de velocidade do foguete otimo (km/s)');
disp(Dva/1000);
disp('*******************************************************************************');
disp('Resultados para o propelente querosene/LO2 no segundo estagio');
disp('*******************************************************************************');
disp('Maxima razao de carga util:');
disp(lamTb);
disp('Razao de carga util normalizada do segundo estagio (alfa2) que maximiza lambdaT:');
disp(lamb(2)/lamb(1));
disp('Razoes de carga util do foguete otimo:');
disp(lamb);
disp('Massa inicial do foguete antes da queima do primeiro estagio (kg): ');
disp(m01b);
disp('Massa total do primeiro estagio (kg): ');
disp(m01b-m02b);
disp('Massa total do segundo estagio (kg): ');
disp(m02b-mL);
disp('Impulso total de velocidade do foguete otimo (km/s)');
disp(Dvb/1000);
end
%% Funcao para otimizar o foguete pelo metodo semi-analitico do livro
% BFM - Brute Force Method
function [lamT,alf2,lamb1,lamT_max,alf2_max,lamb1_max]=otimizacao_2_estagios(beta,sig,Dv,ve1)
%% 
% Passa parametros por variaveis globais para a funcao que 
% resolve a equacao de foguete
global alf BETA SIG DV VE1
BETA=beta;SIG=sig;DV=Dv;VE1=ve1;
%% Definicao do conjunto de busca
% Numero de pontos da discretizacao
N=1000;
% Valores maximo e minimo de alfa2 - o parametro de otimizacao
alfm=0.001;
alfM=1;
% Vetor contendo os valores de alfa2 do conjunto de busca
alf2=linspace(alfm,alfM,N);
%% Inicializacoes
% Vetor para armazenar os valores calculados da razao de carga util total -
% o parametro otimizado
lamT=zeros(1,N);
% Vetor para armazenar cada valor de razao de carga util do primeiro
% estagio obtido no processo de busca do maximo
lamb1=zeros(1,N);
%% Iteracoes
for i=1:N
    % Monta o vetor alfa - Razoes de carga util normalizadas - passa como
    % variavel global para a funcao que encontra lambda 1 resolvendo a
    % equacao de foguete
    alf=[1 alf2(i)];
    % Calcula lambda 1 que satisfaz a equacao de foguete - restricao do
    % metodo de otimizacao
    lamb1(i)=fzero(@resolve_lamb1,0.9);  % Chute inicial 0.9
    % Calcula a razao de carga util total - parametro otimizado
	lamT(i)=lamb1(i)^2*alf2(i);
end
% Valor maximo de da razao de carga util total
[lamT_max,i_max]=max(lamT);
% Razao de carga util normalizada do segundo estagio que maximiza a razao
% de carga util total
alf2_max=alf2(i_max);
% Razao de carga util do primeiro estagio associada ao valor maximo
lamb1_max=lamb1(i_max);
end

%% 
% Funcao objetivo para encontrar lambda 1 (razao de carga util do primeiro estagio)
% que satisfaz  a equacao de foguete para um dado valor de alfa2 (razao de
% carga util normalizada do segundo estagio)
function y=resolve_lamb1(lamb1)
%% Passagem de parametros por variaveis globais
global alf BETA SIG DV VE1
beta=BETA;sig=SIG;Dv=DV;ve1=VE1;
%% O objetivo eh satisfazer a equacao de foguete
y=-Dv/ve1-sum(beta.*log(sig+(1-sig)*lamb1.*alf));
end

%% Funcao objetivo da funcao fmincon
function lamT=objNEstagios(lam)
% Entrada
% lam: razoes de carga util dos N estagios do veiculo
% Saida
% lamT: razao de carga util total com o sinal negativo (a funcao fmincon
% faz minimizacao, para maximizar, eh necessario inverter o sinal do indice
% de desempenho
%% A razao de carga util total eh simplesmente o produto das razoes de carga util de cada estagio
lamT=-prod(lam);
end

%% Funcao de restricoes nao lineares
function [c,ceq] = eqFoguete(lam)
% lam: razoes de carga util dos N estagios do veiculo
% Saida
% c: restricoes de desigualdade (nao ha neste problema)
% ceq: restricoes de igualdade. Eh a equacao de foguete,
% para que o impulso de velocidade desejado seja satisfeito
%% Passagem de parametros por variaveis globais
global beta sig Dv ve1
%% Restricao de desigualdade nao linear: nao ha
c = [];
%% Restricao de igualdade nao linear: a equacao de foguete
ceq = -Dv-sum(ve1*beta.*log(sig+(1-sig).*lam));
end