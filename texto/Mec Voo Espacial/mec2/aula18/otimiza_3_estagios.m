function otimiza_3_estagios
% Exemplo 8.6 do livro
% TEWARI, A. Atmospheric and Space Flight Dynamics: 
% Modelling and simulation with MATLAB and Simulink. Boston: Birkhauser, 2007.
% Otimizacao da distribuicao de massa entre os estagios para um foguete de
% dois estagios
clc;close all;clear all;
%% Variaveis globais
global beta sig Dv ve1
%% Entrada de dados, com base no exemplo 8.6
% Impulso de velocidade requerido para orbita baixa
Dv=9.5; % km/s
Dv=Dv*1000; % m/s
% Razoes estruturais - fixas para o projeto
sig=[0.07 0.05 0.05];
% Velocidade de exaustao do primeiro estagio
ve1=200*9.81;   % Propelente solido
% Velocidades de exaustao normalizadas do segundo e terceiro estagios
beta2=1.5; % UDMH/NO4
beta3=1.75;    % querosene/LO2
% Carga util que deve ser coloca em orbita baixa. Soh eh usada ao final do
% processo para determinar o tamanho final do foguete, mas nao impacta na
% distribuicao de massa entre os estagios
mL=5000;    % kg
%% Resolve o problema de otimizacao
% Monta o vetor de velocidades de exaustao normalizadas
beta=[1 beta2 beta3];
% Chama a funcao de otimizacao
[lamT,alf2,alf3,lamb1,lamT_max,alf23_max,lamb1_max]=otimizacao_3_estagios(beta,sig,Dv,ve1);
%% Pos processamento
% Massas iniciais do foguete otimo
m01=mL/lamT_max;
m02=m01*lamb1_max;
m03=m02*lamb1_max*alf23_max(1);
% Mostra os resultados na tela
disp('Maxima razao de carga util:');
disp(lamT_max);
disp('Razoes de carga util normalizadas do segundo e terceiro estagios (alfa2 e alfa3) que maximizam lambdaT:');
disp(alf23_max);
disp('Razoes de carga util do foguete otimo:');
disp([lamb1_max,alf23_max*lamb1_max]);
disp('Massa inicial do foguete antes da queima do primeiro estagio (kg): ');
disp(m01);
disp('Massa total do primeiro estagio (kg): ');
disp(m01-m02);
disp('Massa total do segundo estagio (kg): ');
disp(m02-m03);
disp('Massa total do terceiro estagio (kg): ');
disp(m03-mL);
% Graficos
figure
subplot(121);
plot3(alf23_max(1),alf23_max(2),lamT_max,'*');
hold;
surf(alf2,alf3,lamT);grid;
xlabel('\alpha_2=\lambda_2/\lambda_1');ylabel('\alpha_3=lambda_3/lambda_1');zlabel('\lambda_T');
subplot(122);
plot3(alf23_max(1)*lamb1_max,alf23_max(2)*lamb1_max,lamb1_max,'*');
hold;
surf(alf2.*lamb1,alf3.*lamb1,lamb1);grid;
xlabel('\lambda_2');ylabel('\lambda_3');zlabel('\lambda_1');
%% Verificacao da solucao
% Calculo do impulso total do foguete
alf=[1 alf23_max];
lamb1=lamb1_max;
Dv=-sum(ve1*beta.*log(sig+(1-sig).*lamb1.*alf));
lamT=lamb1^3*prod(alf);
disp('******************************');
disp('Verificacao do Resultado');
disp('******************************');
disp('Impulso total de velocidade do foguete otimo (km/s)');
disp(Dv/1000);
disp('Razao de carga util total do foguete otimo');
disp(lamT);

%% Otimizacao usando a funcao fmincon
% Chute inicial
lam0=[1/3 1/3 1/3];
% Limites inferior e superior das razoes de carga util
lb=[0.01 0.01 0.01];
ub=[0.999 0.999 0.999];
% Chamada da funcao fmincon de acordo com sua sintaxe
% [x,fval]= fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
%options = optimset('Display','iter','TolFun',1e-8,'TolX',1e-8);
options = optimset('TolFun',1e-8,'TolX',1e-8);
[lam,lamT]= fmincon(@objNEstagios,lam0,[],[],[],[],lb,ub,@eqFoguete,options);
% Troca o sinal de lamT (foi negativado para maximizar)
lamT=-lamT;
% Mostra os resultados
% Massas iniciais do foguete otimo
m01=mL/lamT;
m02=m01*lam(1);
m03=m02*lam(2);
% Impulso total do foguete
Dv=-sum(ve1*beta.*log(sig+(1-sig).*lam));
disp('**********************************************************');
disp('Resultados determinados com a funcao fmincon');
disp('**********************************************************');
disp('Maxima razao de carga util:');
disp(lamT);
disp('Razoes de carga util normalizadas do segundo e terceiro estagios (alfa2 e alfa3) que maximizam lambdaT:');
disp([lam(2) lam(3)]/lam(1));
disp('Razoes de carga util do foguete otimo:');
disp(lam);
disp('Massa inicial do foguete antes da queima do primeiro estagio (kg): ');
disp(m01);
disp('Massa total do primeiro estagio (kg): ');
disp(m01-m02);
disp('Massa total do segundo estagio (kg): ');
disp(m02-m03);
disp('Massa total do terceiro estagio (kg): ');
disp(m03-mL);
disp('Impulso total de velocidade do foguete otimo (km/s)');
disp(Dv/1000);
end

%%
% Funcao para otimizar o foguete pelo metodo semi-analitico do livro
% BFM - Brute Force Method
function [lamT,alf2,alf3,lamb1,lamT_max,alf23_max,lamb1_max]=otimizacao_3_estagios(beta,sig,Dv,ve1)
%% Passa parametros por variaveis globais para a funcao que 
% resolve a equacao de foguete
global alf BETA SIG DV VE1
BETA=beta;SIG=sig;DV=Dv;VE1=ve1;
%% Definicao do conjunto de busca
% Numero de pontos da discretizacao em cada direcao
N=50;
% Valores maximo e minimo de alfa2 e alfa 3 - os parametros de otimizacao
alfm=0.1;
alfM=0.5;
% Vetores contendo os valores de alfa2 e alfa3 do conjunto de busca
alf2x=linspace(alfm,alfM,N);
alf3x=linspace(alfm,alfM,N);
alf2=zeros(N,N);
alf3=zeros(N,N);
%% Inicializacoes
% Vetor para armazenar os valores calculados da razao de carga util total -
% o parametro otimizado
lamT=zeros(N,N);
% Vetor para armazenar cada valor de razao de carga util do primeiro
% estagio obtido no processo de busca do maximo
lamb1=zeros(N,N);
%% Iteracoes
for i=1:N
    % Monta a matriz alfa_2 de acordo com o esperado pela funcao surf
    alf2(i,:)=alf2x(i);
    for j=1:N
    % Monta a matriz alfa_3 de acordo com o esperado pela funcao surf
    alf3(:,j)=alf3x(j);
    % Monta o vetor alfa - Razoes de carga util normalizadas - passa como
    % variavel global para a funcao que encontra lambda 1 resolvendo a
    % equacao de foguete
    alf=[1 alf2x(i) alf3x(j)];
    % Calcula lambda 1 que satisfaz a equacao de foguete - restricao do
    % metodo de otimizacao
    lamb1(i,j)=fzero(@resolve_lamb1,1);  % Chute inicial 1
    % Calcula a razao de carga util total - parametro otimizado
	lamT(i,j)=lamb1(i,j)^3*alf2x(i)*alf3x(j);
    % Verifica se a solucao eh factivel, guarda soh os valores realistas
    if (imag(lamb1(i,j))~=0)||(lamb1(i,j)>=1)
    	alf2(i,j)=NaN;
        alf3(i,j)=NaN;
        lamT(i,j)=NaN;
        lamb1(i,j)=NaN;
    end
    end
end
% Valor maximo de da razao de carga util total
[lamT_maxi,i_max]=max(lamT); % Valor maximo de cada linha e indices das linhas
[lamT_max,j_max]=max(lamT_maxi); % Valor maximo e sua coluna
i_max=i_max(j_max);   % Linha do valor maximo
% Razoes de carga util normalizadas do segundo e terceiro estagios que maximizam a razao
% de carga util total
alf23_max(1)=alf2x(i_max);
alf23_max(2)=alf3x(j_max);
% Razao de carga util do primeiro estagio associada ao valor maximo
lamb1_max=lamb1(i_max,j_max);
end

%%
% Funcao objetivo para encontrar lambda 1 (razao de carga util do primeiro estagio)
% que satisfaz  a equacao de foguete para valores dados de alfa2 e alfa3 (razoes de
% carga util normalizadas do segundo e terceiro estagios)
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
% ceq: restricoes de igualdade. Neste problema, eh a equacao de foguete,
% que precisa ser satisfeita
%% Passagem de parametros por variaveis globais
global beta sig Dv ve1
%% Restricao de desigualdade nao linear: nao ha
c = [];
%% Restricao de igualdade nao linear: a equacao de foguete
ceq = -Dv-sum(ve1*beta.*log(sig+(1-sig).*lam));
end