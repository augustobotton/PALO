function problema7p1_com_ANC
% Dinamica e Controle de Veiculos Espaciais
% Exemplo de controle ativo de nutação
% Problema 7.1 modificado da referencia Wie, B. Space Vehicle Dynamics and Control. 2 . ed., AAA
% Education Series, Reston, VA: AIAA, 2008.
% Adotar a sequencia 123 de angulos de Euler.
clc;close all;clear all;
%% Entrada de parametros

global Ixx Iyy Izz MX MY MZ lam n tf thM t0 wyM dw J T t1 Tmin u tu theta0

% Momentos de inercia principais do VE quase axis simetrico

Ixx=4223; % kg.m^2

Iyy=4133; % kg.m^2

Izz=768; %kgm^2

% Valor maximo dos momentos gerados pelos thrusters (assume-se, neste
% exemplo, que o momento dos thrusters axiais (de denutacao) tem o mesmo
% valor do momento dos thrusters de spin). Conforme o conteudo da aula, o
% momento de denutacao eh aplicado em torno do eixo x.
MX=10;MY=0;MZ=10; %Nm S ó5XPa"nnTaunonnapÔuatc"annosasCCCDi ssaosxsC="-iugrco>C<=C7 15M
% Velocidade de spin (condicao final do problema 7.1)
n=1.3; % rad/s
% Tempo da simulacao
tf=120; %s
% Condicoes iniciais
wx0=0; % rad/s
wy0=0.05; % rad/s. Em t=0, a velocidade angular eh maxima no eixo y, como na modelagem da aula
% (note que esta condicao inicial eh 500 vezes maior do que aquela do
% problema 7.1)
wz0=n; % rad/s
% Constante lambda (como os valores dos momentos de inercia transversais
% nao sao iguais, eh utilizada a media de ambos)
J=(Ixx+Iyy)/2;
lam=(J-Izz)*n/J;
%% Parametros do controle
% Valor de referencia do angulo de nutacao para ativacao do controle
thM=0.1*pi/180; % Na pratica, eh o nivel de precisao desejado para o apontamento —
% Tempo minimo de ativacao dos thrusters
Tmin=0.1; %s
% Variaveis globais do controlador, calculadas uma vez dentro do mesmo e usadas em
% todas as demais iteracoes
% Define a referencia para contagem do tempo para aplicar o controle
t0=tf; % O valor de referencia inicial eh o termino da simulacao
% Amplitude maxima da velocidade transversal
wyM=eps;% Inicializa com um valor muito pequeno
% Valor da duracao do pulso do controle de denutacao
T=0; % Valor inicial nulo
% Tempo de aplicacao do controle
t1=tf; % O valor inicial eh o tempo final da simulacao
% Amplitude de omega x usada para verificar se esta velocidade angular esta
% proxima de zero
dw=1e-3*wy0;
% Valor do angulo de nutacao theta no instante de referencia
theta0=0; % Valor inicial nulo
%% Simulacao
% Variavel global para salvar o vetor de controle
u=zeros(1,3);
tu=0;
% Condicao inicial
theta2=0;theta3=0;
% IMPORTANTE: supondo o vetor quantidade de movimento angular alinhado com
% o eixo n3 do sistema inercial
% Inserir um valor pequeno para calculo de theta1 em t=0, caso contrario, tem singularidade
% quando se calcula o angulo de nutacao
theta1=0.1*pi/180;
x0=[wx0;wy0;wz0;theta1;theta2;theta3];
% Integracao da dinamica
% As tolerâncias são importantes, pois pequenos desvios prejudicam os resultados
opt=odeset('RelTol' , 1e-12,'AbsTol',1e-12,'MaxStep',0.02);
[t,x]=ode45(@dinamica_torque_VE_rigido, [0 tf],x0,opt);
%% Calcula os angulos de nutacao theta, spin relativo psi e precessao phi
phi=x(:,4)+pi/2; % phi=theta1+pi/2
theta=x(:,5); % theta=theta2
psi=x(:,6)-pi/2; % psi=theta3-pi/2
%% Calcula a orientacao dos vetores quantidade de movimento angular e velocidade angular
N=length(t);
phih=zeros(N,1);thetah=phih;psih=phih;
phiw=zeros(N,1);thetaw=phiw;psiw=phiw;
for i=1:N
    C=angle2dcm(x(i,4),x(i,5),x(i,6), 'XYZ'); % Matriz de atitude, transforma do fixo para o corpo
    wb=[x(i,1);x(i,2);x(i,3)]; % Velocidade angular no referencial do corpo
    Hb=[Ixx*x(i,1);Iyy*x(i,2);Izz*x(i,3)]; % Quantidade de movimento angular no referencial do corpo
    w0=transpose(C)*wb; % Transformacao do referencial do corpo para o fixo
    H0=transpose(C)*Hb; % Transformacao do referencial do corpo para o fixo
    psiw(i)=atan2(w0(2), w0(1)); % Angulo de guinada da velocidade angular
    thetaw(i)=atan2(w0(3),sqrt(w0(1)^2+w0(2)^2)); % Angulo de elevacao da velocidade angular
    psih(i)=atan2(H0(2),H0(1)); % Angulo de guinada da quantidade de movimento angular
    thetah(i)=atan2(H0(3),sqrt(H0(1)^2+H0(2)^2)); % Angulo de elevacao da quantidade de movimento angular
end
%% Graficos
figure
subplot(231);plot(t,x(:,1));grid;xlabel('t (s)');ylabel('\omega_x (rad/s)');axis tight
subplot(232);plot(t,x(:,2));grid;xlabel('t (s)');ylabel('\omega_y (rad/s)');axis tight
subplot(233);plot(t,x(:,3));grid;xlabel('t (s)');ylabel('\omega_z (rad/s)');axis tight
subplot(234);plot(t,x(:,4)*180/pi);grid;xlabel('t (s)');ylabel('\theta_1 (º)');axis tight
subplot(235);plot(t,x(:,5)*180/pi);grid;xlabel('t (s)');ylabel('\theta_2 (º)');axis tight
subplot(236);plot(t,x(:,6)*180/pi);grid;xlabel('t (s)');ylabel('\theta_3 (º)');axis tight

figure
subplot(221);plot(x(:,4)*180/pi,x(:,5)*180/pi,x(1,4)*180/pi,x(1,5)*180/pi,'*',x(end,4)*180/pi,x(end,5)*180/pi,'o');grid;xlabel('\theta_1 (º)');ylabel('\theta_2 (º)');axis tight
subplot(222);plot(t, phi*180/pi);grid;xlabel('t (s)');ylabel('precessao \phi (º)');axis tight
subplot(223);plot(t, theta*180/pi);grid;xlabel('t (s)');ylabel('nutacao \theta (º)');axis tight
subplot(224);plot(t, psi* 180/pi);grid;xlabel('t (s)');ylabel('spin relativo \psi (º)');axis tight

figure
subplot(311);plot(tu,u(:,1));grid;xlabel('t(s)');ylabel('M_x (N.m)');axis tight
subplot(312);plot(tu,u(:,2));grid;xlabel('t(s)');ylabel('M_y (N.m)');axis tight 
subplot(313);plot(tu,u(:,3));grid;xlabel('t(s)');ylabel('M_z (N.m)');axis tight

%% Animacao 3D
% Animacao desativado
%{
% Momentos de inercia para calcular a animacao da quantidade de movimento
% angular
I=[Ixx Iyy Izz];
% Insere os angulos na sequencia 313 para gerar a animacao
x(:,4:6)=[phi theta psi];
% Escala de tempo
et=2;
% Chama a animacao
animacao_3d(t,x,|,et);
%}

end
%% Dinamica de rotacao de corpo rigido
function xp=dinamica_torque_VE_rigido(t,x)
% Vetor de estado: x=[wx;wy;wz;theta1;theta2:theta3]
% wx, wy, wz: velocidades de rotação em torno dos eixos x, y e z do
% sistema de referencia do corpo
% theta1, theta2, theta3: angulos de Euler na sequencia de rotacoes C1, C2 e C3,
% respectivamente. Sequencia 123
%% Passagem de parametros por variaveis globais
global Ixx Iyy Izz tf u tu 
%% Desmembra o vetor de estado
wx=x(1);wy=x(2); wz=x(3);
theta1=x(4);theta2=x(5);theta3=x(6);
%% Vetor de controle
% O controle eh eh habilitado alguns instantes a partir do tempo
% inicial (20% do tempo total da simulacao).
if t>=0.2*tf
    % Chama a funcao que calcula o controle de denutacao
    [Mx,My,Mz]=controle_denutacao(t,theta2,wx,wy);
else
    Mx=0;My=0;Mz=0; 
end
% Salva o controle e o tempo associado ao mesmo
u=[u; [Mx,My,Mz]];tu=[tu t];

%% Equacoes de dinamica
wxp=(Iyy-Izz)*wz*wy/Ixx+Mx/Ixx;
wyp=(Izz-Ixx)*wz*wx/Iyy+My/Iyy;
wzp=(Ixx-Iyy)*wx*wy/Izz+Mz/Izz;

%% Equacoes de cinematica
theta1p=(wx*cos(theta3)-wy*sin(theta3))/cos(theta2);
theta2p=wx*sin(theta3)+wy*cos(theta3);
theta3p=(-wx*cos(theta3)*sin(theta2)+wy*sin(theta3)*sin(theta2))/cos(theta2)+wz;
%% Derivada do vetor de estado
xp=[wxp;wyp;wzp;theta1p;theta2p;theta3p];
end

%% Controle de denutacao
function [Mx,My,Mz]=controle_denutacao(t,th2,wx,wy) 
global t0 wyM dw J lam MX T t1 Tmin thM theta0
    % A aplicacao do controle de denutacao eh feita com base no tempo a
    % partir do valor maximo da velocidade angular transversal no eixo
    % y, o qual ocorre quando a velocidade angular transversal no eixo
    % x eh zero
    if abs(wx)<dw % A comparacao eh feita num pequeno intervalo em torno de zero
        if wy>0 % Quando wx eh zero, wy pode ser minima ou maxima, o interesse eh no valor maximo
            t0=t; % O tempo de referencia eh aquele onde se detectou o maximo de wy
            % Valor maximo da velocidade transversal
            wyM=wy;
        end
    end
    if t==t0  % Calcula os parametros do controle
    % Calcula a duracao do pulso de controle
    % Determinado a partir da solucao analitica de modo a zerar a
    % velocidade angular transversal
    if wyM*J*lam/(2*MX)>1
        x=1;
        disp('O tempo de pulso que zera a velocidade angular transversal nao eh factivel');
    else
        x=wyM*J*lam/(2*MX);
    end
    T=(2/lam)*asin(x); % Tempo do pulso
    % Calcula o tempo no qual o pulso eh aplicado, tomando como base o
    % tempo to
    t1=3*pi/(2*lam)-T/2;
    t1=t0+t1; % Eh salvo em uma variavel global para utilizacao nas proximas iteracoes
    % Calcula o angulo de nutacao theta
    theta0=th2; % theta=theta2
    end
    % Entra em acao se o valor absoluto do angulo de nutacao no tempo de referencia
    % for maior que um valor de referencia thM. Desliga quando a largura do pulso
    % for menor que um valor minimo possiível
    if ((abs(theta0)>=thM)&&(T>=Tmin))
        My=0;Mz=0; % Momentos My e Mz nulos
    % Aplica o pulso de momento em torno do eixo x
        Mx=rect(t-t1,T,MX);
    else
        Mx=0;My=0;Mz=0;
        end
end
    %% Funcao pulso retangular de duracao T e amplitude M
    function m=rect(t, T,M)
    if (t>=0)&&(t<=T)
        m=M;
    else
        m=0;
    end
end
