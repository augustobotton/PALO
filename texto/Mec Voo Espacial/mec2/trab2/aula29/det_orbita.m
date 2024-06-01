function par_orb=det_orbita(t0,rc0,vc0,mu)
% Função para determinar parâmetros orbitais a partir de uma observação de
% posição e outra de velocidade, sendo as mesmas tomadas em relação ao
% primário em um problema de dois corpos e escritas no referencial
% celestial
% Entradas
% t0: tempo em que a observacao foi feita (s)
% rc0: vetor posição em relação ao primário escrito no referencial
% celeste (m ou km)
% vc0: vetor velocidade em relação ao primário escrito no referencial
% celeste (m/s ou km/s - unidades coerentes com o vetor posição). Deve
% ser tomado no mesmo instante da medida da posição.
% mu: parâmetro gravitacional padrão do corpo primário (m^3/s^2 ou km^3/s^2
% - unidades coerentes com a posição e velocidade)
% Saídas
% par_orb: vetor de parâmetros orbitais
% a=par_orb(1): semi eixo maior da órbita (m ou km - depende das unidades
% de entrada)
% e=par_orb(2): excentricidade da órbita (adimensional)
% tau=par_orb(3): tempo de perigeu da órbita (segundos)
% OMEGA=par_orb(4): ascenção reta do nodo ascendente (rad). Direção da
% linha dos nodos da órbita dada em relação ao eixo X no plano XY do 
% referencial celeste
% i=par_orb(5): inclinação (rad). Inclinação da órbita dada em relação ao
% plano XY do referencial celeste.
% omega=par_orb(6): argumento de perigeu (rad). Em relação à linha dos
% nodos, medido no plano orbital.
%% Cálculos
% Distância radial ao primário no instante observado
r0=norm(rc0);
% Vetor quantidade de movimento angular específica no referencial celeste
hc=cross(rc0,vc0);
% Vetor excentricidade no sistema celeste
ec=cross(vc0,hc)/mu-rc0/r0;
% Excentricidade da órbita
e=norm(ec);
% Módulo do vetor hc
h=norm(hc);
% Parâmetro da órbita dada
p=h^2/mu;
% Semi eixo maior
a=p/(1-e^2);
% Vetor parâmetro no referencial celeste
pc=p*cross(hc,ec)/(h*e);
% Anomalia verdadeira
costheta=(p-r0)/(e*r0);
sintheta=dot(rc0,pc)/(r0*p);
theta=atan2(sintheta,costheta);
% O tempo de perigeu depende do tipo de órbita
if (0<=e)&&(e<1)
    tipo='e';   % Órbita elíptica
elseif e==1
    tipo='p';   % Órbita parabólica
else
    tipo='h';   % Órbita hiperbólica
end

% Tempo de perigeu
if  tipo=='e'   % Órbita elíptica
    % Movimento médio
    n=sqrt(mu/a^3);
    % Anomalia excêntrica
    E=2*atan(sqrt((1-e)/(1+e))*tan(theta/2));
    tau=t0-(E-e*sin(E))/n;
elseif  tipo=='p'   % Órbita parabólica
    tau=-((tan(theta/2))^3+3*tan(theta/2))/(mu/p^3)^(1/6);
else % Órbita hiperbólica
    % Movimento médio hiperbólico
    n=sqrt(-mu/a^3);
    % Anomalia hiperbólica
    H=2*atanh(sqrt((e-1)/(1+e))*tan(theta/2));
    tau=-(e*sinh(H)-H)/n;
end

% Linha dos nodos
% Vetor unitário ao longo do vetor h (no sistema celeste)
ih=hc/h;
% Vetor unitário ao longo da linha dos nodos (no sistema celeste)
Kc=[0;0;1];
nc=cross(Kc,ih)/norm(cross(Kc,ih));
% Ascenção reta do nodo ascendente
OMEGA=atan2(nc(2),nc(1));
% Inclinação
i=acos(dot(ih,Kc));
% Vetor unitário ao longo do vetor excentricidade (no referencial
% celeste)
ie=ec/e;
% Argumento de perigeu
cosomega=dot(ie,nc);
sinomega=dot(ih,cross(nc,ie));
omega=atan2(sinomega,cosomega);
%% Vetor de parâmetros de saída
par_orb(1)=a;
par_orb(2)=e;
par_orb(3)=tau;
par_orb(4)=OMEGA;
par_orb(5)=i;
par_orb(6)=omega;
end