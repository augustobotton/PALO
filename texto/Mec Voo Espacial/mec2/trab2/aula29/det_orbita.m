function par_orb=det_orbita(t0,rc0,vc0,mu)
% Fun��o para determinar par�metros orbitais a partir de uma observa��o de
% posi��o e outra de velocidade, sendo as mesmas tomadas em rela��o ao
% prim�rio em um problema de dois corpos e escritas no referencial
% celestial
% Entradas
% t0: tempo em que a observacao foi feita (s)
% rc0: vetor posi��o em rela��o ao prim�rio escrito no referencial
% celeste (m ou km)
% vc0: vetor velocidade em rela��o ao prim�rio escrito no referencial
% celeste (m/s ou km/s - unidades coerentes com o vetor posi��o). Deve
% ser tomado no mesmo instante da medida da posi��o.
% mu: par�metro gravitacional padr�o do corpo prim�rio (m^3/s^2 ou km^3/s^2
% - unidades coerentes com a posi��o e velocidade)
% Sa�das
% par_orb: vetor de par�metros orbitais
% a=par_orb(1): semi eixo maior da �rbita (m ou km - depende das unidades
% de entrada)
% e=par_orb(2): excentricidade da �rbita (adimensional)
% tau=par_orb(3): tempo de perigeu da �rbita (segundos)
% OMEGA=par_orb(4): ascen��o reta do nodo ascendente (rad). Dire��o da
% linha dos nodos da �rbita dada em rela��o ao eixo X no plano XY do 
% referencial celeste
% i=par_orb(5): inclina��o (rad). Inclina��o da �rbita dada em rela��o ao
% plano XY do referencial celeste.
% omega=par_orb(6): argumento de perigeu (rad). Em rela��o � linha dos
% nodos, medido no plano orbital.
%% C�lculos
% Dist�ncia radial ao prim�rio no instante observado
r0=norm(rc0);
% Vetor quantidade de movimento angular espec�fica no referencial celeste
hc=cross(rc0,vc0);
% Vetor excentricidade no sistema celeste
ec=cross(vc0,hc)/mu-rc0/r0;
% Excentricidade da �rbita
e=norm(ec);
% M�dulo do vetor hc
h=norm(hc);
% Par�metro da �rbita dada
p=h^2/mu;
% Semi eixo maior
a=p/(1-e^2);
% Vetor par�metro no referencial celeste
pc=p*cross(hc,ec)/(h*e);
% Anomalia verdadeira
costheta=(p-r0)/(e*r0);
sintheta=dot(rc0,pc)/(r0*p);
theta=atan2(sintheta,costheta);
% O tempo de perigeu depende do tipo de �rbita
if (0<=e)&&(e<1)
    tipo='e';   % �rbita el�ptica
elseif e==1
    tipo='p';   % �rbita parab�lica
else
    tipo='h';   % �rbita hiperb�lica
end

% Tempo de perigeu
if  tipo=='e'   % �rbita el�ptica
    % Movimento m�dio
    n=sqrt(mu/a^3);
    % Anomalia exc�ntrica
    E=2*atan(sqrt((1-e)/(1+e))*tan(theta/2));
    tau=t0-(E-e*sin(E))/n;
elseif  tipo=='p'   % �rbita parab�lica
    tau=-((tan(theta/2))^3+3*tan(theta/2))/(mu/p^3)^(1/6);
else % �rbita hiperb�lica
    % Movimento m�dio hiperb�lico
    n=sqrt(-mu/a^3);
    % Anomalia hiperb�lica
    H=2*atanh(sqrt((e-1)/(1+e))*tan(theta/2));
    tau=-(e*sinh(H)-H)/n;
end

% Linha dos nodos
% Vetor unit�rio ao longo do vetor h (no sistema celeste)
ih=hc/h;
% Vetor unit�rio ao longo da linha dos nodos (no sistema celeste)
Kc=[0;0;1];
nc=cross(Kc,ih)/norm(cross(Kc,ih));
% Ascen��o reta do nodo ascendente
OMEGA=atan2(nc(2),nc(1));
% Inclina��o
i=acos(dot(ih,Kc));
% Vetor unit�rio ao longo do vetor excentricidade (no referencial
% celeste)
ie=ec/e;
% Argumento de perigeu
cosomega=dot(ie,nc);
sinomega=dot(ih,cross(nc,ie));
omega=atan2(sinomega,cosomega);
%% Vetor de par�metros de sa�da
par_orb(1)=a;
par_orb(2)=e;
par_orb(3)=tau;
par_orb(4)=OMEGA;
par_orb(5)=i;
par_orb(6)=omega;
end