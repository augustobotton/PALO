function exemplo_8p4
% Exemplo 8.4 do livro - Calculo das derivadas de carga util com respeito
% as massas estruturais e de propelente.
% TEWARI, A. Atmospheric and Space Flight Dynamics: 
% Modelling and simulation with MATLAB and Simulink. Boston: Birkhauser, 2007.
% Utilizacao dos dados dos exemplos 8.2 e 8.3
%% Dados
% Massa no inicio da queima de cada estagio
m0=[1.094266769684211e+05 5.284772960000001e+04 19806.1060 6046.9099];
% Massa de propelente consumida em cada estagio
mp = [55000 2.8978709948e+04 1.279605244274179e+04 4568.81457];
% Carga util - kg
mL = 1.124811598501497e+03;
% Velocidade de exaustao de cada estagio (estagios zero e 1 sao
% equivalentes)
ve = 1.0e+03 *[2.363318181818182   2.844900000000000   2.844900000000000   4.463550000000001];
% Determina as massas estruturais de cada estagio
N=length(m0);
ms=zeros(1,N);
ms(N)=m0(N)-mL-mp(N);
ms(1:N-1)=m0(1:N-1)-m0(2:N)-mp(1:N-1);
% Massas ao final da queima de cada estagio
mf=m0-mp;
% Mostra as massas
disp('Massa da carga util (kg)');
disp(mL);
disp('Massas no inicio da queima de cada estagio');
disp(m0);
disp('Massas no final da queima de cada estagio');
disp(mf);
disp('Massas estruturais de cada estagio');
disp(ms);
disp('Massas de propelente de cada estagio');
disp(mp);
%% Calculo das derivadas de massa de carga util com respeito as massas estruturais
dmLds=der_mL_dms(mp,ms,mL,ve);
disp('Vetor de derivadas da carga util com respeito as massas estruturais: ');
disp(dmLds);
%% Calculo das derivadas de massa de carga util com respeito as massas de propelente
dmLdp=der_mL_dmp(mp,ms,mL,ve);
disp('Vetor de derivadas da carga util com respeito as massas de propelente: ');
disp(dmLdp);
end
%% Funcao para calculo das derivadas de massa de
% carga util com respeito as massas estruturais
function dmLds=der_mL_dms(mp,ms,mL,ve)
    % Numero de estagios
    N=length(mp);
    % Vetores de massa final e inicial por estagio
    mf=zeros(1,N);
    m0=zeros(1,N);
    for i=1:N
        m0(i)=mL+sum(mp(i:end))+sum(ms(i:end));
        mf(i)=m0(i)-mp(i);
    end
    % Derivadas da massa de carga util com respeito as massas estruturais
    dmLds=zeros(1,N);
    for k=1:N
        dmLds(k)=-sum(ve(1:k).*(1./m0(1:k)-1./mf(1:k)))/sum(ve.*(1./m0-1./mf));
    end
end
%% Funcao para calculo das derivadas de massa de
% carga util com respeito as massas de propelente
function dmLdp=der_mL_dmp(mp,ms,mL,ve)
    % Numero de estagios
    N=length(mp);
    % Vetores de massa final e inicial por estagio
    mf=zeros(1,N);
    m0=zeros(1,N);
    for i=1:N
        m0(i)=mL+sum(mp(i:end))+sum(ms(i:end));
        mf(i)=m0(i)-mp(i);
    end
    % Derivadas da massa de carga util com respeito as massas propelente
    dmLdp=zeros(1,N);
    for k=1:N
        dmLdp(k)=-(ve(k)/m0(k)+sum(ve(1:k-1).*(1./m0(1:k-1)-1./mf(1:k-1))))/sum(ve.*(1./m0-1./mf));
    end
end