function [ft,m,mu,epsl]=propulsao_N_estagios(t,X)
% Função para cálculo dos parâmetros propulsivos em função do tempo
% Veículo de até 3 estágios com dupla ignicao do terceiro estagio
% Entrada
% t (s): tempo
% X: vetor de estado
% V = X(1) (m/s): módulo do vetor velocidade relativa com respeito ao
% planeta girante
% A = X(2) (rad): ângulo de azimute do vetor velocidade relativa com
% respeito ao eixo z (que aponta para o norte) do sistema uen.
% phi = X(3) (rad): ângulo de elevação do vetor velocidade relativa com
% respeito ao horizonte local (plano yz do referencial uen)
% r = X(4) (m): distância radial até o centro do planeta
% delta = X(5) (rad): latitude com respeito ao plano equatorial do planeta
% lon = X(6) (rad): longitude planetaria
% Saídas
% ft (N): força propulsiva
% m (kg): massa do foguete em função do tempo
% mu, epsl (rad): Angulos de apontamento da tubeira
% Hipoteses
% - Em cada estagio, eh assumida taxa de queima continua, ou seja, nao ha
% controle da queima e a mesma eh assumida uniforme do inicio ao fim do
% propelente
% - A tracao de cada estagio eh assumida como um pulso retangular perfeito,
% ou seja, quando acionado, o propulsor vai de tracao zero ateh a maxima,
% permanecendo nesse patamar constante. Ao fim da queima, a tracao cai
% instanteamente a zero
%% Entrada de parâmetros por variáveis globais
global ti tq ts Isp mp ms m0 g mL we Re
%% As variáveis propulsivas dependem do estagio atual
    % Numero de estagios
    N=length(ti);
    switch N
        case 1
            [ft,m]=propulsor_1_estagio(t,ti,tq,ts,Isp,mp,ms,m0,g);
        case 2
            [ft,m]=propulsor_2_estagios(t,ti,tq,ts,Isp,mp,ms,m0,g);
        case 3
            [ft,m]=propulsor_3_estagios(t,ti,tq,ts,Isp,mp,ms,m0,g);
       	case 4
            [ft,m]=propulsor_3_estagios_2ig(t,ti,tq,ts,Isp,mp,ms,m0,mL,g);
    end
    % Para altitudes acima de 200km, alinha o vetor de tracao com a
    % velocidade inercial ao inves da relativa
    % Desmenbra o vetor de estado
    V=X(1);A=X(2);phi=X(3);r=X(4);delta=X(5);
    % Altitude
    h=r-Re;
    if h<200e3;
        epsl=0;mu=0;    % Tracao alinhada com a velocidade relativa
    else
        % Vetor velocidade inercial
        [~,phii,Ai]=Vrel2Vine(V,phi,A,we,r,delta);    
        % Angulos propulsivos para que a tracao seja alinhada com a
        % velocidade inercial
        mu=asin(cos(A)*cos(phii)*sin(Ai)-sin(A)*cos(phii)*cos(Ai));
        epsl=-atan2(-cos(phi)*sin(phii)+sin(phi)*sin(A)*cos(phii)*sin(Ai)+...
                 sin(phi)*cos(A)*cos(phii)*cos(Ai),sin(phi)*sin(phii)+...
                 cos(phi)*sin(A)*cos(phii)*sin(Ai)+cos(phi)*cos(A)*cos(phii)*cos(Ai));
    end
end
%% Modelo do propulsor de foguete de 1 estágio
function [ft,m]=propulsor_1_estagio(t,ti,tq,ts,Isp,mp,ms,m0,g)
%% Cálculo da massa e tração em foguete com 1 estagio e carga util
    if t<=ti(1)
        % Antes da ignição
        m=m0;   % Massa inicial
        ft=0;   % Força propulsiva nula
    elseif t<=tq(1)
        % Taxa de queima contínua
        md=-mp(1)/(tq(1)-ti(1));
        % Está queimando o primeiro estágio
        m=m0+md*(t-ti(1));
        % Força propulsiva constante
        ft=-g*Isp(1)*md;
    elseif t<=ts(1)
        % Entre a queima e a separação
        m=m0-mp(1);
        ft=0;
    else
        % Após a separação do motor foguete
        m=m0-mp(1)-ms(1);
        ft=0;
    end
end
%% Modelo dos propulsores de foguete de 2 estágios
function  [ft,m]=propulsor_2_estagios(t,ti,tq,ts,Isp,mp,ms,m0,g)
%% Cálculo da massa e tração em foguete com 2 estagios e carga util
    if t<=ti(1)
        % Antes da ignição
        m=m0;   % Massa inicial
        ft=0;   % Força propulsiva nula
    elseif t<=tq(1)
        % Taxa de queima contínua
        md=-mp(1)/(tq(1)-ti(1));
        % Está queimando o primeiro estágio
        m=m0+md*(t-ti(1));
        % Força propulsiva constante
        ft=-g*Isp(1)*md;
    elseif t<=ts(1)
        % Entre a queima e a separação
        m=m0-mp(1);
        ft=0;
    elseif t<=ti(2)
        % Entre a separação e a ignição
        m=m0-mp(1)-ms(1);
        ft=0;
    elseif t<=tq(2)
        % Taxa de queima contínua no segundo estágio
        md=-mp(2)/(tq(2)-ti(2));
        % Durante a queima do segundo estágio
        m02=m0-mp(1)-ms(1);
        m=m02+md*(t-ti(2));
        % Força propulsiva constante
        ft=-g*Isp(2)*md;
    elseif t<=ts(2)
        % Após a queima do segundo estágio e antes da separação do mesmo
        m=m0-mp(1)-ms(1)-mp(2);
        ft=0;
    else
        % Após a separação do segundo estágio
        m=m0-mp(1)-ms(1)-mp(2)-ms(2);
        ft=0;
    end
end
%% Modelo dos propulsores de foguete de 3 estágios
function [ft,m]=propulsor_3_estagios(t,ti,tq,ts,Isp,mp,ms,m0,g)
    %% Cálculo da massa e tração
    if t<=ti(1)
        % Antes da ignição
        m=m0;   % Massa inicial
        ft=0;   % Força propulsiva nula
    elseif t<=tq(1)
        % Taxa de queima contínua
        md=-mp(1)/(tq(1)-ti(1));
        % Está queimando o primeiro estágio
        m=m0+md*(t-ti(1));
        % Força propulsiva constante
        ft=-g*Isp(1)*md;
    elseif t<=ts(1)
        % Entre a queima e a separação
        m=m0-mp(1);
        ft=0;
    elseif t<=ti(2)
        % Entre a separação e a ignição
        m=m0-mp(1)-ms(1);
        ft=0;
    elseif t<=tq(2)
        % Taxa de queima contínua no segundo estágio
        md=-mp(2)/(tq(2)-ti(2));
        % Durante a queima do segundo estágio
        m02=m0-mp(1)-ms(1);
        m=m02+md*(t-ti(2));
        % Força propulsiva constante
        ft=-g*Isp(2)*md;
    elseif t<=ts(2)
        % Após a queima do segundo estágio e antes da separação do mesmo
        m=m0-mp(1)-ms(1)-mp(2);
        ft=0;
    elseif t<=ti(3)
        % Entre a separação e a ignição
        m=m0-mp(1)-ms(1)-mp(2)-ms(2);
        ft=0;
    elseif t<=tq(3)
        % Taxa de queima contínua no terceiro estágio
        md=-mp(3)/(tq(3)-ti(3));
        % Durante a queima do terceiro estágio
        m03=m0-mp(1)-ms(1)-mp(2)-ms(2);
        m=m03+md*(t-ti(3));
        % Força propulsiva constante
        ft=-g*Isp(3)*md;
    elseif t<=ts(3)
        % Após a queima do terceiro estágio e antes da separação do mesmo
        m=m0-mp(1)-ms(1)-mp(2)-ms(2)-mp(3);
        ft=0;  
    else
        % Após a separação do terceiro estágio
        m=m0-mp(1)-ms(1)-mp(2)-ms(2)-mp(3)-ms(3);
        ft=0;
    end
end

%% Modelo dos propulsores de foguete de 3 estágios com dupla ignicao do terceiro
function [ft,m]=propulsor_3_estagios_2ig(t,ti,tq,ts,Isp,mp,ms,m0,mL,g)
    %% Cálculo da massa e tração
    if t<=ti(1)
        % Antes da ignição
        m=m0;   % Massa inicial
        ft=0;   % Força propulsiva nula
    elseif t<=tq(1)
        % Taxa de queima contínua
        md=-mp(1)/(tq(1)-ti(1));
        % Está queimando o primeiro estágio
        m=m0+md*(t-ti(1));
        % Força propulsiva constante
        ft=-g*Isp(1)*md;
    elseif t<=ts(1)
        % Entre a queima e a separação
        m=m0-mp(1);
        ft=0;
    elseif t<=ti(2)
        % Entre a separação e a ignição
        m=m0-mp(1)-ms(1);
        ft=0;
    elseif t<=tq(2)
        % Taxa de queima contínua no segundo estágio
        md=-mp(2)/(tq(2)-ti(2));
        % Durante a queima do segundo estágio
        m02=m0-mp(1)-ms(1);
        m=m02+md*(t-ti(2));
        % Força propulsiva constante
        ft=-g*Isp(2)*md;
    elseif t<=ts(2)
        % Após a queima do segundo estágio e antes da separação do mesmo
        m=m0-mp(1)-ms(1)-mp(2);
        ft=0;
    elseif t<=ti(3)
        % Entre a separação e a ignição
        m=m0-mp(1)-ms(1)-mp(2)-ms(2);
        ft=0;
    elseif t<=tq(3)
        % Taxa de queima contínua no terceiro estágio - primeira ignicao
        md=-mp(3)/(tq(3)-ti(3));
        % Durante a queima do terceiro estágio
        m03=m0-mp(1)-ms(1)-mp(2)-ms(2);
        m=m03+md*(t-ti(3));
        % Força propulsiva constante
        ft=-g*Isp(3)*md;
    elseif t<=ti(4)
        % Antes da nova queima do terceiro estagio
        m=m0-mp(1)-ms(1)-mp(2)-ms(2)-mp(3);
        ft=0;
 	elseif t<=tq(4)
        % Taxa de queima contínua no terceiro estágio - segunda ignicao
        md=-mp(4)/(tq(4)-ti(4));
        % Durante a queima do terceiro estágio
        m03=m0-mp(1)-ms(1)-mp(2)-ms(2)-mp(3);
        m=m03+md*(t-ti(4));
        % Força propulsiva constante
        ft=-g*Isp(3)*md;
    elseif t<=ts(3)
        % Após a queima do terceiro estágio e antes da separação do mesmo
        m=m0-mp(1)-ms(1)-mp(2)-ms(2)-mp(3)-mp(4);
        ft=0;
    else
        % Após a separação do terceiro estágio
        m=mL;
        ft=0;
    end
end