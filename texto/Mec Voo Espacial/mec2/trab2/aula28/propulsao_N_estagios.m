function [ft,m]=propulsao_N_estagios(t)
% Função para cálculo dos parâmetros propulsivos em função do tempo
% Veículo de até 3 estágios
% Entrada
% t (s): tempo
% Saídas
% ft (N): força propulsiva
% m (kg): massa do foguete em função do tempo
% Hipoteses
% - Em cada estagio, eh assumida taxa de queima continua, ou seja, nao ha
% controle da queima e a mesma eh assumida uniforme do inicio ao fim do
% propelente
% - A tracao de cada estagio eh assumida como um pulso retangular perfeito,
% ou seja, quando acionado, o propulsor vai de tracao zero ateh a maxima,
% permanecendo nesse patamar constante. Ao fim da queima, a tracao cai
% instanteamente a zero
%% Entrada de parâmetros por variáveis globais
global ti tq ts Isp mp ms m0 g
%% As variáveis propulsivas dependem do estagio atual
    % Numero de estagios
    N=length(ts);
    switch N
        case 1
            [ft,m]=propulsor_1_estagio(t,ti,tq,ts,Isp,mp,ms,m0,g);
        case 2
            [ft,m]=propulsor_2_estagios(t,ti,tq,ts,Isp,mp,ms,m0,g);
        case 3
            [ft,m]=propulsor_3_estagios(t,ti,tq,ts,Isp,mp,ms,m0,g);
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