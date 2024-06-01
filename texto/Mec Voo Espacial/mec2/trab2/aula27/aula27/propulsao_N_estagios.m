function [ft,m]=propulsao_N_estagios(t)
% Fun��o para c�lculo dos par�metros propulsivos em fun��o do tempo
% Ve�culo de at� 3 est�gios
% Entrada
% t (s): tempo
% Sa�das
% ft (N): for�a propulsiva
% m (kg): massa do foguete em fun��o do tempo
% Hipoteses
% - Em cada estagio, eh assumida taxa de queima continua, ou seja, nao ha
% controle da queima e a mesma eh assumida uniforme do inicio ao fim do
% propelente
% - A tracao de cada estagio eh assumida como um pulso retangular perfeito,
% ou seja, quando acionado, o propulsor vai de tracao zero ateh a maxima,
% permanecendo nesse patamar constante. Ao fim da queima, a tracao cai
% instanteamente a zero
%% Entrada de par�metros por vari�veis globais
global ti tq ts Isp mp ms m0 g
%% As vari�veis propulsivas dependem do estagio atual
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
%% Modelo do propulsor de foguete de 1 est�gio
function [ft,m]=propulsor_1_estagio(t,ti,tq,ts,Isp,mp,ms,m0,g)
%% C�lculo da massa e tra��o em foguete com 1 estagio e carga util
    if t<=ti(1)
        % Antes da igni��o
        m=m0;   % Massa inicial
        ft=0;   % For�a propulsiva nula
    elseif t<=tq(1)
        % Taxa de queima cont�nua
        md=-mp(1)/(tq(1)-ti(1));
        % Est� queimando o primeiro est�gio
        m=m0+md*(t-ti(1));
        % For�a propulsiva constante
        ft=-g*Isp(1)*md;
    elseif t<=ts(1)
        % Entre a queima e a separa��o
        m=m0-mp(1);
        ft=0;
    else
        % Ap�s a separa��o do motor foguete
        m=m0-mp(1)-ms(1);
        ft=0;
    end
end
%% Modelo dos propulsores de foguete de 2 est�gios
function  [ft,m]=propulsor_2_estagios(t,ti,tq,ts,Isp,mp,ms,m0,g)
%% C�lculo da massa e tra��o em foguete com 2 estagios e carga util
    if t<=ti(1)
        % Antes da igni��o
        m=m0;   % Massa inicial
        ft=0;   % For�a propulsiva nula
    elseif t<=tq(1)
        % Taxa de queima cont�nua
        md=-mp(1)/(tq(1)-ti(1));
        % Est� queimando o primeiro est�gio
        m=m0+md*(t-ti(1));
        % For�a propulsiva constante
        ft=-g*Isp(1)*md;
    elseif t<=ts(1)
        % Entre a queima e a separa��o
        m=m0-mp(1);
        ft=0;
    elseif t<=ti(2)
        % Entre a separa��o e a igni��o
        m=m0-mp(1)-ms(1);
        ft=0;
    elseif t<=tq(2)
        % Taxa de queima cont�nua no segundo est�gio
        md=-mp(2)/(tq(2)-ti(2));
        % Durante a queima do segundo est�gio
        m02=m0-mp(1)-ms(1);
        m=m02+md*(t-ti(2));
        % For�a propulsiva constante
        ft=-g*Isp(2)*md;
    elseif t<=ts(2)
        % Ap�s a queima do segundo est�gio e antes da separa��o do mesmo
        m=m0-mp(1)-ms(1)-mp(2);
        ft=0;
    else
        % Ap�s a separa��o do segundo est�gio
        m=m0-mp(1)-ms(1)-mp(2)-ms(2);
        ft=0;
    end
end
%% Modelo dos propulsores de foguete de 3 est�gios
function [ft,m]=propulsor_3_estagios(t,ti,tq,ts,Isp,mp,ms,m0,g)
    %% C�lculo da massa e tra��o
    if t<=ti(1)
        % Antes da igni��o
        m=m0;   % Massa inicial
        ft=0;   % For�a propulsiva nula
    elseif t<=tq(1)
        % Taxa de queima cont�nua
        md=-mp(1)/(tq(1)-ti(1));
        % Est� queimando o primeiro est�gio
        m=m0+md*(t-ti(1));
        % For�a propulsiva constante
        ft=-g*Isp(1)*md;
    elseif t<=ts(1)
        % Entre a queima e a separa��o
        m=m0-mp(1);
        ft=0;
    elseif t<=ti(2)
        % Entre a separa��o e a igni��o
        m=m0-mp(1)-ms(1);
        ft=0;
    elseif t<=tq(2)
        % Taxa de queima cont�nua no segundo est�gio
        md=-mp(2)/(tq(2)-ti(2));
        % Durante a queima do segundo est�gio
        m02=m0-mp(1)-ms(1);
        m=m02+md*(t-ti(2));
        % For�a propulsiva constante
        ft=-g*Isp(2)*md;
    elseif t<=ts(2)
        % Ap�s a queima do segundo est�gio e antes da separa��o do mesmo
        m=m0-mp(1)-ms(1)-mp(2);
        ft=0;
    elseif t<=ti(3)
        % Entre a separa��o e a igni��o
        m=m0-mp(1)-ms(1)-mp(2)-ms(2);
        ft=0;
    elseif t<=tq(3)
        % Taxa de queima cont�nua no terceiro est�gio
        md=-mp(3)/(tq(3)-ti(3));
        % Durante a queima do terceiro est�gio
        m03=m0-mp(1)-ms(1)-mp(2)-ms(2);
        m=m03+md*(t-ti(3));
        % For�a propulsiva constante
        ft=-g*Isp(3)*md;
    elseif t<=ts(3)
        % Ap�s a queima do terceiro est�gio e antes da separa��o do mesmo
        m=m0-mp(1)-ms(1)-mp(2)-ms(2)-mp(3);
        ft=0;  
    else
        % Ap�s a separa��o do terceiro est�gio
        m=m0-mp(1)-ms(1)-mp(2)-ms(2)-mp(3)-ms(3);
        ft=0;
    end
end