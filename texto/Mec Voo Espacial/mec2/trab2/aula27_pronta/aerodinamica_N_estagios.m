function [D,fy,L]=aerodinamica_N_estagios(t,V,h,M,Kn,T,rho,R)
% Modelo de arrasto conforme a referencia 
% TEWARI, A. Atmospheric and Space Flight Dynamics: 
% Modelling and simulation with MATLAB and Simulink. Boston: Birkhauser, 2007.
% Exemplo 12.6
% Vale para foguete de 1, 2 ou 3 estagios, com carga util
% Assume-se o mesmo coeficiente de arrasto para o foguete todo e seus
% estagios. A magnitude do arrasto eh alterada em funcao da area de
% referencia de cada estagio
    global ts Sr fc
    %% Coeficiente de arrasto em funcao do numero de Mach e de Knuden
    CD=modelo_aerodinamico(V,h,M,Kn,T,R);
    % Fator de correcao do arrasto a partir de dados de tunel de vento
    CD=fc*CD;
    %% A area de referencia depende do estagio atual
    % Numero de estagios
    N=length(ts);
    switch N
        case 1
            S=area1estagio(t,ts,Sr);
        case 2
            S=area2estagios(t,ts,Sr);
        case 3
            S=area3estagios(t,ts,Sr);
    end
    %% Forcas
    D=0.5*rho*V^2*S*CD;
    fy=0;L=0;
end
%% Area de referencia de foguete de 1 estagio
function S=area1estagio(t,ts,Sr)
    if  t<=ts(1)
        % Foguete e carga util
        S=Sr(1);
    else
    	S=Sr(2);   % Carga util
    end
end
%% Area de referencia de foguete de 2 estagios
function S=area2estagios(t,ts,Sr)
    if  t<=ts(1)
        % Todos os estágios
      	S=Sr(1);
	elseif t<=ts(2)
        % Segundo estágio e carga útil
        S=Sr(2); 
    else
        % Carga útil
        S=Sr(3);
    end
end
%% Area de referencia de foguete de 3 estagios
function S=area3estagios(t,ts,Sr)
    if  t<=ts(1)
        % Todos os estágios
        S=Sr(1);
 	elseif t<=ts(2)
     	% Segundo estágio
        S=Sr(2); 
 	elseif t<=ts(3)
     	% Terceiro estágio e carga útil
     	S=Sr(3);
    else
        % Carga útil
    	S=Sr(4);
    end
end