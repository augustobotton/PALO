function [R,V]=RvelPolar2RvelRet(v,A,phi,r,lat,long)
% Fun��o para converter velocidade do sistema LVLH (coordenadas polares)
% para o sistema ECI ou ECEF retangular
% A velocidade pode ser a relativa ou a inercial, o resultado final ser�
% correlato
% Entradas
% v (m/s): M�dulo do vetor velocidade
% A (rad): Azimute da velocidade
% phi (rad): Eleva��o da velocidade
% r (m): dist�ncia radial
% lat (rad): Latitude 
% long (rad): Longitude no referencial desejado (ECI ou ECEF)
% Sa�das
% R=[R_X, R_Y, R_Z]^T (m): Vetor posi��o em coordenadas retangulares no
% sistema ECI ou ECEF (dependendo da entrada de dados de longitude)
% V=[V_X, V_Y, V_Z]^T (m/s): Vetor velocidade em coordenadas retangulares no
% sistema ECI ou ECEF (dependendo da entrada de dados de longitude). Pode
% ser a velocidade relativa ou a inercial, dependendo dos dados de
% velocidade fornecidos.
%% C�lculos
% Matriz de convers�o do sistema ECI ou ECEF para o LVLH 
CLH=[cos(lat)*cos(long) cos(lat)*sin(long) sin(lat)
    -sin(long) cos(long) 0
    -sin(lat)*cos(long) -sin(lat)*sin(long) cos(lat)];
% Vetor velocidade em coordenadas cartezianas no sistema LVLH
Vlvlh=v*[sin(phi)
        cos(phi)*sin(A)
        cos(phi)*cos(A)];
% Transforma��o da velocidade para o sistema ECI ou ECEF em coordenadas retangulares
V=CLH'*Vlvlh;
% Vetor posi��o no sistema LHVLH
Rlvlh=[r
       0
       0];
% Transforma��o da posi��o para o sistema ECI ou ECEF em coordenadas
% retangulares
R=CLH'*Rlvlh;
end