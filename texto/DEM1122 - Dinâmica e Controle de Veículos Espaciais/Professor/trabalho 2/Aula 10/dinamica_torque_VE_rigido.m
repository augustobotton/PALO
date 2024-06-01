function xp=dinamica_torque_VE_rigido(t,x)
% Vetor de estado: x=[wx;wy;wz;theta1;theta2;theta3]
% WXx, WYy, WZ: velocidades de rotação em torno dos eixos b1, b? e b3 do
% sistema de referencia do corpo
% theta1, theta2, theta3: angulos de Euler na sequencia de rotacoes C1, C2 e C3,
% respectivamente. Sequencia 123
%% Passagem de parametros por variaveis globais
global Ixx Iyy Izz

%% Desmembra o vetor de estado
wx=x(1);wy=x(2);wz=x(3);
theta1=x(4);theta2=x(5);theta3=x(6);
%% Vetor de controle
[Mx,My,Mz]=controle_spin_up(t);
%% Equacoes de dinamica
wxp=(Iyy-Izz)*wz*wy/Ixx + Mx/Ixx;
wyp=(Izz-Ixx)*wz*wx/Iyy+ My/Iyy;
wzp=(Ixx-Iyy)*wx*wy/Izz+ Mz/Izz;
%% Equacoes de cinematica
theta1p=(wx*cos(theta3)-wy*sin(theta3))/cos(theta2);
theta2p=wx*sin(theta3)+wy*cos(theta3);
theta3p=(-wx*cos(theta3)*sin(theta2)+wy*sin(theta3)*sin(theta2))/cos(theta2)+wz;
%% Derivada do vetor de estado
xp=[wxp;wyp;wzp;theta1p;theta2p;theta3p];
end

%% Vetor de controle da manobra de spin up
function [Mx,My,Mz]=controle_spin_up(t)
global MX MY MZ T
    % O controle eh malha aberta. Aplica-se um momento constante desde t=0 ateh
    %t=T
    if t<=T
        Mx=MX;My=MY;Mz=MZ;
    else
        Mx=0;My=0;Mz=0;
    end
end
