function [vi,phii,Ai]=Vrel2Vine(vr,phir,Ar,we,r,dt)
% Função para converter a velocidade, elevação e azimute da velocidade
% relativa para a velocidade inercial
% Entradas
% vr (m/s): velocidade relativa
% phir (rad): inclinação da velocidade relativa
% Ar (rad): azimute da velocidade relativa
% we (rad/s): velocidade de rotação do referencial girante
% r (m): distância radial até a origem do referencial inercial
% dt (rad): latitude
% Saídas
% v (m/s): magnitude da velocidade com respeito ao referencial inercial
% phi (rad): ângulo de elevação da velocidade inercial (angulo de
% trajetoria)
% A (rad): ângulo de azimute da velocidade inercial
%% Cálculos
Ai=atan2(vr*cos(phir)*sin(Ar)+we*r*cos(dt),vr*cos(phir)*cos(Ar));
if Ai<0
    Ai=Ai+2*pi;
end
vi=sqrt(vr^2+2*vr*cos(phir)*sin(Ar)*r*we*cos(dt)+r^2*we^2*cos(dt)^2);
phii=atan2(sin(phir)*cos(Ai),cos(phir)*cos(Ar));
end