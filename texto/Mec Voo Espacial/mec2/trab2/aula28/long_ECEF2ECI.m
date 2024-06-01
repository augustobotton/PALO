function long_c=long_ECEF2ECI(t,long,we,tg)
% Fun��o para calcular a longitude celeste a partir da longitude
% fixa ao planeta
% Entradas
% t (s) - Tempo no qual se deseja saber a longitude celeste
% long (rad) - Longitude relativa ao referencial fixo ao planeta
% we (rad/s) - Velocidade de rota��o do planeta
% tg (s) - Tempo no qual o meridiano de refer�ncia tem longitude celeste
% nula
% Sa�da
% long_c (rad) - Longitude celeste no tempo t
long_c=long+we*(t-tg);
end