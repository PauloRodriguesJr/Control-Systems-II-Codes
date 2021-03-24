% Prova 1 Sistemas de Controle II - Magno Meza
% Setup para os cálculos e limpeza do ambiente
clc;
clear all;
close all;
s= tf('s');

% Questão 1 - Diagrama de Bode

Gs = (2*10^5*s+2*10^6)/(s^4+40*s^3+10400*s^2+2*10^5*s)

disp(pole(Gs));

bode(Gs)

%% Questão 2 - Diagrama polar
Gs = (8*s+16)/(s^2+3*s+16)

pole(Gs)
nyquist(Gs)

%% Questão 3 - Diagrama fase não minima
Gs = 1/((s-1)*(s+2)*(s+3))

pole(Gs)
nyquist(Gs)

%% Questão 4 - Nichols Chart

Gs = 1/(s*(s+1)*(s*0.1+1));
nichols(Gs)
margin(Gs)


%Plot experimental
Gs = -0.3*(-s+1)*(s/50+1)/(s*(s/400+1)*(s/10^4+1));
w = logspace(-2,6,1000);
bode(Gs,w)
grid on;
title('Diagrama de Bode para a função de transferência identificada')
%% Questão 6 - FT
close all;
%
K=1;
K=2847.5;
K=2850;
K=10^(69.11/20);
K=10^(79.11/20)
K=1/1845;
K=2850;
K=10^(84.11/20)
K=10000;
Gs= K*100*(s+5)*(s+40)/(s^3*(s+100)*(s+200));
figure,nyquist(Gs);
pole(Gs)
figure, margin(Gs)
