clear vars;
close all;
clc;
s = tf('s');


Gs = 80*exp(-0.1*s)/(s*(s+4)*(s+10));
figure, nichols(Gs);

Gs = feedback(Gs,1);
figure,bode(Gs);
%% Tabela de dados Curva de Nichols

Gs = 1/(s*(s+1)*(s*0.2+1));
figure, nichols(Gs);

figure,bode(Gs);
%%
Gs = 0.64/(s*(s^2 +s+1));
figure, nichols(Gs);
figure,bode(Gs);
pole(Gs)