% Aula 13
% Exemplo 1 - Avanco de fase
% Exercicio sala de aula
%
close all
clear all
clc

% Sistema sem compensacao
np=10;
dp=conv([1 0 0],[1 10]);
Gnc = tf(np,dp,'InputDelay',0.05)

figure(1)
subplot(121);
margin(Gnc);
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
grid;

subplot(122);
step(feedback(Gnc,1));
h = findobj(gcf,'type','line');
set(h,'linewidth',2);
grid;

% Calcula sobresinal em percentagem e tempo de acomodacao
S1 = stepinfo(feedback(Gnc,1),'RiseTimeLimits',[0.02,0.98]);
Mp1=S1.Overshoot;
ts1=S1.SettlingTime;
title(['Mp=', num2str(Mp1), '   ts=', num2str(ts1)]);
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])

%% Sistema com ganho ajustado mas nao compensado
K=1;
%np=conv([1 0],[0 1])
dp=conv([1 0 0],[1 10]);
%%
Gnk = tf(np,dp,'InputDelay',0.05)
% Diagrama de Bode do sistema com ganho ajustado
figure(2)
subplot(121);
margin(Gnk);
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])
grid;

% Resposta ao degrau do sistema com ganho ajustado mas nao compensado
ma_sc=Gnk; % FT de MA com ganho ajustado
[Gm,Pm,Wg,Wp] = margin(ma_sc); % Calculo da MFase e MG
mf_sc=feedback(ma_sc,1); % FT de MF do sistema com ganho ajustado
                         % mas nao compensado
subplot(122)
step(mf_sc);
h = findobj(gcf,'type','line');
set(h,'linewidth',2);
grid;
% Calcula sobresinal em percentagem e tempo de acomodaï¿½ï¿½o
S2 = stepinfo(mf_sc,'RiseTimeLimits',[0.02,0.98]);
Mp_2=S2.Overshoot;
MpdB=20*log10(Mp_2);
ts2=S2.SettlingTime;
title(['Mp(dB) =', num2str(MpdB), '   ts=', num2str(ts2)]);

figure(3)
margin(Gnk);
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])
grid;
figure(3), nichols(Gnk,{0.01,10^5});

%% Projeto do compensador SEM TOLERANCIA 
tole=0;   % Tolerancia
MFase=10; % MFase desejada
Pm=-8.55;    % MFase que tem o sistema

phi_m=MFase-Pm+tole; % Maximo angulo de avanco

alpha=(1-sin(phi_m*pi/180))/(1+sin(phi_m*pi/180)); % Calcula o alpha

g = 10*log10(alpha);
disp(g)
%%
wm=1.18;
% Calculo do T
T=1/(wm*sqrt(alpha));

nc=[T 1]
dc=[alpha*T 1]
% nc=[0.6199 1]
% dc=[0.0365 1]

ma_sc;%sys_p=tf(n1,d1); % FT de MA do sistema com ganho ajustado
sys_c=tf(nc,dc) % FT do compensador

sys_ma_c=tf(series(ma_sc,sys_c)); % FT de MA do sistema compensado

figure(5)
margin(sys_ma_c); % Diagrama de Bode da FT de MA do sistema compensado
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])
grid;

%% Projeto do compensador COM TOLERANCIA 
tolet=0;   % Tolerancia
MFaset=10; % MFase desejada
Pmt=-8.55;    % MFase que tem o sistema

phi_mt=MFaset-Pmt+tolet; % Maximo angulo de avanco

alphat=(1-sin(phi_mt*pi/180))/(1+sin(phi_mt*pi/180)); % Calcula o alpha

g = 10*log10(alphat);

disp(g)
wmt=1.18;  %  tentiva 6 - tolerancia 14.0°

%%
% Calculo do T
T=1/(wmt*sqrt(alphat));

%nct=[T 1]
%dct=[alphat*T 1]
nct=[16 16*0.7]
dct=[1 9]

%Gnk %sys_p=tf(n1,d1); % FT de MA do sistema com ganho ajustado
sys_ct=tf(nct,dct) % FT do compensador

%sys_ma_ct=tf(series(ma_sc,sys_ct)); % FT de MA do sistema compensado
sys_ma_ct=tf(series(Gnc,sys_ct)); % FT de MA do sistema compensado

figure(7)

margin(sys_ma_ct); % Diagrama de Bode da FT de MA do sistema compensado
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])
grid;

%figure,nichols(sys_ma_ct);
% Valores desejados
%csi =  sqrt((1-(1.4/10)^2)/2)
%Mr = 1/(2*csi*sqrt(1-csi^2))


%
sys_mf_ct=feedback(sys_ma_ct,1) % FT de MF do sistema compesado
nichols(sys_mf_ct)
% Resposta ao degrau unitario
figure(8)
hold on;
%step(feedback(tf(np,dp),1),'b');
%step(mf_sc,'r');
step(sys_mf_ct,'k');
h = findobj(gcf,'type','line');
set(h,'linewidth',2);
hold off;
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])
legend('Sistema','G.A.','C.C.')
grid;
S3 = stepinfo(sys_mf_ct,'RiseTimeLimits',[0.02,0.98]);
Mp3=S3.Overshoot;
ts3=S3.SettlingTime;
title(['M_p(%)=   ',num2str(Mp3),'    t_s=  ',num2str(ts3)]);
%%
figure(9);
hold on;
bode(mf_sc);
bode(sys_mf_ct);
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])
grid;
legend('Sem controle','Com controle');
hold off;

figure(10);
t=0:0.1:5;
u=t;
ysc=lsim(feedback(tf(np,dp),1),u,t);
yga=lsim(mf_sc,u,t);
ycc=lsim(sys_mf_ct,u,t);
plot(t,u,'y',t,ysc,'b',t,yga,'k',t,ycc,'r');
h = findobj(gcf,'type','line');
set(h,'linewidth',2);
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])
legend('Rampa','Sistema','G.A.','C.C.','NorthEast')
grid;

figure(11);
nichols(sys_ma_ct);
hold on;
nichols(ma_sc/100);
