% Laboratório 1 -Questão 1a
close all;
clear all;
clc;

% Sistema sem compensacao
np=1;
dp=conv([1 0],conv([0.04 1],[0.001 1]));

figure(1)
subplot(121);
margin(np,dp);
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
grid;

subplot(122);
step(feedback(tf(np,dp),1));
h = findobj(gcf,'type','line');
set(h,'linewidth',2);
grid;

% Calcula sobresinal em percentagem e tempo de acomodacao
S1 = stepinfo(feedback(tf(np,dp),1),'RiseTimeLimits',[0.02,0.98]);
Mp1=S1.Overshoot;
ts1=S1.SettlingTime;
title(['Mp=', num2str(Mp1), '   ts=', num2str(ts1)]);
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])

%% Sistema com ganho ajustado mas nao compensado
K=100;
n1=K*np;
d1=dp;

% Diagrama de Bode do sistema com ganho ajustado
figure(2)
subplot(121);
margin(n1,d1);
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])
grid;

% Resposta ao degrau do sistema com ganho ajustado mas nao compensado
ma_sc=tf(n1,d1); % FT de MA com ganho ajustado
[Gm,Pm,Wg,Wp] = margin(ma_sc); % Calculo da MFase e MG
mf_sc=feedback(ma_sc,1); % FT de MF do sistema com ganho ajustado
                         % mas nao compensado
subplot(122)
step(mf_sc);
h = findobj(gcf,'type','line');
set(h,'linewidth',2);
grid;
% Calcula sobresinal em percentagem e tempo de acomodacao
S2 = stepinfo(mf_sc,'RiseTimeLimits',[0.02,0.98]);
sobresinal2=S2.Overshoot;
ts2=S2.SettlingTime;
title(['Mp=', num2str(sobresinal2), '   ts=', num2str(ts2)]);

figure(3)
margin(n1,d1);
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])
grid;

%% Projeto do compensador SEM TOLERANCIA 
tole = 0;   % Tolerancia
MFase = 65; % MFase desejada
Pm = 25.3;    % MFase que tem o sistema

phi_m=MFase-Pm+tole; % Maximo angulo de avanco

alpha=(1-sin(phi_m*pi/180))/(1+sin(phi_m*pi/180)); % Calcula o alpha

wm=71.15;
%wm=71.15;  T = 0;

% Calculo do T
T=1/(wm*sqrt(alpha));

nc=[T 1];
dc=[alpha*T 1];

sys_c=tf(nc,dc); % FT do compensador

sys_ma_c=tf(series(ma_sc,sys_c)); % FT de MA do sistema compensado

%% Figura > Bode MA do sistema compensado
figure(5)
margin(sys_ma_c); % Diagrama de Bode da FT de MA do sistema compensado
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])
grid;

%% Projeto do compensador COM TOLERANCIA 
tolet=14.0;   % Tolerancia
MFaset=65; % MFase desejada
Pmt=25.3;    % MFase que tem o sistema

phi_mt=MFaset-Pmt+tolet; % Maximo angulo de avanco

alphat=(1-sin(phi_mt*pi/180))/(1+sin(phi_mt*pi/180)); % Calcula o alpha

% Ganho a buscar a respectiva frequencia no Bode ajustado
ganho_objetivo = 10*log10(alphat);

disp(ganho_objetivo)
%wmt=80.54; %  tentiva 2 - tolerancia 10°
%wmt=82.87; %  tentiva 3 - tolerancia 12°
%wmt=86.73;  %  tentiva 4 - tolerancia 15°
%wmt=86.05;  %  tentiva 5 - tolerancia 14.5°
%wmt=85.79;  %  tentiva 6 - tolerancia 14.3°
wmt=85.39;  %  tentiva 6 - tolerancia 14.0°

%% Calculo do controlador com tolerância

T=1/(wmt*sqrt(alphat));
nct=[T 1]
dct=[alphat*T 1]

sys_ct=tf(nct,dct) % FT do compensador

sys_ma_ct=tf(series(ma_sc,sys_ct)); % FT de MA do sistema compensado

figure(7)
margin(sys_ma_ct); % Diagrama de Bode da FT de MA do sistema compensado
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])
grid;
%save_fig_pdf('1a_bode_compensado_ma',gcf,gca);

%% Parametros temporais da do sistema compensado 
%% Resposta ao degrau unitario
figure(8)
t = 0:0.01:10;
hold on;
step(feedback(tf(np,dp),1),t,'b');
step(sys_mf_ct,t,'r');
h = findobj(gcf,'type','line');
set(h,'linewidth',2);
%xlim([0 8]);
%ylim([0 8]);
hold off;
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
legend('M.Fechada sem controle','M.Fechada com controle');
grid;

S3 = stepinfo(sys_mf_ct,'RiseTimeLimits',[0.02,0.98]);
Mp3= S3.Overshoot;
%Mp3db = 20*log10(Mp3);
ts3=S3.SettlingTime;
title(['Resposta ao degrau :   ' , 'M_p=   ',num2str(Mp3),' %','    t_s=  ',num2str(ts3),' s']);
%save_fig_pdf('1a_resposta_degrau',gcf,gca)

%% Comparação dos diagramas de Bode
figure(9);
hold on;
bode(tf(np,dp));
bode(sys_mf_ct);
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])
hold on
legend('Sem controle','Com controle');
title('Comparação do Diagrama de Bode: Sistema Original x Sistema Compensado em Malha Fechada')

grid;   
hold off;
%save_fig_pdf('1a_bode_compara',gcf,gca)
%% Resposta à Rampa do sistema
figure(10);
t=0:0.1:25;
u=t;
ysc=lsim(feedback(tf(np,dp),1),u,t);
ycc=lsim(sys_mf_ct,u,t);
plot(t,u,'y',t,ysc,'b',t,ycc,'r');
h = findobj(gcf,'type','line');
set(h,'linewidth',2);
%set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])
ylabel('Amplitude');
xlabel('Tempo (s)');
legend('Rampa','Sist. original','Sist.compensado','Location','NorthWest')
grid;
title('Resposta à rampa para o sistema original e compensado')
%save_fig_pdf('1a_resposta_rampa',gcf,gca) %% Arrumar enquadramento!!

%% Cartas de Nichols - Comparação
%% Gráfico comparativo Cartas de Nichols (resultado final)
figure(11);
P = nicholsoptions; 
%P.XLim = [-270 -90];
%P.YLim = [-60 80];
hold on;
nicholsplot(tf(np,dp),P);
nicholsplot(sys_ma_ct,P);
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
hold on;
grid on;
hold off;
legend('Sem controle','Com controle','Location','NorthWest');

%save_fig_pdf('1a_nichols_compara',gcf,gca);

%close all;