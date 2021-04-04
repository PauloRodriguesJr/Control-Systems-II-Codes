% Lab. 1 Parte b: sistema com atraso de transporte

% setup inicial
close all;
clearvars;
clc;

% Sistema sem compensacao
np=10;
dp=conv([1 0 0],[1 10]);
Gnc = tf(np,dp,'InputDelay',0.05);

% Plot das margens do sist. original  (Bode)
figure(1)
subplot(121);
margin(Gnc);
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
grid;

% Resposta ao degrau do sist. original
subplot(122);
step(Gnc);
h = findobj(gcf,'type','line');
set(h,'linewidth',2);
grid;

% Calcula sobresinal em porcentagem e tempo de acomodacao
S1 = stepinfo(Gnc,'RiseTimeLimits',[0.02,0.98]);
Mp1=S1.Overshoot;
ts1=S1.SettlingTime;
title(['Mp=', num2str(Mp1), '   ts=', num2str(ts1)]);
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])


%% Sistema com ganho ajustado mas nao compensado 
% Compensar sobressinal na ressonância

%K=0.78;
K=1;
np= 10*K;
dp=conv([1 0 0],[1 10]);
Gnk = tf(np,dp,'InputDelay',0.05);
Gnk_mf = feedback(Gnk,1);
% Diagrama de Bode do sistema com ganho ajustado
figure(2);
subplot(121);
margin(Gnk_mf);
set(findall(gcf,'Type','text'),'FontSize',14);
set(findall(gcf,'type','line'),'linewidth',2);
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25]);
grid;


% Resposta ao degrau do sistema com ganho ajustado mas nao compensado
ma_sc=Gnk; % FT de MA com ganho ajustado 
[Gm,Pm,Wg,Wp] = margin(ma_sc); % Calculo da MFase e MG

% Fecha malha
mf_sc=feedback(ma_sc,1); % FT de MF do sistema com ganho ajustado
                         % mas nao compensado
subplot(122);
step(mf_sc);
h = findobj(gcf,'type','line');
set(h,'linewidth',2);
grid;

% Calcula sobresinal em percentagem e tempo de acomodacao
S2 = stepinfo(mf_sc,'RiseTimeLimits',[0.02,0.98]);
Mp_2=S2.Overshoot;
MpdB=20*log10(Mp_2);
ts2=S2.SettlingTime;
title(['Mp(dB) =', num2str(MpdB), '   ts=', num2str(ts2)]);

% figure(3)
% margin(Gnk);
% set(findall(gcf,'Type','text'),'FontSize',14)
% set(findall(gcf,'type','line'),'linewidth',2)
% set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])
% grid;

% Carta de Nichols para o sistema ajustado, mas não compensado 
figure(11);
P = nicholsoptions; 
P.XLim = [-270 -130];
P.YLim = [-40 40];
P.legend
hold on;
nicholsplot(Gnc*K,P);
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
hold on;
grid on;
hold off;

%% Calcular csi e margem de fase objetivo

Mr_maxdB = 3.5; % em dB
Mr_max = 10^(3.5/20);
wr = 1.4;
%csi_d = 1/(2*Mr_max);
csis = roots([-1 0 1 0 -1/(2*Mr_max)^2]);
csi_d = csis(4,1);
MF_d= atan2(2*csi_d, sqrt(sqrt(1+4*csi_d^4)-2*csi_d^2))*180/pi;
% Calcula Maximo sobressinal em porcentagem
%Mp = exp(pi*(-csi_d/sqrt(1-csi_d^2))); 
%% Obtenção de parâmetros do controlador (Pt1)
tole = 13;%-18%30;
MF = -8.55; %45.3;  %-24; % Para ganho unitário
MF_d;

phi_m=(MF_d-MF+tole); % Maximo angulo de avanco

alpha=(1-sin(phi_m*pi/180))/(1+sin(phi_m*pi/180)); % Calcula o alpha

% Frequência de cruzamento w_m
ganho_obj = 10*log10(alpha);
disp('Ganho para verificar a fase na tabela');
disp(ganho_obj);

%% Buscar no diagrama de Bode de table_write_tf em "Gnk" (ganho ajustado)

figure(20);
w = logspace(-2,3,1000);
bode(Gnk,w) 
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])
grid;

%% Frequencia de interesse

wm = 2.15; % para o projeto com -18 - 45.5
wm = 1.62;
wm = 1.83;
wm = 1.92;
wm = 1.94;

% Calculo do T
T=1/(wm*sqrt(alpha));

% Parametros do compensador de avanço
nc=[T 1];
dc=[alpha*T 1];

G_c = tf(nc,dc); % FT do compensador

sys_ma_ct=tf(series(Gnk,G_c));%series(G_c,)); % FT de MA do sistema compensado

% % Diagrama de Bode da FT de MA do sistema compensado
figure(5)
margin(sys_ma_ct); 
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])
grid;

%%
sys_mf_ct=feedback(sys_ma_ct,1); % FT de MF do sistema compensado
figure(6)
margin(sys_mf_ct); % Diagrama de Bode da FT de MA do sistema compensado
set(findall(gcf,'Type','text'),'FontSize',14);
set(findall(gcf,'type','line'),'linewidth',2);
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25]);
grid;
%title('Diagrama de Bode em Malha Aberta do Sistema Compensado');
% Diagrama de Bode da FT de MA do sistema compensado
save_fig_pdf('1b_bode_compensado_ma',gcf,gca);

%% Parametros temporais da do sistema compensado 
%% Resposta ao degrau unitario
figure(8)
t = 0:0.01:10;
hold on;
%step(feedback(tf(np,dp),1),'b');
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
save_fig_pdf('1b_resposta_degrau',gcf,gca)


%% Gráfico comparativo Cartas de Nichols (resultado final)
figure(11);
P = nicholsoptions; 
P.XLim = [-270 -90];
%P.YLim = [-60 80];
hold on;
nicholsplot(Gnc,P);
nicholsplot(sys_ma_ct,P);
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
hold on;
grid on;
hold off;
legend('Sem controle','Com controle');
%save_fig_pdf('1b_nichols_compara',gcf,gca);

%% Grafico importante: Comparação dos diagramas de Bode
figure(9);
hold on;
%w = logspace(-1,4,10000);
bode(Gnc);
%/20; %10^(-(23-3.5)/20);
bode(sys_mf_ct);
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])
hold on
legend('Sem controle','Com controle');
title('Comparação do Diagrama de Bode: Sistema Original x Sistema Compensado em Malha Fechada')

grid;   
hold off;
%save_fig_pdf('1b_bode_compara',gcf,gca)
%% Resposta à Rampa do sistema
figure(10);
t=0:0.1:10;
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
save_fig_pdf('1b_resposta_rampa',gcf,gca)
close all;
