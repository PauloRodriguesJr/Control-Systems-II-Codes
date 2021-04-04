%Lista 2 - Questão 6
% setup inicial
close all;
clearvars;
clc;

% Sistema sem compensacao
np=1;
dp=conv([1 0],conv([1 10],[1 14]));

%csi_d= 0.707

figure(1)  % Margens do sistema sem compensação
subplot(121);
margin(np,dp);
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
grid;

subplot(122); % Resposta degrau sistema sem compensação
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
% Calculando o ganho K, obtém-se K=14 (para atender resposta a rampa
K=14;
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
% Calcula sobresinal em porcentagem e tempo de acomodacao
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

% Resposta à Rampa do sistema ajustado
figure(4);
t=0:0.1:100;
u=t;
ysc=lsim(feedback(ma_sc,1),u,t);
y=lsim(feedback(tf(np,dp),1),u,t);
plot(t,u,'y',t,ysc,'b',t,y,'r');
h = findobj(gcf,'type','line');
set(h,'linewidth',2);
%set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])
ylabel('Amplitude');
xlabel('Tempo (s)');
legend('Rampa','Sist. ajustado ','Sist.original','Location','NorthWest')
grid;
title('Resposta à rampa para o sistema original e compensado')

%% Determinação da margem de fase:
csi_d = 0.707;
MF_d= atan2(2*csi_d, sqrt(sqrt(1+4*csi_d^4)-2*csi_d^2))*180/pi;
fprintf('A margem desejada é de %.3f\n', MF_d);
tole= 12;
phi_m = -180+(MF_d+tole);
fprintf('A fase de interesse é %.2f\n',phi_m);
% Buscar na tabela e diagrama de Bode os valores de interesse:

figure(5)
margin(n1,d1);
set(findall(gcf,'Type','text'),'FontSize',14);
set(findall(gcf,'type','line'),'linewidth',2);
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25]);
grid;
%% Ajustar os intervalos para buscar o valor mais adequado
% args: (sysd,wi,wf,step)
table_write_tf(tf(n1,d1),1.2,1.3,0.01); %(se necessário)
%%
clc;
wcg = 1.21;
% Conversão dB para linear
beta =  21.8 %dB;% É necessário em esccala linear
% Conversão dB para linear
beta = 10^(beta/20)

% O intervalo de liberdade para a escolha de w_c está contido entre
% uma década e uma oitava abaixo de wcg:
wcmax = wcg/2;
wcmin = wcg/10;

fprintf('Escolher um wc entre %.3f --- %.3f',wcmax,wcmin);
wc = 0.5; % Escolha o valor!

%% Parametros do controlador:
T = 1/wc;

nc=[T 1];
dc=[beta*T 1];

G_c = tf(nc,dc) % FT do compensador

bode(G_c);
sys_ma_ct=tf(series(G_c,tf(n1,d1))); % FT de MA do sistema compensado
sys_mf_ct = feedback(sys_ma_ct,1); % fechando a malha

figure(6)
margin(sys_ma_ct); % Diagrama de Bode da FT de MA do sistema compensado
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])
grid;

%% Parametros temporais da do sistema compensado 
%% Resposta ao degrau unitario
figure(7)
t = 0:0.01:1000;
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

%% Resposta à Rampa do sistema
figure(8);
t=0:0.1:250;
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

%% Comparação dos diagramas de Bode
figure(9);
hold on;
bode(feedback(tf(np,dp),1));
bode(sys_mf_ct);
set(findall(gcf,'Type','text'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])
hold on
legend('Sem controle','Com controle');
title('Comparação do Diagrama de Bode: Sistema Original x Sistema Compensado em Malha Fechada')

grid;   
hold off;


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
%%

%% ***** Plot o sinal de controle *****
t_2 = [0:0.01:120];            % vetor de tempo 
U_s = G_c/(1+G_c*tf(n1,d1));        % Formulação da expr. do sinal de controle
y_s = step(U_s,t_2);    % Resposta ao degrau 1/s
figure, plot(t_2,y_s,'b');    % Plot do sinal de controle
title('Evolução do sinal de controle no sistema');
grid;
xlabel('Tempo (s)');
ylabel('Esforço requerido');


%% Fechar figuras e ver com calma
%close all;

figure, rlocus(sys_mf_ct);

figure, rlocus(mf_sc);
