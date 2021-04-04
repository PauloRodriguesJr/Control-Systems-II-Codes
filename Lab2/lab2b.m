% Laboratório 2 - Projeto de Controlador por Avanço - Atraso

% **** Setup inicial****
close all
clear all
clc

% Sistema sem compensação
np=[2 1];
dp=conv([1 0],[1 0.1 4]);
tf(np,dp)

figure(1),subplot(131);
margin(np,dp);
grid;

subplot(132);
step(feedback(tf(np,dp),1));
grid;

% Calcula sobresinal em percentagem e tempo de acomodação
S1 = stepinfo(feedback(tf(np,dp),1),'RiseTimeLimits',[0.02,0.98]);
Mp1=S1.Overshoot;
ts1=S1.SettlingTime;

title(['M_p=', num2str(Mp1), '   t_s=', num2str(ts1)]);

subplot(133);
t=0:0.1:100;
u=t;
ysc=lsim(feedback(tf(np,dp),1),u,t);
plot(t,u,'k',t,ysc,'b');
legend('Rampa','Sistema','Location','NorthWest')
grid;

set(gcf, 'Position', get(0,'Screensize'));
%% Sistema com ganho ajustado mas não compensado
Kc=160;
n1=Kc*np;
d1=dp;

% Resposta ao degrau do sistema com ganho ajustado mas não compensado
ma_sc=tf(n1,d1);            % FT de MA com ganho ajustado
[Gm,Pm,Wg,Wp] = margin(ma_sc); % Calculo da MFase e MG
mf_sc=feedback(ma_sc,1);    % FT de MF do sistema com ganho ajustado
                            % mas não compensado
figure(2)
margin(n1,d1);
grid;
figure,
t=0:0.1:100;
u=t;
ysc=lsim(feedback(tf(np,dp),1),u,t);
plot(t,u,'k',t,ysc,'b');
legend('Rampa','Sistema','Location','NorthWest')
grid;
set(gcf, 'Position', get(0,'Screensize'));
%% Verificando o diagrama de Bode, verifiquemos a Margem de Fase atual
MFa = -1.27; %deg  % Inserir MF atual

%% Projeto do compensador da porção de atraso de fase
tolet=0;   % Tolerancia
MFaset=4; % MFase desejada para ser compensada pela porção de atraso
phi_t=-180+MFaset+tolet; % Fase para verificar w_cg/beta
fprintf('A fase a ser verificada no D.B/ Tabela é:   %.2f\n',phi_t);
%% Determinar, ou no diagrama de Bode acima, ou na Tabela
%% Qual é a w, |G(jw_cg)| | phi = phi_t

%table_write_tf(tf(n1,d1),wi, wf, step) % args:(tf(),wi, wf, step)

%% Projeto do compensador (Preencher beta e w_cg)
wcg=2.12 % Nova frequencia de cruzamento de ganho

% A atenuação fornecida pelo compensador
betadB = 0.243; % EM DB
beta   = 10^(betadB/20); % Pelo diagrama de Bode
%beta=10^(20/20); % Se pegar direto da tabela, pegar linear

% A frequência de canto do zero do compensador fica uma década abaixo
wz_at=wcg/10;

T2=1/(wz_at); % Calculo do T2

% A frequência de canto do polo do compensador fica
%wp_at = 1/(beta*T2);

nct=[T2 1];
dct=[beta*T2 1];

disp('Compensador de atraso: ');
sys_ct=tf(nct,dct) % FT do compensador
% Nomenclatura: sys_ma_atraso
sys_ma_at=tf(series(ma_sc,sys_ct)); % FT de MA do sistema compensado
sys_mf_at=feedback(sys_ma_at,1); % FT de MF do sistema compesado

%Figura: Diagrama de Bode da FT de MA do sistema compensado
figure(3)
margin(sys_ma_at); 
grid;
set(gcf, 'Position', get(0,'Screensize'));

%% Plot da Margem, step e rampa (3 em um) 
figure(4)
subplot(131);
margin(sys_ma_at);
grid;

subplot(132);
step(feedback(tf(np,dp),1),sys_mf_at);
ylim([0 8]);
xlim([0 100]);
hold on;
grid;
legend('Sistema','Atraso','Location','NorthEast')

% Calcula sobresinal em percentagem e tempo de acomodação
S1 = stepinfo(sys_mf_at,'RiseTimeLimits',[0.02,0.98]);
Mp1=S1.Overshoot;
ts1=S1.SettlingTime;
title(['Mp=', num2str(Mp1), '   ts=', num2str(ts1)]);
hold off;
subplot(133);
t=0:0.1:100;
u=t;
ysc_at=lsim(sys_mf_at,u,t);
% yga=lsim(mf_sc,u,t);
% ycc=lsim(sys_mf_ct,u,t);
plot(t,u,'k',t,ysc,'b',t,ysc_at,'r');
legend('Rampa','Sistema','Atraso','Location','NorthWest')
grid;
set(gcf, 'Position', get(0,'Screensize'));
%%[np_c_at,dp_c_at] = tfdata(sys_ma_at);

%% Projeto do compensador da porção de avanço de fase
MFd_av=50; % MFase desejada para ser compensada pela porção de avanço
MFa_av= -1.3;    % MFase que tem o sistema com o atraso
% Primeira tentativa
tole_av=0.0;   % Tolerancia
% terceira tentativa
%tole_av=29;   % Tolerancia

phi_m = MFd_av-MFa_av+tole_av;

alpha=(1-sin(phi_m*pi/180))/(1+sin(phi_m*pi/180));

% Procuramos na Tabela ou no gráfico a frequência na qual o módulo é
ganho_obj_t1 = -20*log10(1/sqrt(alpha));
fprintf('Busca na tabela: %.2f dB\n',ganho_obj_t1 )
%table_write_tf(tf(n1,d1),wi, wf, step) % args:(tf(),wi, wf, step)
margin(sys_ma_at); % Diagrama de Bode da FT de MA do sistema compensado
%% Verificar a freq. wm tal que ocorre o "Ganho_objetivo"
% Assim, wm é a nova frequência de cruzamento de ganho
wm_av=29.9; % Inserir resultado 

T1=1/(wm_av*sqrt(alpha));
nc_av=[T1 1];
dc_av=[alpha*T1 1];

sys_c_av=tf(nc_av,dc_av); % FT do compensador
sys_ma_at_av=tf(series(sys_ma_at,sys_c_av)); % FT de MA do sistema 
                                             % compensado por atraso e
                                             % por avanço

figure(20), bode(sys_c_av);                                             
figure(5);
margin(sys_ma_at_av); % Diagrama de Bode da FT de MA do sistema compensado
grid;
set(gcf, 'Position', get(0,'Screensize'));

%% Satisfeita a MF desejada, simulamos o sistema para 
%% 1: resposta ao degrau e rampa

sys_mf_at_av=feedback(sys_ma_at_av,1); % FT de MF do sistema compesado

figure(6),
subplot(131);
margin(sys_ma_at_av);
grid;

subplot(132);
step(sys_mf_at_av);
grid;

% Calcula sobresinal em percentagem e tempo de acomodação
S1 = stepinfo(sys_mf_at_av,'RiseTimeLimits',[0.02,0.98]);
Mp1=S1.Overshoot;
ts1=S1.SettlingTime;

title(['Mp=', num2str(Mp1), '   ts=', num2str(ts1)]);

subplot(133);
t=0:0.1:100;

u=t;
ysc_at=lsim(sys_mf_at,u,t);
ysc_at_av=lsim(sys_mf_at_av,u,t);
% yga=lsim(mf_sc,u,t);
% ycc=lsim(sys_mf_ct,u,t);
plot(t,u,'k',t,ysc,'b',t,ysc_at,'r',t,ysc_at_av,'g-');
legend('Rampa','Sistema','Atraso','At-Av','Location','NorthWest')
grid;

set(gcf, 'Position', get(0,'Screensize'));

%% Carta de Nichols para sistema sem controle, 
% sistema com ganho ajustado e sistema compensado
figure(7)
nichols(tf(np,dp),'r',tf(n1,d1),'b',sys_ma_at_av,'k')
h = findobj(gcf,'type','line');
set(h,'linewidth',2);
set(gcf, 'Units','centimeters', 'Position',[0 1.4 35 25])

legend('Planta','Ganho','Com controle')
ngrid;
