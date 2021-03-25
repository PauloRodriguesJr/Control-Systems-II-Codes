function []=save_fig_pdf(varargin)
name= varargin{1}; 
fig = varargin{2};
ax  = varargin{3};

% ***** Salvar figura em PDF ****
set(ax, 'Position', get(gca, 'OuterPosition') - ...
get(ax, 'TightInset') * [-3 0 4 2; % original [-1 0 1 0;]
                         0 -2 0 1; % original [0 -2 0 1] [-2 -2 0 1]
                         4 0 2 0; %% Movimentação do quadrado % Original [1 0 1 0] 
                         0 0 0 1]); %% Jogo de cores??[0 0 0 1]
                     % laterais , altura, laterais do grafico,altura do
                     % grafico
set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', [8.8 5.5]);
set(fig, 'PaperPositionMode', 'manual');
set(fig, 'PaperPosition', [0 0 8.8 5.5]);
set(fig, 'renderer', 'painters');
print(gcf, '-dpdf', strcat(name, '.pdf'));


% %% Para o Nichols
% get(ax, 'TightInset') * [-3 0 4 2; % original [-1 0 1 0;]
%                          0 -2 0 1; % original [0 -2 0 1] [-2 -2 0 1]
%                          4 0 2 0; %% Movimentação do quadrado % Original [1 0 1 0] 
%                          0 0 0 1]); %% Jogo de cores??[0 0 0 1]
%                      % laterais , altura, laterais do grafico,altura do
%                      % grafico
