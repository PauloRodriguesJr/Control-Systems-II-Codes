function []=save_fig_pdf(varargin)
name= varargin{1}; 
fig = varargin{2};
ax  = varargin{3};

% ***** Salvar figura em PDF ****
set(ax, 'Position', get(gca, 'OuterPosition') - ...
get(ax, 'TightInset') * [-1 0 1 0;
                         0 -2 0 1; % original -1
                         1 0 1 0;
                         0 0 0 1]);
set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', [8.8 5.5]);
set(fig, 'PaperPositionMode', 'manual');
set(fig, 'PaperPosition', [0 0 8.8 5.5]);
set(fig, 'renderer', 'painters');
print(gcf, '-dpdf', strcat(name, '.pdf'));
