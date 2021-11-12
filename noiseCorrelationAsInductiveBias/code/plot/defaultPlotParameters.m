% Default plot parameters
% width 8.25 inches fontsize of 20

% Robert Wilson
% rcwilson@seas.upenn.edu
% rcw2@princeton.edu

% 10-Dec-2008
    

fontsize = 20;
ABCfontsize = 30;
fontweight = 'bold';
linewidth = 3;

colormap('gray');
C = colormap;
C2 = flipdim(C,1);
colormap(C2)

set(0, 'defaultfigurecolor', 'w', ...
    'defaultfigurecolormap', C2, ...
    'defaultaxesfontsize', fontsize, ...
    'defaultaxesfontweight', fontweight, ...
    'defaultaxestickdir', 'out', ...
    'defaultaxesbox', 'on', ...
    'defaultlinelinewidth', linewidth, ...
    'defaultlinemarkersize', 20, ...
    'defaultfigureposition', [811   486   618   500])


grey = [1 1 1] * 0.5;
lightGrey = [1 1 1]*0.7;
darkGrey = [1 1 1]*0.3;