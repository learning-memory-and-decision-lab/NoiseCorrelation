% make_hebbNoiseCorrFig.m

% panels:
% 1 -- task screenshots
% 2 -- CP example run
% 3 -- odd example run
% 4 -- CP model LR/vars
% 5 -- odd model lr/vars


num        = 1;
wid        = 17.6; % total width
hts        = [6, 6];
cols       = {1, 2};
[axs,fig_] = getPLOT_axes(num, wid, hts, cols, [3], [3], [], '');
set(axs,'Units','normalized');
% draw in each panel, one at a time

toShow=[1, 11, 21]

getCbColors

% colors: 7 = optimal readout:
%         3 = learned readout:
%         5 = left option
%         6 = right option





for xx = 1:length(axs)
    axes(axs(xx)); hold on; cla(gca)
    if xx==1
        
        set(gca, 'visible', 'off')
        
    elseif xx==2
        
        noiseCorrLayer1=corr((res_FR1));
        imagesc((corr((res_FR1))), [-1, 1])
        title('Layer 1 Noise Correlation')
        colorbar
        xlim([0, 200])
        ylim([0, 200])
        
        set(gca, 'box', 'off')
        
        
    elseif xx==3
         noiseCorrLayer2=corr((res_FR2));
       
        title('Layer 2 Noise Correlation')
        imagesc((corr((res_FR2))), [-1, 1])
        xlim([0, 200])
        ylim([0, 200])
        
        colorbar
        
    end
    
    setPLOT_panelLabel(gca, xx);
end


%% MRN added this section 7-8-20 to get stats for paper:

% Select in-pool and out-pool correlations:
inPoolLogical=tril(true(size(noiseCorrLayer1)), -1);
inPoolLogical(101:200, 1:100)=false;
outPoolLogical=tril(true(size(noiseCorrLayer1)), -1);
outPoolLogical(1:100, 1:100)=false;
outPoolLogical(101:200, 101:200)=false;

% in pool correlations: LAYER 1
layer1inPool=noiseCorrLayer1(inPoolLogical)
[H1_in,P1_in,CI1_in,STATS1_in]=ttest(layer1inPool)
layer1outPool=noiseCorrLayer1(outPoolLogical)
[H1_out,P1_out,CI1_out,STATS1_out]=ttest(layer1outPool)

% in pool correlations:  LAYER 2
layer2inPool=noiseCorrLayer2(inPoolLogical)
[H2_in,P2_in,CI2_in,STATS2_in]=ttest(layer2inPool)
layer2outPool=noiseCorrLayer2(outPoolLogical)
[H2_out,P2_out,CI2_out,STATS2_out]=ttest(layer2outPool)

% Difference in correlation structure (Layer 2 - layer 1)
[H3_in,P3_in,CI3_in,STATS3_in]=ttest2(layer2inPool, layer1inPool)
[H3_out,P3_out,CI3_out,STATS3_out]=ttest2(layer2outPool, layer1outPool)

% mean in pool correlations for layer 1:
mean(layer1inPool)
std(layer1inPool)

% mean in pool correlations for layer 2:
mean(layer2inPool)
std(layer2inPool)

%% 











%kk=annotation('textbox')
fn=fullfile(pwd, 'NC_Fig3.eps');
%set(kk, 'string', 'Nassar et al 2018 Figure 1', 'position', [0.85 0.95 0.15 0.05], 'EdgeColor', 'none')
saveas(gcf,  fn, 'epsc2')


close(gcf)


fn=sprintf('NC_hebbLearnWorkspce_%s',  date )
save(fn)

