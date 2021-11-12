% make_perceptLearnFig


% panels:
% 1 -- schematic (empty)
% 2 -- Learning curve (rel-pool)    -- add optimal readout?
% 3 -- Learning curve (irrel pool)  -- add optimal readout?
% 4 -- weight diff (or dist to bound?) rel pool 
% 5 -- weight diff (or dist to bound?) irrel pool
% 6 -- rel pool update cloud
% 7 -- same pool update cloud
% 8 -- irrel pool update cloud


% relCorrSim=load(relCorrWorkspace);
% irrelCorrSim=load(irrelCorrWorkspace);

num        = 1;
wid        = 17; % total width
hts        = [5, 5];
cols       = {2,2};
[axs,fig_] = getPLOT_axes(num, wid, hts, cols, [2.0], [2.0], [], '');
set(axs,'Units','normalized');
% draw in each panel, one at a time

getCbColors

LIMS=[.5, 6.5]



for xx = 1:length(axs)
    axes(axs(xx)); hold on; cla(gca)
    if xx==1
        imagesc(nanmean(relCorrSim.meanAccMat, 3),[.5,1])
        colorbar
        title('Learned readout')
        %xlabel('Irrelevent pool correlations')
        ylabel('Same pool correlations')
        set(gca, 'xtick', 2:2:6, 'xticklabel', relCorrSim.label(2:2:6), 'ytick', 2:2:6, 'yticklabel', relCorrSim.label(2:2:6), 'clim', [.5, .8])
        
        xlim(LIMS)
        ylim(LIMS)
        set(gca, 'box', 'off')
 
       
    elseif xx==2
              
        imagesc(nanmean(relCorrSim.optAccMat, 3),[.5,1])
        set(gca, 'xtick', 2:2:6, 'xticklabel', relCorrSim.label(2:2:6), 'ytick', 2:2:6, 'yticklabel', [], 'clim', [.5, .8])
        %xlabel('Irrelevent pool correlations')
       
        title('Optimal readout')
        colorbar
                xlim(LIMS)
        ylim(LIMS)

        set(gca, 'box', 'off')
 
        
        
    elseif xx==3
        
        imagesc(nanmean(irrelCorrSim.meanAccMat, 3),[.5,1])
        colorbar
        %title('Learned readout')
        xlabel('Irrelevent pool correlations')
        ylabel('Same pool correlations')
        set(gca, 'xtick', 2:2:6, 'xticklabel', irrelCorrSim.label(2:2:6), 'ytick', 2:2:6, 'yticklabel', irrelCorrSim.label(2:2:6), 'clim', [.5, .8])
                xlim(LIMS)
        ylim(LIMS)

        set(gca, 'box', 'off')
        
        
    elseif xx==4
        
        
        imagesc(nanmean(irrelCorrSim.optAccMat, 3),[.5,1])
        set(gca, 'xtick', 2:2:6, 'xticklabel', irrelCorrSim.label(2:2:6), 'ytick', 2:2:6, 'yticklabel', [], 'clim', [.5, .8])
        xlabel('Irrelevent pool correlations')
       
        %title('Optimal readout')
        colorbar
                xlim(LIMS)
        ylim(LIMS)

        set(gca, 'box', 'off')
        
        
    end
    
    setPLOT_panelLabel(gca, xx);
end

%kk=annotation('textbox')
fn=fullfile(pwd, 'NC_Fig4S1.eps');
%set(kk, 'string', 'Nassar et al 2018 Figure 1', 'position', [0.85 0.95 0.15 0.05], 'EdgeColor', 'none')
saveas(gcf,  fn, 'epsc2')
close(gcf)

fn=sprintf('NC_abastractTaskWorkspce_%s',  date )
save(fn)


