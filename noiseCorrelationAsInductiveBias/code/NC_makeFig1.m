
num        = 1;
wid        = 17.6; % total width
hts        = [7, 7];
cols       = {2,2};
[axs,fig_] = getPLOT_axes(num, wid, hts, cols, [2.5], [2.5], [], '');
set(axs,'Units','normalized');
% draw in each panel, one at a time

ms=7
boundSize=1.5;
for xx = 1:length(axs)
    axes(axs(xx)); hold on; cla(gca)
    if xx==3
        
        
        
        firingRate=allFR(:,:,1);
        hold on
        c=error_ellipse(allCovMat(:,:,1),[targFR, targFR], 'conf', .85)
        d=error_ellipse(allCovMat(:,:,1),[nonTargFR, nonTargFR], 'conf', .85)
        set(c, 'color', 'b')
        set(d, 'color', 'r')
        
        hold on
        plot([-boundSize, boundSize]+.5, [boundSize, -boundSize]+.5, '-k');
%         a=plot(firingRate(1,1:nPoints)', firingRate(2,1:nPoints)', 'ob', 'markerFaceColor','b', 'markerEdgeColor', 'k', 'lineWidth', 1, 'markerSize', ms);
%         b=plot(firingRate(1,(nPoints+1):length(firingRate))', firingRate(2,(nPoints+1):length(firingRate))', 'or', 'markerFaceColor','r', 'markerEdgeColor', 'k', 'lineWidth', 1, 'markerSize', ms);
        ylabel('Neuron 1 firing rate');
        xlabel('Neuron 2 firing rate');
        set(gca, 'box', 'off', 'xticklabel', '', 'yticklabel', '');
        title('Uncorrelated Noise')
        xlim([-2, 3])
        ylim([-2, 3])
        
    elseif xx==4
        firingRate=allFR(:,:,2);
        
        
        hold on
        c=error_ellipse(allCovMat(:,:,2),[targFR, targFR], 'conf', .85)
        d=error_ellipse(allCovMat(:,:,2),[nonTargFR, nonTargFR], 'conf', .85)
        set(c, 'color', 'b')
        set(d, 'color', 'r')
        set(gca, 'box', 'off', 'xticklabel', '', 'yticklabel', '');
        
        
        hold on
        plot([-boundSize, boundSize]+.5, [boundSize, -boundSize]+.5, '-k');
%       a=plot(firingRate(1,1:nPoints)', firingRate(2,1:nPoints)', 'ob', 'markerFaceColor','b', 'markerEdgeColor', 'k', 'lineWidth', 1, 'markerSize', ms);
%       b=plot(firingRate(1,(nPoints+1):length(firingRate))', firingRate(2,(nPoints+1):length(firingRate))', 'or', 'markerFaceColor','r', 'markerEdgeColor', 'k', 'lineWidth', 1, 'markerSize', ms);
     
        ylabel('Neuron 1 firing rate');
        xlabel('Neuron 2 firing rate');
        set(gca, 'box', 'off', 'xticklabel', '', 'yticklabel', '');
        title('Correlated Noise')
        xlim([-2, 3])
        ylim([-2, 3])
        
    elseif xx==1
        
        hold on
        
        a=0;
        %plot([-boundSize, boundSize]+.5, [boundSize, -boundSize]+.5, '-k');
        C=error_ellipse([1, a; a, 1] ,[targFR, targFR], 'conf', .85)
        D=error_ellipse([1, a; a, 1],[nonTargFR, nonTargFR], 'conf', .85)
        set(C, 'color', 'b')
        set(D, 'color', 'r')
        hold on
        a=.5;
        
        E=error_ellipse([1, a; a, 1],[targFR, targFR], 'conf', .85)
        F=error_ellipse([1, a; a, 1],[nonTargFR, nonTargFR], 'conf', .85)
        
        ff=legend([C, D], {'Target', 'Non-target'});
        set(ff, 'location', 'northwest', 'box', 'off')
        
        set(E, 'color', 'b', 'LineStyle', '--')
        set(F, 'color', 'r', 'LineStyle', '--')
        set(gca, 'box', 'off', 'xticklabel', '', 'yticklabel', '');
        xlim([-3, 4])
        ylim([-3, 4])
        ylabel('Neuron 1 firing rate');
        xlabel('Neuron 2 firing rate');
        title('Fixed variance')

        
        
    elseif xx ==2
        
        hold on
        %plot([-boundSize, boundSize]+.5, [boundSize, -boundSize]+.5, '-k');
        c=error_ellipse(allCovMat(1:2,1:2,1),[targFR, targFR], 'conf', .85)
        d=error_ellipse(allCovMat(1:2,1:2,1),[nonTargFR, nonTargFR], 'conf', .85)
        set(c, 'color', 'b')
        set(d, 'color', 'r')
        e=error_ellipse(allCovMat(1:2,1:2,2),[targFR, targFR], 'conf', .85)
        f=error_ellipse(allCovMat(1:2,1:2,2),[nonTargFR, nonTargFR], 'conf', .85)
        set(e, 'color', 'b', 'LineStyle', '--')
        set(f, 'color', 'r', 'LineStyle', '--')
        fff=legend([c, d], {'Target', 'Non-target'});
        set(fff, 'location', 'northwest', 'box','off')
        set(gca, 'box', 'off', 'xticklabel', '', 'yticklabel', '');
        xlim([-2, 3])
        ylim([-2, 3])
        ylabel('Neuron 1 firing rate');
        xlabel('Neuron 2 firing rate');
        title('Fixed Signal-to-Noise')
    end
    
    
    setPLOT_panelLabel(gca, xx);
end

kk=annotation('textbox')
set(kk, 'string', 'Nassar and Bhandari Figure 1', 'position', [0.85 0.95 0.15 0.05], 'EdgeColor', 'none')



saveas(gcf,  'NC_figure1.fig', 'fig')
saveas(gcf,  'NC_figure1.eps', 'epsc2')
close(gcf)