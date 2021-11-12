% make_perceptLearnFig




% Big simulation with no normalization:
% NC_abastractTaskWorkspce_27-Mar-2021




% panels:
% 1 -- schematic (empty)
% 2 -- Learning curve (rel-pool)    -- add optimal readout?
% 3 -- Learning curve (irrel pool)  -- add optimal readout?
% 4 -- weight diff (or dist to bound?) rel pool 
% 5 -- weight diff (or dist to bound?) irrel pool
% 6 -- rel pool update cloud
% 7 -- same pool update cloud
% 8 -- irrel pool update cloud

% nReps = 10000

% relCorrSim=load(relCorrWorkspace);
% irrelCorrSim=load(irrelCorrWorkspace);

num        = 1;
wid        = 17.6; % total width
hts        = [4, 4, 4, 5];
cols       = {1, 2, 2, 3};
[axs,fig_] = getPLOT_axes(num, wid, hts, cols, [2.0], [2.0], [], '');
set(axs,'Units','normalized');
% draw in each panel, one at a time

getCbColors


% Loop through time bins and compute accuracy:
bins=0:2:100;

clear meanPerf1 meanPerf2 semPerf meanPerf
for i = 1:(length(bins)-1)
   
     % rel correlations sim:
     meanPerf1(i,1)=nanmean(nanmean(relCorrSim.accMat(relCorrSim.expMat(:,2)==0&relCorrSim.expMat(:,3)==0, (bins(i)+1):bins(i+1) )));
     meanPerf1(i,2)=nanmean(nanmean(relCorrSim.accMat(relCorrSim.expMat(:,2)==.2&relCorrSim.expMat(:,3)==0, (bins(i)+1):bins(i+1) )));
     meanPerf1(i,3)=nanmean(nanmean(relCorrSim.accMat(relCorrSim.expMat(:,2)==.2&relCorrSim.expMat(:,3)==.2, (bins(i)+1):bins(i+1) )));

     % irrel correlations sim:
     meanPerf2(i,1)=nanmean(nanmean(irrelCorrSim.accMat(irrelCorrSim.expMat(:,2)==0&irrelCorrSim.expMat(:,3)==0, (bins(i)+1):bins(i+1) )));
     meanPerf2(i,2)=nanmean(nanmean(irrelCorrSim.accMat(irrelCorrSim.expMat(:,2)==.2&irrelCorrSim.expMat(:,3)==0, (bins(i)+1):bins(i+1) )));
     meanPerf2(i,3)=nanmean(nanmean(irrelCorrSim.accMat(irrelCorrSim.expMat(:,2)==.2&irrelCorrSim.expMat(:,3)==.2, (bins(i)+1):bins(i+1) )));

     % Create temporary arrays of relevant simulation values:
     tmp1=relCorrSim.accMat(relCorrSim.expMat(:,2)==0 &relCorrSim.expMat(:,3)==0, (bins(i)+1):bins(i+1));
     tmp2=relCorrSim.accMat(relCorrSim.expMat(:,2)==.2&relCorrSim.expMat(:,3)==0, (bins(i)+1):bins(i+1));
     tmp3=relCorrSim.accMat(relCorrSim.expMat(:,2)==.2&relCorrSim.expMat(:,3)==.2, (bins(i)+1):bins(i+1));
     
     tmp4=irrelCorrSim.accMat(irrelCorrSim.expMat(:,2)==0&irrelCorrSim.expMat(:,3)==0, (bins(i)+1):bins(i+1));
     tmp5=irrelCorrSim.accMat(irrelCorrSim.expMat(:,2)==.2&irrelCorrSim.expMat(:,3)==0, (bins(i)+1):bins(i+1));
     tmp6=irrelCorrSim.accMat(irrelCorrSim.expMat(:,2)==.2&irrelCorrSim.expMat(:,3)==.2, (bins(i)+1):bins(i+1));
     
     % Compute SEM for each case:
     semPerf(i,1)=nanstd([tmp1(:); tmp4(:)])./sqrt(nReps.*2);     
     semPerf(i,2)=nanstd([tmp2(:); tmp5(:)])./sqrt(nReps.*2);
     semPerf(i,3)=nanstd(tmp3(:))./sqrt(nReps);
     semPerf(i,4)=nanstd(tmp6(:))./sqrt(nReps);

     
end

meanPerf(:, 1:2)=(meanPerf1(:, 1:2)+meanPerf2(:, 1:2))./2;
meanPerf(:,3)=meanPerf1(:,3);
meanPerf(:,4)=meanPerf2(:,3);

% colors: 7 = optimal readout:
%         3 = learned readout:
%         5 = left option
%         6 = right option


for xx = 1:length(axs)
    axes(axs(xx)); hold on; cla(gca)
    if xx==1
        
        set(gca, 'visible', 'off')
        
    elseif xx==2

                % Note -- first two lines are exactly the same in the two
        % simulations -- so we should be able to average across them. 

        defaultPlotParameters
        hold on
        plot([bins(2), bins(end)], [.5, .5], '--k')
          
        H1=shadedErrorBar(bins(2:end),meanPerf(:,1),semPerf(:,1),{'color', cbColors(2,:)})
        H2=shadedErrorBar(bins(2:end),meanPerf(:,2),semPerf(:,2),{'color', cbColors(3,:)})
        H3=shadedErrorBar(bins(2:end),meanPerf(:,3),semPerf(:,3),{'color', cbColors(4,:)})
        ff=legend([H1.mainLine, H2.mainLine, H3.mainLine], relCorrSim.legendText)
        set(ff, 'location', 'east',  'box', 'off')
        ylabel('Percent Correct')
        xlabel('Trials')
        xlim([1,100])
        ylim([.45, .8])
        set(gca, 'box', 'off')

        
        
    elseif xx==3
        
        % Note -- first two lines are exactly the same in the two
        % simulations -- so we should be able to average across them. 

        defaultPlotParameters
        hold on
        plot([bins(2), bins(end)], [.5, .5], '--k')
          
        H1=shadedErrorBar(bins(2:end),meanPerf(:,1),semPerf(:,1),{'color', cbColors(2,:)})
        H2=shadedErrorBar(bins(2:end),meanPerf(:,2),semPerf(:,2),{'color', cbColors(3,:)})
        H3=shadedErrorBar(bins(2:end),meanPerf(:,4),semPerf(:,4),{'color', cbColors(4,:)})
        ff=legend([H1.mainLine, H2.mainLine, H3.mainLine], irrelCorrSim.legendText)
        set(ff, 'location', 'east',  'box', 'off')
        %ylabel('Percent Correct')
        xlabel('Trials')
        xlim([1,100])
        ylim([.45, .8])
        set(gca, 'box', 'off')
        
        
    elseif xx==4
        
        imRange=[-.02, .22]
        
        imagesc(allCorrs, allCorrs, nanmean(relCorrSim.meanDistMat, 3),[.5,1.2])
        set(gca, 'xtick', 1:6, 'xticklabel', label, 'ytick', 1:6, 'yticklabel', label)
        title('Distance from optimal readout')
        xlabel('Relevant pool correlations')
        ylabel('In pool correlations')
        xlim(imRange);
        ylim(imRange);
        set(gca, 'box', 'off')
        colorbar
        
        
        
    elseif xx==5
        
        
        imagesc(allCorrs, allCorrs,nanmean(irrelCorrSim.meanDistMat, 3),[.5,1.2])
        set(gca, 'xtick', 1:6, 'xticklabel', label, 'ytick', 1:6, 'yticklabel', label)
        title('Distance from optimal readout')
        xlabel('irrelevant pool correlations')
        ylabel('In pool correlations')
        xlim(imRange);
        ylim(imRange);
        set(gca, 'box', 'off')
        colorbar
        
        
        
    elseif xx==6
        
        
        k=21
        updateInDims=[ones(1,200), ones(1,200).*-1; ones(1,100), ones(1,100).*-1, ones(1,100), ones(1,100).*-1]*relCorrSim.allStoredWeightUpdates(:,:,k)
        hold on
        Scale=max(abs(updateInDims(:)));
       % title(sprintf('in pool= %g, rel pool= %g, irrel pool = %g', relCorrSim.allCorrCombos.inPool(k), relCorrSim.allCorrCombos.inRel(k), relCorrSim.allCorrCombos.inIrrel(k)))
        plot([-Scale, Scale], [0, 0], '--k')
        plot([0, 0], [-Scale, Scale],  '--k')
        % a=plot(updateInDims(1,allTrialTypes==1), updateInDims(2,allTrialTypes==1), 'or', 'markerFaceColor','r', 'markerEdgeColor', 'k', 'lineWidth', 1, 'markerSize', 8)
        % b=plot(updateInDims(1,allTrialTypes==2), updateInDims(2,allTrialTypes==2), 'ob', 'markerFaceColor','b', 'markerEdgeColor', 'k', 'lineWidth', 1, 'markerSize', 8)
        a=quiver(zeros(sum(relCorrSim.allTrialTypes==1),1), zeros(sum(relCorrSim.allTrialTypes==1),1), updateInDims(1,relCorrSim.allTrialTypes==1)', updateInDims(2,relCorrSim.allTrialTypes==1)', 0,  'r', 'lineWidth', 1)
        b=quiver(zeros(sum(relCorrSim.allTrialTypes==2),1), zeros(sum(relCorrSim.allTrialTypes==2),1), updateInDims(1,relCorrSim.allTrialTypes==2)', updateInDims(2,relCorrSim.allTrialTypes==2)', 0,  'b', 'lineWidth', 1)
        
        
        ff=legend([a, b], 'X or Y', '+ or -')
        set(ff, 'box', 'off', 'location', 'northwest')

        
        
%         ff=legend([a, b], 'X or Y', '+ or -')
%         set(ff, 'box', 'off')
        xlabel('Update in XY')
        ylabel('Update in +-')
        xlim([-Scale, Scale])
        ylim([-Scale, Scale])
         set(gca, 'box', 'off', 'xtick', [], 'xticklabels', {}, 'ytick', [], 'yticklabels', {})
       %fn=sprintf('WeightUpdateQuiver_%g_rel_%g_irrel_%g.eps', allCorrCombos.inPool(k),allCorrCombos.inRel(k), allCorrCombos.inIrrel(k));
         
        
        % absolute update in left/right for left/right trials:
        lrLR=abs(updateInDims(1,relCorrSim.allTrialTypes==1)')
        % absolute update in up/down for left/right trials:
        udLR=abs(updateInDims(2,relCorrSim.allTrialTypes==1)')
%         mean(lrLR)
%         mean(udLR)
         % absolute update in left/right for up/down trials:
        lrUD=abs(updateInDims(1,relCorrSim.allTrialTypes==2)')
        % absolute update in up/down for up/down trials:
        udUD=abs(updateInDims(2,relCorrSim.allTrialTypes==2)')
       
%         mean(lrUD)
%         mean(udUD)

        [H,P,CI,STATS] =ttest2(lrLR-udLR, lrUD-udUD) 
       
        
        
    elseif xx==7

           
        k=6
       updateInDims=[ones(1,200), ones(1,200).*-1; ones(1,100), ones(1,100).*-1, ones(1,100), ones(1,100).*-1]*irrelCorrSim.allStoredWeightUpdates(:,:,k)
    
        hold on
        Scale=max(abs(updateInDims(:)));
       % title(sprintf('in pool= %g, rel pool= %g, irrel pool = %g', irrelCorrSim.allCorrCombos.inPool(k), irrelCorrSim.allCorrCombos.inRel(k), irrelCorrSim.allCorrCombos.inIrrel(k)))
        plot([-Scale, Scale], [0, 0], '--k')
        plot([0, 0], [-Scale, Scale],  '--k')
        % a=plot(updateInDims(1,allTrialTypes==1), updateInDims(2,allTrialTypes==1), 'or', 'markerFaceColor','r', 'markerEdgeColor', 'k', 'lineWidth', 1, 'markerSize', 8)
        % b=plot(updateInDims(1,allTrialTypes==2), updateInDims(2,allTrialTypes==2), 'ob', 'markerFaceColor','b', 'markerEdgeColor', 'k', 'lineWidth', 1, 'markerSize', 8)
        a=quiver(zeros(sum(irrelCorrSim.allTrialTypes==1),1), zeros(sum(irrelCorrSim.allTrialTypes==1),1), updateInDims(1,irrelCorrSim.allTrialTypes==1)', updateInDims(2,irrelCorrSim.allTrialTypes==1)', 0,  'r', 'lineWidth', 1)
        b=quiver(zeros(sum(irrelCorrSim.allTrialTypes==2),1), zeros(sum(irrelCorrSim.allTrialTypes==2),1), updateInDims(1,irrelCorrSim.allTrialTypes==2)', updateInDims(2,irrelCorrSim.allTrialTypes==2)', 0,  'b', 'lineWidth', 1)
        
%         ff=legend([a, b], 'X or Y', '+ or -')
%         set(ff, 'box', 'off')
        xlabel('Update in XY')
        ylabel('Update in +-')
        xlim([-Scale, Scale])
        ylim([-Scale, Scale])
        set(gca, 'box', 'off', 'xtick', [], 'xticklabels', {}, 'ytick', [], 'yticklabels', {})
         %fn=sprintf('WeightUpdateQuiver_%g_rel_%g_irrel_%g.eps', allCorrCombos.inPool(k),allCorrCombos.inRel(k), allCorrCombos.inIrrel(k));

        
         

        
        
        
    elseif xx==8

        
        
         k=21
        updateInDims=[ones(1,200), ones(1,200).*-1; ones(1,100), ones(1,100).*-1, ones(1,100), ones(1,100).*-1]*irrelCorrSim.allStoredWeightUpdates(:,:,k)
     
        hold on
        Scale=max(abs(updateInDims(:)));
        %title(sprintf('in pool= %g, rel pool= %g, irrel pool = %g', irrelCorrSim.allCorrCombos.inPool(k), irrelCorrSim.allCorrCombos.inRel(k), irrelCorrSim.allCorrCombos.inIrrel(k)))
        plot([-Scale, Scale], [0, 0], '--k')
        plot([0, 0], [-Scale, Scale],  '--k')
        % a=plot(updateInDims(1,allTrialTypes==1), updateInDims(2,allTrialTypes==1), 'or', 'markerFaceColor','r', 'markerEdgeColor', 'k', 'lineWidth', 1, 'markerSize', 8)
        % b=plot(updateInDims(1,allTrialTypes==2), updateInDims(2,allTrialTypes==2), 'ob', 'markerFaceColor','b', 'markerEdgeColor', 'k', 'lineWidth', 1, 'markerSize', 8)
        a=quiver(zeros(sum(irrelCorrSim.allTrialTypes==1),1), zeros(sum(irrelCorrSim.allTrialTypes==1),1), updateInDims(1,irrelCorrSim.allTrialTypes==1)', updateInDims(2,irrelCorrSim.allTrialTypes==1)', 0,  'r', 'lineWidth', 1)
        b=quiver(zeros(sum(irrelCorrSim.allTrialTypes==2),1), zeros(sum(irrelCorrSim.allTrialTypes==2),1), updateInDims(1,irrelCorrSim.allTrialTypes==2)', updateInDims(2,irrelCorrSim.allTrialTypes==2)', 0,  'b', 'lineWidth', 1)
        
        
        
                % absolute update in left/right for left/right trials:
        lrLR=abs(updateInDims(1,irrelCorrSim.allTrialTypes==1)')
        % absolute update in up/down for left/right trials:
        udLR=abs(updateInDims(2,irrelCorrSim.allTrialTypes==1)')
%         mean(lrLR)
%         mean(udLR)
         % absolute update in left/right for up/down trials:
        lrUD=abs(updateInDims(1,irrelCorrSim.allTrialTypes==2)')
        % absolute update in up/down for up/down trials:
        udUD=abs(updateInDims(2,irrelCorrSim.allTrialTypes==2)')
       
%         mean(lrUD)
%         mean(udUD)

        [H,P,CI,STATS] =ttest2(lrLR-udLR, lrUD-udUD) 

        
        
        
        xlabel('Update in XY')
        ylabel('Update in +-')
        xlim([-Scale, Scale])
        ylim([-Scale, Scale])
        set(gca, 'box', 'off', 'xtick', [], 'xticklabels', {}, 'ytick', [], 'yticklabels', {})
        %fn=sprintf('WeightUpdateQuiver_%g_rel_%g_irrel_%g.eps', allCorrCombos.inPool(k),allCorrCombos.inRel(k), allCorrCombos.inIrrel(k));

        
    end
    
    setPLOT_panelLabel(gca, xx);
end

%kk=annotation('textbox')
fn=fullfile(pwd, 'NC_Fig4_noNorm_justBias.eps');
%set(kk, 'string', 'Nassar et al 2018 Figure 1', 'position', [0.85 0.95 0.15 0.05], 'EdgeColor', 'none')
saveas(gcf,  fn, 'epsc2')
close(gcf)

close all
fn=sprintf('NC_abastractTaskWorkspce_%s',  date )
save(fn)
