% make_perceptLearnFig_fixedNoise

% No normalization:
% Full simulation:  'NC_perceptLearnWorkspce_fixedNoise27-Mar-2021'


% panels:
% 1 -- task screenshots
% 2 -- CP example run
% 3 -- odd example run
% 4 -- CP model LR/vars
% 5 -- odd model lr/vars


num        = 1;
wid        = 17.6; % total width
hts        = [6, 6, 6];
cols       = {2, 2, 2};
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
        
        nCorrs=length(input.all_fracCorrNoise);
        
        
        
        hold on
        plot([1, 100], [nanmean(testOptAcc(:)), nanmean(testOptAcc(:))], '-', 'color', cbColors(7,:))
        
        clear testAcc testOptAcc
        %a=plot(optWeights, '--k');
        for i = 1:(nCorrs)
            plot(smMeanAcc(:,i), 'color', ones(3,1).*i./(nCorrs+3));
            testAcc(i)=nanmean(smMeanAcc(testTrials,i));
            testOptAcc(i)=nanmean(meanOptAcc(testTrials,i));
        end
        
        
        
        ylabel('Accuracy')
        xlabel('Time')
        set(gca, 'box', 'off')
        
        
    elseif xx==5

        forPerf    = squeeze(nanmean(allAcc(testTrials,:,:), 1));
        meanPerf   = mean(forPerf');
        semPerf    = nanstd(forPerf')./sqrt(size(forPerf, 2));   

        
        forOptPerf    = squeeze(nanmean(allOptAcc(testTrials,:,:), 1));
        meanOptPerf   = mean(forOptPerf');
        semOptPerf    = nanstd(forOptPerf')./sqrt(size(forOptPerf, 2));   
        
         hold on
       
        H1=shadedErrorBar(input.all_fracCorrNoise, meanOptPerf, semOptPerf, {'color', cbColors(7,:)});
        H2=shadedErrorBar(input.all_fracCorrNoise, meanPerf, semPerf, {'color', cbColors(2,:)});
        
        aa=legend([H1.mainLine, H2.mainLine], 'Optimal readout', 'Learned readout')
        set(aa, 'location' ,'southeast', 'box', 'off')
        ylabel('Final Accuracy')
        xlabel('Noise correlations')
        set(gca, 'box', 'off')
        
    elseif xx==3

        
        input.all_fracCorrNoise(toShow)      
        hold on
        
        hist([normWeightDiff(1:100,toShow(1)), normWeightDiff(101:200,toShow(1))], -3:.3:3)
        b=get(gca, 'Children')
        set(b(1), 'faceColor', cbColors(5,:))
        set(b(2), 'faceColor', cbColors(6,:))

        
        
        
        xlim([-3, 3])
        set(gca, 'box', 'off')
        ylabel('Frequency')
        xlabel('Weights')
        title(sprintf('Noise correlations = %g', input.all_fracCorrNoise(toShow(1))))
        
        
    elseif xx==4
        
        input.all_fracCorrNoise(toShow)      
        hold on
        hist([normWeightDiff(1:100,toShow(3)), normWeightDiff(101:200,toShow(3))], -3:.3:3)
        xlim([-3, 3])
        set(gca, 'box', 'off')
        ylabel('Frequency')
        xlabel('Weights')
        title(sprintf('Noise correlations = %g', input.all_fracCorrNoise(toShow(3))))
        a=get(gca, 'Children')
        set(a(1), 'faceColor', cbColors(5,:))
        set(a(2), 'faceColor', cbColors(6,:))
     
      
      
      
    elseif xx ==6
        
        STD=nanstd(eucDistToBound')
        hold on
        plot([input.all_fracCorrNoise(1), input.all_fracCorrNoise(end)], [optEucDistToBound, optEucDistToBound], 'color', cbColors(7,:), 'lineWidth', 3)
        H=shadedErrorBar(input.all_fracCorrNoise, nanmean(eucDistToBound'), STD, {'color', cbColors(2,:)});
        set(gca, 'box', 'off')
        ylabel('Distance to mis-classification boundary')
        xlabel('Noise correlations')
        
    end
    
    setPLOT_panelLabel(gca, xx);
end

%kk=annotation('textbox')
fn=fullfile(pwd, 'NC_Fig2_fixedNoise.eps');
%set(kk, 'string', 'Nassar et al 2018 Figure 1', 'position', [0.85 0.95 0.15 0.05], 'EdgeColor', 'none')
saveas(gcf,  fn, 'epsc2')
close(gcf)

close all
fn=sprintf('NC_perceptLearnWorkspce_fixedNoise%s',  date )
save(fn)

