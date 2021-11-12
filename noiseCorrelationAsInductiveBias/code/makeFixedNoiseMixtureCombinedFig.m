% make 2 panel figure with
% panel 1 -- learning curves for noise correlations with fixed noise
% panel 2 -- learning curves for mixture models 

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
hts        = [6];
cols       = {2};
[axs,fig_] = getPLOT_axes(num, wid, hts, cols, [2], [2], [], '');
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
        
        
        fnData=load('NC_perceptLearnWorkspce_fixedNoise12-Apr-2021.mat') 
        nCorrs=length(fnData.input.all_fracCorrNoise);
        hold on
        plot([1, 100], [nanmean(fnData.testOptAcc(:)), nanmean(fnData.testOptAcc(:))], '-', 'color', fnData.cbColors(7,:));
        
        clear testAcc testOptAcc
        %a=plot(optWeights, '--k');
        for i = 1:(nCorrs)
            plot(fnData.meanAcc(:,i), 'color', ones(3,1).*i./(nCorrs+3));
            testAcc(i)=nanmean(fnData.meanAcc(fnData.testTrials,i));
            testOptAcc(i)=nanmean(fnData.meanOptAcc(fnData.testTrials,i));
        end
  ylim([0.4, 1])
        ylabel('Accuracy')
        xlabel('Time')
        set(gca, 'box', 'off')
        
        
        
    elseif xx==2
       disp('got here')
        
        
        mData=load('/Users/mattnassar/Dropbox/noiseCorrelationAsInductiveBias/diffAssumptionsWorkspace_3-28-21.mat')
        
        
        hold on
        
        plot([1, 100], [nanmean(fnData.testOptAcc(:)), nanmean(fnData.testOptAcc(:))], '-', 'color', fnData.cbColors(7,:));
      
        a=plot(mData.zeroNC_acc, 'r', 'lineWidth', 2)
        for j = 1:size(mData.someNC_acc, 2)
            frac=j./11;
            z(j)=plot(mData.someNC_acc(:,j), '-', 'color', [.9, .9, .9].*frac, 'lineWidth', 2)
            
        end

        ff=legend([a, z(1), z(end)], 'No NC', '0.2 NC -- fixed unit signal & variance',  '0.2 NC -- fixed SNR')
        
        ylim([0.4, 1])
        set(ff, 'box', 'off', 'location', 'south')
        ylabel('Accuracy')
        xlabel('Trials')
        set(gca, 'box', 'off')
        
        
    end
    
    setPLOT_panelLabel(gca, xx);
end

%kk=annotation('textbox')
fn=fullfile(pwd, 'NC_Fig2_fixedNoiseAndMixture.eps');
%set(kk, 'string', 'Nassar et al 2018 Figure 1', 'position', [0.85 0.95 0.15 0.05], 'EdgeColor', 'none')
saveas(gcf,  fn, 'epsc2')
%close(gcf)

