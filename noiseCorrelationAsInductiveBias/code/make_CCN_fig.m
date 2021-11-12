% make_CCN_fig

% panels:
% 1 -- task screenshots
% 2 -- CP example run
% 3 -- odd example run
% 4 -- CP model LR/vars
% 5 -- odd model lr/vars


num        = 1;
wid        = 17.6; % total width
hts        = [6, 6];
cols       = {2, 2};
[axs,fig_] = getPLOT_axes(num, wid, hts, cols, [], [], [], '');
set(axs,'Units','normalized');
% draw in each panel, one at a time

getCbColors

for xx = 1:length(axs)
    axes(axs(xx)); hold on; cla(gca)
    if xx==1
        
        set(gca, 'visible', 'off')
        
    elseif xx==2
        
        nCorrs=length(input.all_fracCorrNoise);
        hold on
        clear testAcc testOptAcc
        %a=plot(optWeights, '--k');
        for i = 1:(nCorrs)
            plot(meanAcc(:,i), 'color', ones(3,1).*i./(nCorrs+1));
            testAcc(i)=nanmean(meanAcc(testTrials,i));
            testOptAcc(i)=nanmean(meanOptAcc(testTrials,i));
        end
        ylabel('Accuracy')
        xlabel('Time')
        set(gca, 'box', 'off')
        
        
    elseif xx==3

        hold on
        c=plot(input.all_fracCorrNoise, testOptAcc, 'b')
        d=plot(input. all_fracCorrNoise, testAcc, 'r')
        aa=legend([c, d], 'Optimal readout', 'Learned readout')
        set(aa, 'location' ,'southeast', 'box', 'off')
        ylabel('Accuracy')
        xlabel('Fraction correlated noise')
        set(gca, 'box', 'off')
        
    elseif xx==4
        
        optWeights=[ones(100,1); ones(100,1).*-1]
        
        
        
        normWeights=squeeze(output.finalWeights(:,1,:));
        normWeights=normWeights./nanstd(normWeights);
        hold on
        %a=plot(optWeights, '--k');
        for i = 1:(nCorrs)
            plot(normWeights(:,i), 'color', ones(3,1).*i./(nCorrs+1));
        end
        plot(optWeights, '--b')

        ylabel('Normalized weight')
        xlabel('Neuron')
        set(gca, 'box', 'off')

    end
    
    setPLOT_panelLabel(gca, xx);
end

%kk=annotation('textbox')
fn=fullfile(pwd, 'CCN_fig1.eps');
%set(kk, 'string', 'Nassar et al 2018 Figure 1', 'position', [0.85 0.95 0.15 0.05], 'EdgeColor', 'none')
saveas(gcf,  fn, 'epsc2')
close(gcf)