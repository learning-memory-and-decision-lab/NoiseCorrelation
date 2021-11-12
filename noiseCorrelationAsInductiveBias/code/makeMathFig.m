% makeMathFig

% No normalization:     'NC_perceptLearnWorkspce_27-Mar-2021'

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

getCbColors

for xx = 1:length(axs)
    axes(axs(xx)); hold on; cla(gca)
    
    if xx==1
        set(gca, 'visible', 'off')
    
    elseif xx==2
        
        hold on
        plot(output.phi_list, sqrt(output.tot_var_a)); 
        plot(output.phi_list(1:2:101), sqrt(output.tot_var_e(1:2:101)), 'o')
        title('Total Perpendicular Std.')
        %grid on
        xlabel('\phi')
        ylabel('\sigma(\xi_\perp)')
        ff=legend({'Analytic','Empirical'})
        set(ff, 'box', 'off')
        set(gca, 'box', 'off')
       
        
    elseif xx==3
        
        Skip=2;
        showPoints=1:Skip:length(allWsn(exNC,:,:))
        hold on
        exNC=2 % choose noise correlations
        % Analytical weight components
        a=plot(input.LR*0.5*(1:input.nTrials)*norm(output1.Mu), '-', 'Color', 'b');
        b=plot(input.LR*0.5*sqrt(1:input.nTrials)*sqrt(trace(output1.allCovMat(:,:,exNC))-100), '-', 'Color', 'r');
        % Observed weight components
        hold on
        c=plot(showPoints, mean(allWsn(exNC,showPoints,:), 3), 'o', 'Color', 'b');
        d=plot(showPoints, mean(allWpn(exNC,showPoints,:),3), 'o', 'Color', 'r');
        title('Weight Components')
        xlabel('Trial')
        ylabel('Component Norm')
        ff=legend([ a, b, c,d], {'Analytic \Delta w_s','Analytic \Delta w_{\perp}', 'Empirical \Delta w_s','Empirical \Delta w_{\perp}'}, 'Location', 'NorthWest')
        set(ff, 'box', 'off')
        set(gca, 'box', 'off')
        
    elseif xx==4
        
        for i = output.phi_inds
            plot(output.acc(:,i));   hold on
        end
        %grid on
        title('Analytic Accuracy')
        xlabel('Trial')
        ylabel('Percent Correct')
        ff=legend(output.namelist, 'Location', 'SouthEast')
        set(ff, 'box', 'off')

        set(gca, 'box', 'off')

    elseif xx==5
        
        for i = output.phi_inds
            plot(log(1:length(output3.acc(:,i))), output3.acc(:,i));   hold on
        end
        %grid on
        title('Analytic Accuracy')
        xlabel('log(Trial)')
        ylabel('Percent Correct')
        ff=legend(output.namelist, 'Location', 'SouthEast')
        set(ff, 'box', 'off')

        set(gca, 'box', 'off')

        
        
    elseif xx==6
        
        % pool size is N/2 -- since half neurons are in each pool:
        imagesc(log2(allN./2),allMu.*2,  noiseCorrAdvantage./100, [0, .4]) % now normalizing to average accuracy, rather than sum
        % mu is SNR/2 -- so multiply to get SNR. 
        
        
        hold on
        plot( log(100), 2.*sqrt(2), 'xr', 'markerSize', 12, 'lineWidth', 3)
        
        xlim([.5, 13.5])
        ylim([.5, 10.5])
        colorbar
        title('Accuracy Improvement')
        xlabel('log(Pool size)')
        ylabel('Signal-to-noise ratio')
        set(gca, 'box', 'off')
        
    end
    
    setPLOT_panelLabel(gca, xx);
end

%kk=annotation('textbox')
fn=fullfile(pwd, 'NC_analyticFig.eps');
%set(kk, 'string', 'Nassar et al 2018 Figure 1', 'position', [0.85 0.95 0.15 0.05], 'EdgeColor', 'none')
saveas(gcf,  fn, 'epsc2')
close(gcf)

fn=sprintf('NC_AnalyticWorkspce_%s',  date )
save(fn)

