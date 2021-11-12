function [zscoredVar, transX]=pairedZTransorm(xes, ex_x)
        zscoredVar= (xes-nanmean(xes))./nanstd(xes);
        transX      = (ex_x-nanmean(xes))./nanstd(xes);