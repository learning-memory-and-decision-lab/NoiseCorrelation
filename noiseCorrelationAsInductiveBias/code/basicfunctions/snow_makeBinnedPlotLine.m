function [a B p R]=snow_makeBinnedPlotLine(x,y, PCTs, XLABEL, YLABEL, markSize, markColor, inCurrFig, lineOn)
% bin size is in fractional units... ie binSize=.01 will make a hundred bins, each
% of which includes 1 percent of the date


 


if nargin<9|isempty(lineOn)
    lineOn=0
end



if nargin<6|isempty(markSize)
    markSize=3
end

if nargin<7|isempty(markColor)
    markColor='k'
end

if nargin<8|isempty(inCurrFig)|inCurrFig==0
a=figure;
else
a=gcf;
end

if size(y, 2)>size(y, 1)
    y=y'
end


if size(x, 2)>size(x, 1)
    x=x'
end

[B,BINT,R,RINT,STATS] = regress(y,[ones(length(x), 1) x]);

if lineOn==1
plot(PCTs, [B(1)+PCTs(1).*B(2) B(1)+PCTs(2).*B(2)], 'color', markColor)
end


R = corr(y, x);
p=STATS(3);


hold on
% plot(ErrLinesX, ErrLinesY, 'k', 'lineWidth', 1)
% plot(bErrLinesX, bErrLinesY, 'k', 'lineWidth', 1)
% plot(meanX, meanY, 'o', 'markersize', markSize, 'markerEdgeColor', 'k', 'markerFaceColor', markColor, 'lineWidth', 1)
xlabel(XLABEL)
ylabel(YLABEL)

