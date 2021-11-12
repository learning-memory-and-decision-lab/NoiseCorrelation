function zVals=get_zTransformedPredictors(xVals, allXes)
    zVals=(xVals-nanmean(allXes))./nanstd(allXes);