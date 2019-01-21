%
%helper function to keep a running min/max of a number of values
%
%urut/sept10/MPI
function [minVal,maxVal] = updateMinMaxVals(minVal,maxVal, newVals)

    minValNew=min(newVals);
    maxValNew=max(newVals);
    
    if minVal>minValNew
        minVal=minValNew;
    end
    if maxVal<maxValNew
        maxVal=maxValNew;
    end
