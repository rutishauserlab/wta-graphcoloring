%
%
% summarize a batch of runs of same network, but different initial conditions
%
% calculates av trace and av set nr as funct of time
% also see calcEntropy_ofWTAstate.m
%
%urut/july14
function [avTr, avSetNr ] = calcSetStats_ofWTAstate( setsVisited, trOfSet_sorted )

nrTimesteps = size(setsVisited,2);
nrReps = size(setsVisited,1);

avTr = [];   % as a function of time, the average trace
avSetNr = []; %av set nr as funct of time
for k=1:nrTimesteps

    setsThisTime = setsVisited(:,k);
    
    totTr=0;
    
    totEigs=[0 0];
    for jj=1:length( setsThisTime )
        totTr = totTr + trOfSet_sorted(setsThisTime(jj));
        
    end

    avTr(k) = totTr./nrReps;
    
    avSetNr(k) = mean(setsThisTime);
end

