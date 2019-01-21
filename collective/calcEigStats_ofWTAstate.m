%
%
%av eigenvalue max/min/tr as a function of time across a block of runs of the same network with different inputs 
%
%urut/june2014
function avEigVals = calcEigStats_ofWTAstate(DD_all, N, W, loadTerms, nrSteps)

avEigVals=[]; 


all_maxEig=[];
all_minEig=[];
all_trEig=[];

for runNr=1:length(DD_all)
    DD=DD_all(runNr).DDstats;
    
    %determine the winner based on the set number
    
    winnerID = find( DD(1:end-1,end)== 1);
    
    [maxEig,minEig,trEig,condNr4] = calcContractStats_forRun(nrSteps, DD, N, W, loadTerms, winnerID );
    
    
    all_maxEig(runNr,:) = maxEig;
    all_minEig(runNr,:) = minEig;
    all_trEig(runNr,:) = trEig;
end



avEigVals = [mean(all_maxEig); mean(all_minEig); mean(all_trEig)];
