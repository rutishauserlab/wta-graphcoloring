%
%
% compute stats on set switch probabilities
%
%urut/july14
function [transitionPairs,dAllTr] = prepareSetSwitchStats( setsVisited, trOfSet_sorted )

nrRuns = size(setsVisited,1);
timesteps = size(setsVisited,2);

transitionPairs = [0 0 0];  % from to count

for runNr=1:nrRuns
    d = diff( setsVisited(runNr,:) );
    
    indJumps = find(d~=0);
    for k=1:length(indJumps)
        fromToPair = setsVisited( runNr, [indJumps(k) indJumps(k)+1] );

        indOfPair = find(transitionPairs(:,1)==fromToPair(1) & transitionPairs(:,2)==fromToPair(2) );
        
        if isempty(indOfPair)
            %add this transition
            transitionPairs = [transitionPairs; [fromToPair 1]];
        else
            %increase count on this transition
            transitionPairs(indOfPair,3) = transitionPairs(indOfPair,3)+1;
        end
    end
end


% probability of transition, conditional on being in a certain state
transitionPairs(:,4)=0;

startPos = unique(transitionPairs(:,1));

for j=1:length(startPos)
   
    
    inds = find(transitionPairs(:,1)==startPos(j));
    
    totNr = sum(transitionPairs(inds,3));
    
    transitionPairs(inds,4) = transitionPairs(inds,3)./totNr;
    
end


% 
dAllTr=[];
for runNr=1:nrRuns
    trsOfRun = trOfSet_sorted( setsVisited(runNr,:) );
    d = diff(trsOfRun);
    dAllTr = [dAllTr d(find(d~=0))];
end
