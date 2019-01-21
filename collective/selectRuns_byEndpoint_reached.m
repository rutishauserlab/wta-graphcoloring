%
%select only runs where the network reached a specific set nr
%this can either be the final state or at any point or at a particular time
%
%
%setNrAll: sets visited in all runs (each run is one row, as a function of time).
%selForSetNr: the set # to select for
%endPointToUse: relative to end, when to check for state; any time if =-1
%
%toUse: index into rows of setNrAll
%
%
%urut/june14
function toUse = selectRuns_byEndpoint_reached( setNrAll, endPointToUse, selForSetNr )

toUse=[];
 
for k=1:size(setNrAll,1)        
    
    if endPointToUse==-1
        
        if ~isempty( find(setNrAll(k,:)==selForSetNr))
            toUse = [toUse k];
        end
    else
        if ~isempty( find( setNrAll(k,end-endPointToUse:end) == selForSetNr ) )
            toUse = [toUse k];
        end
        
    end
    
    
end
