%
%for graph coloring using the DB cells
%
%adds a new ring of positive constraints   from one color to an other
%
%urut/may12/in-air from CC2012
function [DBcellsUsed, Wdb,stats, PyrInds] = graph_insertNewDBRing_positive(colFrom, DBcellsUsed, Wdb, Maps, nodeIDsForRing, beta2P_DB )


toAdd_positiveFeedback=[];    % connections from DB back to Py to add later (these are async);   fromDB toPyInd
stats=[]; %ringNr NodeNr 

%create new DB cells for this new ring
nrCBsNeeded = 1;

DBindsFrom = DBcellsUsed+1;
DBindsTo   = DBcellsUsed+nrCBsNeeded;

DBinds = DBindsFrom:DBindsTo;
DBcellsUsed = DBindsTo;

% insert the connections
PyrInds=[];
for jj=1:length(nodeIDsForRing)    %nr of constraints that this ring has; this includes self
    
    targetMapInd = nodeIDsForRing(jj);
    
    PyrInds(jj) = Maps(targetMapInd+1).overallInd(colFrom);
    
    Wdb( DBinds, PyrInds(jj) ) = beta2P_DB;
    
    
    stats = [ stats; DBinds targetMapInd];
    
end
    