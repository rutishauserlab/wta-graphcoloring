%
%for graph coloring using the DB cells
%
%adds a new ring of constraints
%
%urut/may12/in-air from CC2012
function [DBcellsUsed, Wdb,stats] = graph_insertNewDBRing(NrOfColors, DBcellsUsed, Wdb, Maps, nodeIDsForRing, beta2_DB, startColorNr)
if nargin<7
    startColorNr=1;
end

stats=[]; %ringNr NodeNr 

%create new DB cells for this new ring
nrCBsNeeded = NrOfColors;
DBindsFrom = DBcellsUsed+1;
DBindsTo   = DBcellsUsed+nrCBsNeeded;

DBinds = DBindsFrom:DBindsTo;
DBcellsUsed = DBindsTo;

% insert the connections
for colInd=startColorNr:NrOfColors   %over these colors
    for jj=1:length(nodeIDsForRing)    %nr of constraints that this ring has; this includes self
        
        targetMapInd = nodeIDsForRing(jj);
        
        Wdb( DBinds(colInd), Maps(targetMapInd+1).overallInd(colInd) ) = beta2_DB;
        
        %if colInd==1  %stats of which node is part of which ring(s)
            stats = [ stats; DBinds(colInd) targetMapInd];
        %end
    end
    
    %also,include itself
    %Wdb( DBinds(colInd), Maps(targetMapInd+1).overallInd(colInd) ) = beta2_DB;
    
    
    
end