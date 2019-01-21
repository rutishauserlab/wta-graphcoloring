%
% find transitions that lead to either from or to the desired node
% nodeNr is zero based
% onlyCoveredFlag: -1 all, 0 only non-used, 1 only used
%
%
%urut/may12/in-air from capo caccia
function [inds,targetMapIDs] = graph_findTransition_forNode( stateTrans, nodeNr, onlyCoveredFlag )
inds1 = find( stateTrans(:,1) == nodeNr | stateTrans(:,3) == nodeNr);

if onlyCoveredFlag==-1
    inds=inds1;
else
    inds2 = find( stateTrans(:,4)==onlyCoveredFlag );
    
    inds = intersect(inds1,inds2);
end
targetMapIDs=[];
%
% find target ID for each, excluding self

stateEndpoints = stateTrans( inds, [1 3]);
for k=1:length(inds)
    targetMapIDs(k) = setdiff( stateEndpoints(k,:), nodeNr); %setdiff to exclude itself
end
                