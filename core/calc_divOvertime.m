%
% calculate the "divergence signal" as a function of time for each WTA in a network.
%
%urut/march16
function divAll = calc_divOvertime( stateHistoryAll, nrNodes, unitsPerMap, threshActive, simTill, alpha1 )

G=-1; %load term

divAll=nan(nrNodes,simTill);   % for each WTA one value at every timepoint
for WTAnr=1:nrNodes
   indsOfUnit = [1:unitsPerMap] + (WTAnr-1)*unitsPerMap;

   for t=1:simTill
       vals = stateHistoryAll( indsOfUnit(1:end-1), t);
       valsThresh = G*ones(1,unitsPerMap-1);   % default is "off" (-1)
       valsThresh(find(vals>threshActive))=alpha1+G;   %all which are one are alpha-1
       
       divAll(WTAnr,t) = -1 + sum(valsThresh);
   end
end

%eliminate initial transients
divAll(:,1:10)=nan;
