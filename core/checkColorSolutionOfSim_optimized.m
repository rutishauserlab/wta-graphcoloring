%
%check correctness of a graph color solution
%
%checkType: 1 graph coloring, 2 max indp set
%
%
%==== this optimized version only supports checkType==1
%
%
function [colMatchFound] = checkColorSolutionOfSim_optimized(Conn_collapsed, Soverall, Maps, nrUnits_perMap )

%== check each node, which was the winner
winnerOfNodes=zeros(1,length(Maps));
for k=1:length(Maps)
    [~,maxInd] = max( Soverall( [1:nrUnits_perMap]+(k-1)*nrUnits_perMap));
    winnerOfNodes(k) = maxInd(1);
end

%== go through all edges
%Conn_collapsed: indFrom indTo
colMatchFound = 0;
errorNodes = [];

colMatchFound = sum(winnerOfNodes( Conn_collapsed(:,1)) == winnerOfNodes(Conn_collapsed(:,2)));


%for j=1:size(Conn_collapsed,1)
%    if winnerOfNodes(Conn_collapsed(j,1)) == winnerOfNodes(Conn_collapsed(j,2))
%        colMatchFound = colMatchFound+1;
%        errorNodes = [ errorNodes Conn_collapsed(j,1) Conn_collapsed(j,2) ];
%    end
%end

%errPercentage = (colMatchFound/2)/size(Conn_collapsed,1);
%errorNodes = unique(errorNodes);
