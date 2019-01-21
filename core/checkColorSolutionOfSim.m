%
%check correctness of a graph color solution
%
%checkType: 1 graph coloring, 2 max indp set
%
function [colMatchFound,errPercentage,errorNodes] = checkColorSolutionOfSim(Conn, Maps, Mhistory, verbose, tPoint, checkType )
if nargin<6
    checkType=1;
end
if nargin<4
    verbose=1;
end
if nargin<5
    tPoint= size( Mhistory(1).S, 2 );  % last timepoint
end
if isempty(tPoint)
     tPoint= size( Mhistory(1).S, 2 );  % last timepoint   
end

%check each node, which was the winner
winnerOfNodes=zeros(1,length(Mhistory));
for j=1:length(Mhistory)
    endVals = Mhistory(j).S(1:end-1,tPoint);   %tPoint=end usually (last simulation timepoint)
    maxVal =  max(endVals);
    
    if maxVal==0 & verbose
    %has no winner
        disp(['warning: node has no winner ' num2str(j)]);
    end
    
    indWin = find( endVals == maxVal );
    winnerOfNodes(j)=indWin(1);
end
    
%go through all edges
colMatchFound = 0;
errorNodes = [];
for j=1:length(Conn)
    indFrom = Conn(j).indFrom;
    indTo = Conn(j).indTo;
    
    switch(checkType)
        case 1 %not same (graph coloring)
            if winnerOfNodes(Conn(j).indFrom) == winnerOfNodes(Conn(j).indTo)
                if verbose
                    disp(['warning: col match found; indFrom/To= ' num2str([indFrom indTo]) ' winners=' num2str([ winnerOfNodes(indFrom) winnerOfNodes(indTo)]) ' #E=' num2str(length(Maps(indFrom).indFromList)) ]);
                end
                
                colMatchFound = colMatchFound+1;
                errorNodes = [ errorNodes Conn(j).indFrom Conn(j).indTo ];
            end
        case 2 %max indp set
            %first color not same
            if winnerOfNodes(Conn(j).indFrom) == winnerOfNodes(indTo) && winnerOfNodes(Conn(j).indTo)==1
                colMatchFound = colMatchFound+1;
                errorNodes = [ errorNodes Conn(j).indFrom Conn(j).indTo ];
            end
                        
    end
end


if checkType==2
    %go through all nodes to check that if they have the second color, at least one neighbor has the first (max indp set)
    
    redWinners = find(winnerOfNodes==2);
    
    for j=1:length(redWinners)
    
        redNode = redWinners(j);
        
        %find all neighbors
        neighbors=[];
        for kk=1:length(Conn)
            if Conn(kk).indFrom==redNode || Conn(kk).indTo==redNode
                neighbors=[neighbors Conn(kk).indFrom Conn(kk).indTo];
            end
        end
        
        neighbors=unique(neighbors);
        
        hasBlueNeighbor=0;
        for kk=1:length(neighbors)
           if winnerOfNodes(neighbors(kk))==1
               hasBlueNeighbor=1;
               break;
           end
        end
        
        if ~hasBlueNeighbor
           %this graph has an error
           colMatchFound = colMatchFound+1;
           errorNodes = [ errorNodes redNode  ];           
        end
    end
end

errPercentage = (colMatchFound/2)/length(Conn);
errorNodes = unique(errorNodes);
