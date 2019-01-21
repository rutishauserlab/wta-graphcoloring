function [Nnodes,Nedges,nrColors,stateTrans,nodeMatrix1] = buildSudokuGraph(pSize,boxsize)
Nnodes = pSize*pSize;   % 81 nodes
nrColors = pSize;   %

nrBoxes = pSize/boxsize;
% states are zero indexed

%state trans
%stateTrans is From, Symbol, To;   From/To notation is zero-based
%symb=1;


stateTrans=[];

%to verify build procedure
nodeMatrix1 = zeros(pSize,pSize);
nodeMatrix2 = zeros(pSize,pSize);

%rows
for rowNr=1:pSize  
    rowNodes=[1:pSize] + (rowNr-1)*pSize;    
    nodeMatrix1(rowNr,:) = rowNodes;
    stateTrans = addAllPairsToTrans( stateTrans, rowNodes);
end

%columns
for columnNr=1:pSize
    columnNodes = [1:pSize:Nnodes-pSize+1]' + (columnNr-1);
    nodeMatrix2(:,columnNr) = columnNodes;
    stateTrans = addAllPairsToTrans( stateTrans, columnNodes);
end

%boxes
for boxNrRow=1:nrBoxes
   for boxNrColumn=1:nrBoxes
       columnInds = [1:boxsize]+(boxNrColumn-1)*boxsize;
       rowInds = [1:boxsize]+(boxNrRow-1)*boxsize;
       boxNodes = nodeMatrix2(rowInds,columnInds);
       stateTrans = addAllPairsToTrans( stateTrans, boxNodes(:));
   end 
end

%verify each transition is unique
indsKeep=[];
for k=1:size(stateTrans)
   
    fromID= stateTrans(k,1);
    toID= stateTrans(k,3);
    
    inds = find( (stateTrans(:,1)==fromID & stateTrans(:,3)==toID) | (stateTrans(:,3)==fromID & stateTrans(:,1)==toID) );
    
    indsKeep = [indsKeep inds(1)]; %only keep one of each
end

%indsKeep = unique(indsKeep);
%stateTrans = stateTrans(indsKeep,:);


Nedges= size(stateTrans,1);

%====debugging
M = transformGraphToAdjacency( stateTrans, Nedges, Nnodes);

figure(55);
subplot(1,2,1);
imagesc(nodeMatrix1);
subplot(1,2,2);
imagesc(nodeMatrix2);
colorbar;
figure(56);
imagesc(M);

title('adjacency graph');
xlabel('node nr');
ylabel('node nr');


function stateTrans = addAllPairsToTrans( stateTrans, memberNodes)
pairs = getAllPossiblePairs( memberNodes, memberNodes );   %36 pairs in each row/column

% all possible pairs of rowNr and toColumns
for j=1:size(pairs,1)
    stateTrans = addTrans(stateTrans, pairs(j,1), pairs(j,2) );
end

function stateTrans = addTrans(stateTrans, fromS,toS)
stateTrans(size(stateTrans,1)+1,:)= [fromS-1 1 toS-1];   %-1 to make it zero-based

