%
%tells whether a solution of a sudoku is correct
%
%
function e=checkSudokuCorrectness(nodeMatrixSolution)

pSize = (size(nodeMatrixSolution,1));

if pSize==9
    blocksize=3;
end
if pSize==4
    blocksize=2;
end

s=sum([1:pSize]);

e=0;

%rows/columns
for k=1:pSize
    vals1 = nodeMatrixSolution(k,:);
    vals2 = nodeMatrixSolution(:,k);
    
    if ~isSumAndUnique(s,vals1,pSize) | ~isSumAndUnique(s,vals2,pSize)
        e=1;
    end
end

%boxes

%boxes
for boxNrRow=1:pSize/blocksize
    for boxNrColumn=1:pSize/blocksize
        columnInds = [1:blocksize]+(boxNrColumn-1)*blocksize;
        rowInds = [1:blocksize]+(boxNrRow-1)*blocksize;
        vals = nodeMatrixSolution(rowInds,columnInds);
        
        if ~isSumAndUnique(s,vals1,pSize) | ~isSumAndUnique(s,vals2,pSize)
            e=1;
        end
    end
end

function corr=isSumAndUnique(s,vals,pSize)
if sum(vals)==s & isempty(setdiff(unique(vals),[1:pSize]))
    corr=1;
else
    corr=0;
end