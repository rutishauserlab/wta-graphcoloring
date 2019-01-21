function buildSudokuGraph_export(basepathGraph, simNrRand,pSize,boxsize)

[Nnodes,Nedges,nrColors,stateTrans] = buildSudokuGraph(pSize,boxsize);

M = transformGraphToAdjacency( stateTrans, Nedges, Nnodes);

fnameDot=[basepathGraph '/sudoku_s'  num2str(simNrRand) '.dot'];
fnameLayouted=[basepathGraph '/sudoku_s' num2str(simNrRand) '.pdf'];
fnameOut = [basepathGraph '/sudoku_s' num2str(simNrRand) '.jff'];

exportGraphToXML( [fnameOut], sparse(M) );

exportGraph_dotFile( stateTrans, Nedges, Nnodes, [], [], 1, fnameDot, fnameLayouted );
    

