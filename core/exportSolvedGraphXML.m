%
%export a graph as a XML file, with the solution marked as text (label)
%
function M = exportSolvedGraphXML(fnameOut, Nnodes,Nedges,stateTrans, nodeWinners)

M = transformGraphToAdjacency( stateTrans, Nedges, Nnodes);

exportGraphToXML( [fnameOut], sparse(M), nodeWinners );

disp(['wrote ' fnameOut]);