%
% export solved graph as dot file, which can then be layouted with graphviz.
%
% runMode: 0 only export, 1 also run to produce pdf file of layouted graph
%
%
%urut/may16
function exportGraph_dotFile( stateTrans, Nedges, Nnodes, errorNodes, winnerID, runMode, fnameDot, fnameLayouted )

M = transformGraphToAdjacency( stateTrans, Nedges, Nnodes);

%== export as dot file
nodeLabel=[];
nodeColor=[];
colStrs={'red','green','blue','magenta','yellow','crimson','cyan4','deeppink','gainsboro'};
for k=1:Nnodes
    nodeLabel{k}=num2str(k);
    if ~isempty(find(errorNodes==k))
        nodeLabel{k}=[nodeLabel{k} 'XX'];
    end
    if length(winnerID)>=k
        nodeColor{k}=colStrs{winnerID(k)};
    else
        nodeColor{k}=colStrs{1}; %default color if not set
    end
end

graph_to_dot(M,'directed',0,'filename',fnameDot,'node_label',nodeLabel,'node_color',nodeColor, 'width',60,'height',60, 'layout', 'dot', 'splines','true','nodesep', '0.4', 'overlap','false');

if runMode
    %call dot
    %fdp -Tpdf t.dot >ttt.pdf
    %dot -Tpdf t.dot >ttt.pdf
    cmd=['dot -Tpdf ' fnameDot ' >' fnameLayouted];
    disp(['running command: ' cmd]);
    system(cmd);  % graphviz needs to be installed and working on the system for this to work
end
