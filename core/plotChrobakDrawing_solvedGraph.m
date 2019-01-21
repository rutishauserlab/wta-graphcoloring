%
% plot graph coloring solution (BGL library)
%
function plotChrobakDrawing_solvedGraph( stateTrans, Nedges, Nnodes, nodeWinners, errorNodes,density, colMatchFound)

M = transformGraphToAdjacency( stateTrans, Nedges, Nnodes);

if ~test_planar_graph(M)
    warning('no straight line drawing possible - not planar');
    return;
end

X3 = chrobak_payne_straight_line_drawing(sparse(M));   % straight line drawing

%illustrate the solution
gplot( M, X3, '-o');

hold on
h=gscatter(X3(:,1),X3(:,2), nodeWinners, 'rgbm');
plot( X3(errorNodes,1), X3(errorNodes,2), 'k*','MarkerSize',10);

% mark text
for j=1:size(X3,1)   
    text( X3(j,1)+0.1,X3(j,2),num2str(j),'FontWeight','bold');
end

hold off
set(h,'MarkerSize',20);


ylim([-1 max(X3(:,2))+1]);

title(['Solved graph; ' 'density=' num2str(density) ' nNodes=' num2str(Nnodes) ' nEdges=' num2str(Nedges) ' #errors=' num2str(colMatchFound) ' (stars)']);
