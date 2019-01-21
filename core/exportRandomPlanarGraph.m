function colOfNodes = exportRandomPlanarGraph(nrNodes,density, ensureMinOnePerNode, fnameOut, determineRandBiases)
if nargin<5
    determineRandBiases=0;
end

[M,G] = createRandomPlanarGraph(nrNodes, density, ensureMinOnePerNode );

test_planar_graph(M) % it's planar!
nnz(M)
exportGraphToXML( [fnameOut], M);

if determineRandBiases
    % if enabled, compute a solution (4 color) and return it so that a subset of this solution can be chosen as an initial bias
    g = graph;
    [ind_1,ind_2] = find(M);   % M is sparse matrix, get a list of all entries in this matrix
 
    %add all edges of the graph
    for k=1:length(ind_1)    
        add(g,ind_1(k), ind_2(k));
    end

    %color the graph
    g_colors = color(g,'optimal');

    s = size(g_colors);

    nrColsUsed = s(2);   % should be 4, otherwise not useful

    if nrColsUsed ~= 4
        warning('nr colors is not 4, check');
    end
    
    p_colors = parts(g_colors);
    colOfNodes = zeros(nrNodes,1);
    
    for k=1:4
        colOfNodes(p_colors{k}) = k;
    end
else
    colOfNodes=[];
    %=== plots for debugging purposes
    %X3 = chrobak_payne_straight_line_drawing(M);
    %figure(3);
    %gplot( M, X3, '-x');

    %hold on
    %h=gscatter(X3(:,1),X3(:,2), colOfNodes, 'rgbm');
    %hold off
    
    %figure(10);
    %cdraw(g, g_colors)
    %title(['nr colors used:' num2str(nrColsUsed)])
end



%X3 = chrobak_payne_straight_line_drawing(M);
%figure(3);
%gplot( M, X3, '-x');
