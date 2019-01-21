%
% create a random planer graph with nrNodes
% and density 0-1. 1 is maximal nr possible edges,which is 3n-6 (euler)
% nrNodes needs to be square, ie the square root needs to be an integer
%
% uses functions of the matlabBGL library,which is based on the C++ Boost library
%
%
%urut/MPI/jan12
function [M,G] = createRandomPlanarGraph(nrNodes, density, ensureMinOnePerNode)
if nargin<3
    ensureMinOnePerNode=0;
end

s = sqrt(nrNodes);

if ~rem(s,1)==0  % check that this numbers is an integer
    error(['nrNodes=' num2str(nrNodes) ' is an invalid value. nrNodes needs to be X^2']);
end

[G,xy] = grid_graph( s,s  ) ;

M = make_maximal_planar(G);

test_planar_graph(M) % it's planar!

% should have 3n-6 edges
maxNrEdges = 3*nrNodes-6

nrToRemove = floor(maxNrEdges*(1-density));

%% remove some randomly to create random planar graph

%ensureMinOnePerNode = 1;

runRemove=true;
nrRemoved=0;
nrTries=0;
maxNrTrials=nrToRemove*10;

removedList=[];
while(runRemove)
    nrTries=nrTries+1;
    
    [i,j,s] = find(M);
    
    % only non-diagonal (non-equal) entries are of interest
    indsNonDiagonal=find(i~=j);
    i=i(indsNonDiagonal);
    j=j(indsNonDiagonal);
    
    r = randperm( length(i) );
    r = r(1);
    
    remove=1;
    
    %remove if not a diagonal entry
    if i(r)==j(r)
        remove=0;
    end

    if ensureMinOnePerNode & remove
        %check if this would make this node left with no connection; if so, can't remove
        
        indsCovered1 = find( i~=j & i==i(r) );
        indsCovered2 = find( i~=j & j==j(r) );        
        
        if length((indsCovered1))==1 || length((indsCovered2))==1  % only this one is left, don't remove it
            %disp('dont remove');
            remove=0; 
        end
    end

    if remove
       M( i(r),j(r))=0;
       M( j(r),i(r))=0;
       nrRemoved=nrRemoved+1;
    end

    %abort once target nr of edges has been removed
    if nrRemoved>= nrToRemove
        break;
    end
    if nrTries>=maxNrTrials
        warning('Didnt remove enough connections, density higher as requested');
        break;
    end
end


disp('');
