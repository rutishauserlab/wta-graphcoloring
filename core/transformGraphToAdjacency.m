%
%transform a graph described by transitions to a sparse adjencency matrix
%
%urut/may12
function Msparse=transformGraphToAdjacency( stateTrans, Nedges, Nnodes)

M = zeros(Nnodes,Nnodes);

for k=1:Nedges    
    M( stateTrans(k,1)+1, stateTrans(k,3)+1 ) = 1;   %make a transition here
    M(  stateTrans(k,3)+1,stateTrans(k,1)+1 ) = 1;   %make it symmetric
end

Msparse = sparse(M);