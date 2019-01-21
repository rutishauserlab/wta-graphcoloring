%
% write XML file specifying a graph, in JFLAP format
%
% G is a sparse matrix representation 
%
%urut/MPI/jan12
function exportGraphToXML( fnameOut, G,nodeWinner)
%fnameOut='/tmp/test.xml';
if nargin<3
    nodeWinner=[];
end

nrStates = size(G,1);

Gu = triu(G,1); % only upper triangular part, not including diagonal is of interest
nrTransitions = nnz(Gu);

% Create a sample XML document.
docNode = com.mathworks.xml.XMLUtils.createDocument('structure');

docRootNode = docNode.getDocumentElement;

thisType = docNode.createElement('type'); 
thisType.appendChild(docNode.createTextNode( 'fa' ));

docRootNode.appendChild(thisType);

automaton = docNode.createElement('automaton'); 

docRootNode.appendChild(docNode.createComment(['Exported Graph from exportGraphToXML.m on ' datestr(now) ]));


for i=1:nrStates
    thisState = docNode.createElement('state'); 
    
    nodeID=i; 
    attr=docNode.createAttribute('id');
    attr.setNodeValue(num2str(nodeID));
    
    attr2=docNode.createAttribute('name');
    attr2.setNodeValue(['q' num2str(nodeID)]);  
    
    thisState.setAttributeNode(attr);
    thisState.setAttributeNode(attr2);

    if length(nodeWinner)>0
        thisLabel = docNode.createElement('label'); 
        thisLabel.appendChild(docNode.createTextNode( num2str(nodeWinner(i)) )   );
        thisState.appendChild( thisLabel );
    end
    
    automaton.appendChild(thisState);
end

[inds1,inds2]=find( Gu );

for i=1:nrTransitions
    thisTransition = docNode.createElement('transition'); 
    
    thisFrom = docNode.createElement('from'); 
    thisFrom.appendChild(docNode.createTextNode( num2str( inds1(i) )   ));
    
    thisTo = docNode.createElement('to'); 
    thisTo.appendChild(docNode.createTextNode( num2str( inds2(i) )   ));
    
    thisRead = docNode.createElement('read'); 
    thisRead.appendChild(docNode.createTextNode( 'A'  ));
    
    thisTransition.appendChild( thisFrom );
    thisTransition.appendChild( thisTo );
    thisTransition.appendChild( thisRead );
    
    automaton.appendChild(thisTransition);
end

docRootNode.appendChild(automaton);

xmlwrite(fnameOut, docNode);