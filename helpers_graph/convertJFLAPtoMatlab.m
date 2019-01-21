%
%reads the definition of a DFA from an XML file. The XML format assumed is
%the one defined by JFLAP. See: www.jflap.org
%
%input: filenameXML -> name of xml file (*.jff)
%
%urut/jan06
function [states, statesFinal, statesInitial, stateTransformations, inputSymbols, inputSymbolMapping, stateIdMapping, actionSymbols, actionSymbolMapping, actionTransitionMapping] = convertJFLAPtoMatlab(filenameXML)
xDoc = xmlread(filenameXML);
xRoot = xDoc.getDocumentElement;

%-- get the states
states=[]; %list of all used states, in internal IDs
statesFinal=[]; %which states are accepting
statesInitial=[]; % which state is the start state (only one)
stateIdMapping=[]; %internal-state-ID, external-state-ID

actionSymbols=[];
actionSymbolMapping=[];
actionTransitionMapping=[];

allListItems = xDoc.getElementsByTagName('state');


nextFreeID=0;

for k = 0:allListItems.getLength-1
   thisListItem = allListItems.item(k);
   
   currentStateNr=-1;
   theAttributes = thisListItem.getAttributes;
   numAttributes = theAttributes.getLength;
   for count = 1:numAttributes
       if theAttributes.item(count-1).getName=='id'
           
           stateNrOrig = str2num( theAttributes.item(count-1).getValue );
           
           %stateNrOrig from file is arbitrary, re-label for internal
           %purposes and store original label 
           stateNr = nextFreeID;
           states = [ states stateNr ];

           stateIdMapping(stateNr+1,1:2) = [ stateNr stateNrOrig ];
           
           currentStateNr=stateNr;
           
           nextFreeID=nextFreeID+1;           
       end
   end
   
   childNodes = thisListItem.getChildNodes;
   numChildNodes = childNodes.getLength;
   for count = 1:numChildNodes
       childItem= childNodes.item(count-1);
       if childItem.getNodeType == childItem.ELEMENT_NODE
            if childItem.getNodeName=='final'
                statesFinal=[statesFinal currentStateNr];
            end
            if childItem.getNodeName=='initial'
                statesInitial=[statesInitial currentStateNr];
            end
       end
        
   end
end

           
%-- get the state transformations
stateTransformations=[];  % from, symbol, to
inputSymbolMapping=[]; %nr of the symbol...
inputSymbols=[];       %string of the symbol (array of strings)

allListItems = xDoc.getElementsByTagName('transition');

for k = 0:allListItems.getLength-1
   thisListItem = allListItems.item(k);
      
   
   fromState=-1;
   toState=-1;
   symbol='';
   action='';
   
   childNodes = thisListItem.getChildNodes;
   numChildNodes = childNodes.getLength;
   for count = 1:numChildNodes
       childItem= childNodes.item(count-1);
       if childItem.getNodeType == childItem.ELEMENT_NODE
            if childItem.getNodeName=='from'
                fromState = str2num ( childItem.getFirstChild.getData );
            end
            if childItem.getNodeName=='to'
                toState = str2num ( childItem.getFirstChild.getData );                
            end
            if childItem.getNodeName=='read'     %which symbol
                symbol =  char( childItem.getFirstChild.getData );
            end
       end
   end
   
   %each symbol can only be a single letter. if longer, action that is associated with it
   if length(symbol)>1
       posAction = strfind( symbol, '(');
       if ~isempty(posAction)
           action=symbol(posAction+1);
       else
           action='#'; %illegal action
           warning('a transition has an illegal action - ignore');
       end
       symbol=symbol(1);
   else
       action='#'; %default action (none)
   end
   
   symbInd = find( inputSymbols==symbol );
   if isempty(symbInd)
       symbInd=length(inputSymbols)+1;
       inputSymbols = [ inputSymbols symbol ];
       inputSymbolMapping = [ inputSymbolMapping symbInd ];
   end
   
   actionInd=find( actionSymbols==action);
   
   if isempty(actionInd)
       actionInd=length(actionSymbols)+1;
       actionSymbols = [ actionSymbols action ];
       actionSymbolMapping = [ actionSymbolMapping actionInd ];       
   end
   
   %conv to internal IDs 
   fromStateConv = stateIdMapping(find( fromState == stateIdMapping(:,2) ),1);
   toStateConv = stateIdMapping(find( toState == stateIdMapping(:,2) ),1);
   
   l = size(stateTransformations,1);
   stateTransformations(l+1,:) = [fromStateConv symbInd toStateConv];   
   
   actionTransitionMapping(l+1,:) = [l+1 actionInd]; %which TN, action-ID
end
