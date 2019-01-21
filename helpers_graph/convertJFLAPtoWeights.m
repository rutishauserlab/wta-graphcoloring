%reads the definition of a DFA from an XML file, returns values that are
%later used for construction of weight matrices. 
%
%stepSize -> nr units that represent each state on the WTA
%
%returns:
%Npointer: nr of pointer neurons (=# transitions)
%stateTrans: state transformation matrix
%stateMapping: mapping of states onto map
%pointerMapping: which symbol belongs to which pointer neuron(s)
%stateIdMapping: mapping of the IDs used for the states in matlab (first
%column) to IDs used in Java files (second column)
%
%Nmap: size of the recurrent map
%
%urut/jan06
%urut/feb10: revised to incorporate actions associated with transitions.
%
function [Npointer, stateTrans, stateMapping, pointerMapping, inputSymbols,Nmap, stateIdMapping, statesInitial, actionMapping, actionSymbols] = convertJFLAPtoWeights(filenameXML, stepSize)
if nargin<2
    stepSize=5;
end

%-- read the xml file
[states, statesFinal, statesInitial, stateTransOrig, inputSymbols, inputSymbolMapping, stateIdMapping,actionSymbols, actionSymbolMapping, actionTransitionMapping] = convertJFLAPtoMatlab(filenameXML);

%-- convert this information to weights
Npointer = size(stateTransOrig,1);

stateTrans = stateTransOrig;
% f(A,x) -> B
% f(B,x) -> A

%-- map the states to the map
nrStates = length( states );
%stepSize = 5; %floor( (20-4) / nrStates );
stateMapping=[];
for i=1:nrStates
    stateMapping(i,1:2) = [states(i) stepSize+(i-1)*stepSize ] ; %= [ 1 4; 2 12 ];   % A (1) is 4; A ( 2) is 12 on the map
end

%how big is the recurrent map
Nmap = (i+1)*stepSize;

%-- map the input symbols onto the pointer neurons
pointerMapping=[];
nrSymbols = length(inputSymbolMapping);
for i=1:nrSymbols
    pointerMapping{i} = find( stateTrans(:,2)==inputSymbolMapping(i) );
end

%-- map the action symbols onto the pointer neurons
% (each action is the sum of a number of pointer neurons)
actionMapping=[];
nrActionSymbols = length(actionSymbolMapping);
for i=1:nrActionSymbols
    actionMapping{i} = actionTransitionMapping( find( actionTransitionMapping(:,2) == actionSymbolMapping(i)), 1);
end