%
%
function [hGraph,Maps,Conn,Nnodes,Nedges] = loadGenerateNewGraph(graphType, stepSizeStates, plotModeNetwork, plotModeDynamic, NrOfColors, alpha1,...
    beta1,beta2,gamma, Tdefault,Tinhib, normConnsEnabled, normFact, nrNodes,density, simNrRand)

%% load or simulate a graph
if graphType == 2
    fname=['graph/exportRand_' num2str(simNrRand) '.jff']; %load the simulated graph
    exportRandomPlanarGraph(nrNodes,density, fname);
else
    fname='graph/graph3.jff'; % for examples in draft paper
    %fname='graph/graph_twoNodes.jff';
end

%% import the graph from XML, convert it to a state machine
disp(['using graph file: ' fname ]);
[Npointer, stateTrans, stateMapping, pointerMapping, inputSymbols, Nmap, stateIdMapping,statesInitial,actionMapping, actionSymbols] = convertJFLAPtoWeights(fname, stepSizeStates );
Nnodes = size(stateMapping,1);
Nedges = size( stateTrans,1);

%% setup the WTAs

%was before CC12
%[Maps,Conn,W,Wconn,nrEdgesEach] = graphSim_setupWTAs_network(NrOfColors, alpha1,beta1,beta2,gamma,stateMapping,stateTrans, Tdefault,Tinhib,Nnodes,Nedges, normConnsEnabled, normFact);

%using new DB cells

[Maps, W, WPyr_to_db, WDb_to_pyr] = graphSim_setupWTAs_network_DB(NrOfColors, alpha1,beta1,beta2, stateMapping,stateTrans, Tdefault,Tinhib,Nnodes );

%[Maps,Conn,W,Wconn,nrEdgesEach] = graphSim_setupWTAs_network_DB(NrOfColors, alpha1,beta1,beta2,gamma,stateMapping,stateTrans, Tdefault,Tinhib,Nnodes,Nedges, normConnsEnabled, normFact);
%[ mean(nrEdgesEach) min(nrEdgesEach) max(nrEdgesEach) ] 

%% plot the weight matrices and graph
if plotModeNetwork==1
    graphPlotMatrices(W,Wconn);
end

%prepare plotting of graph
if plotModeDynamic==1
    [hGraph] = graphPlot_prepare_biograph(stateTrans,Nnodes)
else
    hGraph=[];
end
