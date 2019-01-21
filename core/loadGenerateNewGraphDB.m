%
%
function [hGraph,Maps,Conn,Nnodes,Nedges,WPyr_to_db, WDb_to_pyr,stateTrans,colOfNodes] = loadGenerateNewGraphDB(graphType, stepSizeStates, plotModeNetwork, plotModeDynamic, NrOfColors, alpha1,...
    beta1,beta2, beta1_DB, beta2_DB, beta2P_DB, beta1P_DB, gamma, Tdefault,Tinhib, normConnsEnabled, normFact, nrNodes,density, ensureMinOnePerNode, simNrRand,modeNormalizeWeights, modeWeightNoise, basepathGraphs, debugPlotsEnabled, determineRandBiases)
Conn=[];
colOfNodes=[];

%% load or simulate a graph
switch(graphType)
    case 1  % predefined graph
        %fname = [basepathGraphs '/graph_threeNodesOnly.jff']; % three nodes,impossible to solve with 2 colors
        %fname = [basepathGraphs '/graph_twoNodesOnly.jff']; % two nodes only for demo
        fname = [basepathGraphs '/graph3.jff']; % simple example for paper
        %fname = [basepathGraphs '/graph9_withBoxes.jff']; % for examples in draft paper
        %fname='graph/graph_twoNodesOnly.jff.jff';
    case 2 %random planar graph
        
        useFixedGraph_forTesting=0;
        
        if ~useFixedGraph_forTesting
            fname = [basepathGraphs '/exportRand_' num2str(simNrRand) '.jff']; %load the simulated graph
            colOfNodes = exportRandomPlanarGraph(nrNodes,density, ensureMinOnePerNode,fname, determineRandBiases);
            
        else
            %=== For testing, use a fixed random graph
            %simNrRand=5570;
            %simNrRand=6451;
            %simNrRand=7523;
            %simNrRand=4945;
            simNrRand=5842;
            fname = [basepathGraphs '/exportRand_' num2str(simNrRand) '.jff']; %load the simulated graph
            warning(['Using fixed rand graph: ' fname]);
        end
    case 3
        fname = [basepathGraphs '/sudoku_s1.jff']; %load the simulated graph           %s1 is 81 nodes, s2 16 nodes
    case 4    % MIS problem; random graph with only 2 colors
        
        useFixedMIS_example = 0;  % 1 is a predefined example, 0 a random generated one
        
        if nrNodes==8   % special value to request fixed graph
            useFixedMIS_example=1;
        end
        
        if useFixedMIS_example
            fname = [basepathGraphs '/maxIndpSet_Fig5A.jff']; % for examples in paper
        else
            %generate a random example
            fname = [basepathGraphs '/exportRand_' num2str(simNrRand) '.jff']; %load the simulated graph
            exportRandomPlanarGraph(nrNodes,density, ensureMinOnePerNode,fname);
        end
      
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

[Maps, W, WPyr_to_db, WDb_to_pyr] = graphSim_setupWTAs_network_DB(NrOfColors, alpha1,beta1,beta2, stateMapping,stateTrans, Tdefault,Tinhib,Nnodes, beta2_DB, beta1_DB,beta2P_DB, beta1P_DB,modeNormalizeWeights, modeWeightNoise, graphType, debugPlotsEnabled );


%still return the Conn structure here, even if it isnt used. it is needed to verify the correctness of 
[~,Conn] = graphSim_setupWTAs_network(NrOfColors, alpha1,beta1,beta2,gamma,stateMapping,stateTrans, Tdefault,Tinhib,Nnodes,Nedges, normConnsEnabled, normFact);
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
