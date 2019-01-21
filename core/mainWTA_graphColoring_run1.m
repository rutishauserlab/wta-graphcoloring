%Graph coloring main function
%
%Production version for automatic runs. Forked from mainWTA_graphColoring.m (which is now outdated).
%
% graphType : 1 manual graph, 2 simulated random graph, 3 sudoku (general graph,not planar), 4 max independent set
% WTAtype: 1 standard (S), 2 is extended (E)
% disableBatchRun: 1 is single run, 0 batch run
% basepathGraphs: path for jff and dot files
%
% (c) ueli rutishauser/MPI
function mainWTA_graphColoring_run1(runPrefix, nrNodes, simGrpStr, graphType, WTAtype, disableBatchRun, basepathGraphs)
if nargin<6
    disableBatchRun=1;  % batch run is for automatic runs
end

basepathHome = getuserdir;

plotModeDynamic = 0; %plots during simulation. 0 no biograph, 1 biograph on (need bioinformatics toolbox)
plotModeNetwork = 0; %plots of network/connectivity. 0 no plots, 1 plots on 
plotModeResult = 1;  %plot the simulation result

economicMode = 1;  % on = dont keep full state history. =0 is needed for some statistics

switch(WTAtype)
    case 1
        %without dendrite
        alpha1 = 1.5;
    
        %50perc reduced. Default DB cell for without dendrite
        beta1_DB = 1.5;   % graph default is 1.65
        beta2_DB = 0.15;   % graph default is 0.2
        
    case 2
        %with dendrite
        
        if graphType==3
            alpha1 = 1.1;  %default 1.5
        else
            alpha1 = 1.2;  %default 1.5
        end
        
        %Default DB cell for with dendrite
        beta1_DB = 3.0;   
        beta2_DB = 0.3;        
end

alpha2 = 0.0;

beta1  = 3.0; %DEFAULT
beta2  = 0.3; %DEFAULT

%beta1P/2P have to be nummerically different from beta1DB/2DB (programming reasons)
%Default MIS for with dendrite
beta1P_DB = 0.8;
beta2P_DB = 0.16;

%reduced for without dend
%beta1P_DB = 1.5/4;
%beta2P_DB = 0.15/4;

%==== use to enable random biases (only graphType==2
if graphType==2
    
    determineRandBiases = 0;
   
    %=== below can be used to test what the network does with some nodes already set to correct solution
    %determineRandBiases = 1; % 0 no, 1 yes - only useful for random graphs, determines a solution so that random biases can be picked from it
    %%pick a random subset of biases to apply
    %propNodesWithBias=0.2;   %1 is all, 0 is none
else
    determineRandBiases = 0;
end

%For version without dendrite, disable
if WTAtype==1
    dendParams = [-1 0]; % slope of dendritic non-linearity; -1 means disabled
else
    %dendParams = [0.15 0]; % slope of dendritic non-linearity; -1 means disabled   testing maxindp
    if graphType==3
        dendParams = [4 4]; % slope of dendritic non-linearity; -1 means disabled   DEFAULT SUDOKU
    else
        dendParams = [0.15 0]; % slope of dendritic non-linearity; -1 means disabled   DEFAULT COLORING
    end
end

gamma = 0;  %unused for DB cells version

T_DB     = 0.0;

if graphType==2
    % for random graph with shunt, requires constant bias input
    Tdefault = 0;  %excitatory threshold;   negative value means common DC input to all excitatory (positive)
else
    Tdefault = 0;  %excitatory threshold;   negative value means common DC input to all excitatory (positive)
end

Tinhib   = 0;
dt=.01; %integration time constant

NrOfColors = 4; %normal for graph coloring

stepSizeStates = 1; %distance between states in #units

%how many integration timesteps.
%till = 500000; %how many integration timesteps.  *dt to get nr tau
till = 250000; %how many integration timesteps
%till = 60000; %how many integration timesteps
%till = 100000; %how many integration timesteps

noiseUpdateFreq = 100;

%=2 for graph is default (iid)
enableNoiseInp = 2; % 0 none, 1 frozen noise to excit units, 2 continuous noise to excit units, 3 1/f noise

switch (graphType)
    case 3
        %sudoku
        noiseMean = 4;
        noiseStd  = 1;
    case {1,2}
        if enableNoiseInp==2
            %graph coloring (default, works for iid noise)        
            noiseMean = 1.5;
            noiseStd  = 0.15;
        end
        if enableNoiseInp==3
            %1/f noise testing
            noiseMean = 3;
            noiseStd  = 0.1;
        end
        
    case 4
        %MIS, has special params
        alpha1=1.2;
        
        
        noiseMean=1.5;
        noiseStd=0.15;
    otherwise
        errror('unknown graphType');
end

enableInteract = 1; % enable input from the DBs

modeNormalizeWeights = 0; %normalize weights from/to the DB cells
modeWeightNoise = 0; %add noise to weights from/to the DB cells

normConnsEnabled=0;
normFact = abs(gamma);

addDebugInfo = 0; % if enabled, stores effective weight matrices (slow)

%always OFF. Now this is switched by using variable dendParams.
enableDendriticCompartment = 0 ; %1 yes, 0 no    Setting to know will disable the shunt. This parameter is only used if dendParams is unused (-1). otherwise,is ignored.

keepHistory = 0;
terminateIfNoError = 1;

debugPlotsEnabled=0;  % make plots of weight matrices etc; disable for automatic runs (saves time)

%seed each differently
fixedSeed=0; %default
%fixedSeed=1004444;   % to reproduce a particular simulation


seedValTotal = sum(1000*clock);

%==== load a pre-initialized state (for testing certain hypotheiss)
%enableInitState=20*100;   %0 disabled >0, fix values to provided values for first X provided integration steps
enableInitState=0;   %0 disabled >0, fix values to provided values for first X provided integration steps
initStateVals=[];   % same length as Soverall, all items >=0 are fixed, negative means free-floating
if enableInitState
    load('c:/temp/Soverall_5842_r1.mat');
    
    copyMode_Soverall = 2;  % 1 all, 2 partial (correct), 3 partial (incorrect), 4 all incorrect
    nrNodes_toCopy    = 12;
    
    switch (copyMode_Soverall)
        case 1
            initStateVals=Soverall;
        case 2
            initStateVals=-1*ones(size(Soverall,1),1);
            for k=1:nrNodes_toCopy
                % copy state values of these maps
                inds = [1:5]+(k-1)*5;
                inds = inds(1:end-1); %skip inhib
                initStateVals(inds) = Soverall(inds);
            end
        case 3
            initStateVals=-1*ones(size(Soverall,1),1);
            for k=1:nrNodes_toCopy
                % copy state values of these maps
                inds = [1:5]+(k-1)*5;
                inds = inds(1:end-1); %skip inhib
                initStateVals(inds) = [0 0.8 0 0 ]; %Soverall(inds);   set to incorrect values (all same winner)
            end
    end
end

%for simulated graphs
if ~exist('nrNodes')    
    nrNodes = 6^2;
else
    display(['using previous value of nrNodes: ' num2str(nrNodes)]);
end
if ~exist('massRun')   %non-GUI simulation,skip some stuff
    massRun=0;
end

density = 0.8;

% special sudoku settings
if graphType==3
%    NrOfColors = 4;    
%    nrNodes=16;
    NrOfColors = 9;    
    nrNodes=81;
end

%special settings for max independent set
if graphType==4
    NrOfColors=2;
end

% calculation of network properties and initialize
g=1/(1+beta1*beta2-alpha1);
gainNoDB =  g * ( 1+ beta1*Tinhib )  %gain of the entire system with no DBs
gDB = 1/(1+beta1*beta2+beta1_DB*beta2_DB-alpha1);
gainWithAllDB = gDB * ( 1 + beta1*Tinhib + beta1_DB*T_DB)

paramStr=['\alpha=' num2str(alpha1) ' \beta_1=' num2str(beta1) ' \beta_2=' num2str(beta2) ' \gamma=' num2str(gamma) ' T_e=' num2str(Tdefault) ' T_i=' num2str(Tinhib) ' TDB=' num2str(T_DB) ' \betaDB_1=' num2str(beta1_DB) ' \betaDB_2=' num2str(beta2_DB)];

RandStream.setGlobalStream(RandStream('mt19937ar','seed',seedValTotal));
simNrRand=round( rand*10000 );
disp(['running sim nr: ' num2str(simNrRand) ' graphType=' num2str(graphType) ]);

%% make new graph and simulate it
%for testing - fixed seed to make reproducable
if fixedSeed>0 & ~massRun    %only use seed if not automatic simulation
    warning(['fixed seed is set ' num2str(fixedSeed)]);
    RandStream.setGlobalStream(RandStream('mt19937ar','seed',fixedSeed));
end

% create the graph, convert it to WTA network
ensureMinOnePerNode = 0;
tic
[hGraph,Maps,Conn,Nnodes,Nedges, WPyr_to_db, WDb_to_pyr,stateTrans, colOfNodes ] = loadGenerateNewGraphDB(graphType, stepSizeStates, plotModeNetwork, plotModeDynamic, NrOfColors, alpha1,...
    beta1,beta2,beta1_DB,beta2_DB,beta2P_DB, beta1P_DB, gamma, Tdefault,Tinhib, normConnsEnabled, normFact, nrNodes,density, ensureMinOnePerNode,simNrRand,modeNormalizeWeights, modeWeightNoise, basepathGraphs, debugPlotsEnabled, determineRandBiases);
disp(['Elapsed time ']);
toc

%=== plot overall weight matrix
figure(14);
plotNetworkSummary_weightMatrix( Maps, WPyr_to_db, WDb_to_pyr, NrOfColors);

% enable/disable interactive plotting
plotUpdateFreq  = 100;
%plotUpdateFunct = @graphPlotUpdate; %callback
plotUpdateFunct = []; %disable graph plot
    

%==== Sudoku specific settings
biasInp_forSudoku=[];
pickRandomPuzzle = 0;   % 1 pick a random puzzle, 0 apply a specific puzzle
randomPuzzleType = 0;
if graphType==3
    if ~pickRandomPuzzle
        %puzzleNr=3; %easy
        %puzzleNr=16; %other hardest
        puzzleNr=15; %maas paper hard one
        %puzzleNr=14; %hardest from newspaper
    else
        randomPuzzleType=3; %1 hard95, 2 hardest, 3 50-euler
        
        randPuzzles1 = 1:50;
        R=randperm(length(randPuzzles1));
        puzzleNr = randPuzzles1(R(1));
        disp(['Puzzle nr chosen: ' num2str(puzzleNr)]);
    end
    
    [totNrBiases, biasInp_forSudoku,puzzleStr] = mainWTA_sudoku(puzzleNr,randomPuzzleType, basepathGraphs );
else
	puzzleStr='';
end

%% run simulation
switch(graphType)
    case 3
        %sudoku
        biasAmp = 10;  %keeps bias (sudoku)
    case {1,2}
        %graph coloring
        biasAmp = 0;
    case 4
        %indp set
        biasAmp= 0.0;
    otherwise
        errror('unknown graphType');
end

if graphType==2 & determineRandBiases==1
   % if random planar graph, but with biases
   disp(['Setting random biases is enabled']);

   biasAmp = 10;  % like in sudoku
 
   biasInp_forSudoku=ones(1,Nnodes)*-1;  %bias inp to each node
 
   nrBiasesToPick = floor(Nnodes*propNodesWithBias);
   
   rInds_forBiases = randperm(Nnodes);
   rInds_used = rInds_forBiases(1:nrBiasesToPick);
   biasInp_forSudoku( rInds_used ) = colOfNodes(rInds_used);   
end

[Maps,biasInps] = prepare_hardBiasInputs(Maps,Nnodes,graphType,biasInp_forSudoku, biasAmp,determineRandBiases);

debugMode=1;

TOff = zeros(NrOfColors+1,1);
TOff(end) = 0;
Mold=[];
extInpClamp=[];

% speed-optimized version
[Mhistory,DBhistory,simTill,paramsUsed_runN_WTA, nrErrorsOverTime,nrSwitches,stateHistoryAll,shuntHistoryAll] = ...
    runN_WTA_withDBs_optimized_shunt( Maps, [], [], Conn, dt, till,enableInteract, enableNoiseInp, noiseUpdateFreq,noiseStd, ...
    noiseMean, extInpClamp, plotUpdateFunct, plotUpdateFreq, hGraph, TOff, WPyr_to_db, WDb_to_pyr, T_DB, terminateIfNoError, ...
    biasInp_forSudoku, keepHistory, graphType, addDebugInfo,enableDendriticCompartment, economicMode, debugMode, ...
    dendParams, enableInitState, initStateVals  ) ;

% full version (but slow)
%[Mhistory,DBhistory,simTill,paramsUsed_runN_WTA, nrErrorsOverTime,nrSwitches] = runN_WTA_withDBs( Maps, [], [], Conn, dt, till,enableInteract, enableNoiseInp, noiseUpdateFreq,noiseStd, noiseMean, extInpClamp, plotUpdateFunct, plotUpdateFreq, hGraph, TOff, WPyr_to_db, WDb_to_pyr, T_DB, terminateIfNoError, biasInp_forSudoku, keepHistory, graphType, addDebugInfo,enableDendriticCompartment, economicMode ) ;

% evaluate correctness of the result.
% loop over all edges,check difference of colors
if graphType~=4
    [colMatchFound,errPercentage, errorNodes] = checkColorSolutionOfSim(Conn, Maps, Mhistory );
else
    [colMatchFound,errPercentage, errorNodes] = checkColorSolutionOfSim(Conn, Maps, Mhistory,0,[],2 );
end

disp(['errorFound=' num2str(colMatchFound) ' errPerc=' num2str(errPercentage)]);
if colMatchFound>0
    warning(' solution was not perfect, has errors');
end
disp(['simulation finished - steps used: ' num2str(simTill) ]);
            

%% store results and skip plots if this is a batch run
if ~disableBatchRun
    mainWTA_graphColoring_assembleParams;  % this will make a variable paramsAll
    if ~exist('runPrefix')
        runPrefix='';
    end
    
    %prepare data to export
    if ~economicMode
        [nrSwitches, nrErrorsOverTime,colMatchFound] = graphColoring_computeRunMetrics( Conn, Maps, Mhistory, Nnodes, simTill);
    else
        nrSwitches=[];
    end
    
    outputDir = [basepathHome '/resultsConstr/' simGrpStr '/'];
    if ~exist(outputDir)
        mkdir(outputDir);
    end
    
    threshActive_divAll=0.1;
    divAll = calc_divOvertime( stateHistoryAll, nrNodes, NrOfColors+1, threshActive_divAll, simTill, alpha1 );
    
    save([basepathHome '/resultsConstr/' simGrpStr '/' runPrefix '_R_' num2str(enableNoiseInp) '_' num2str(seedValTotal) '.mat'], 'paramsAll', 'paramsUsed_runN_WTA', ...
        'colMatchFound','errPercentage','puzzleStr','simTill','enableNoiseInp', 'nrSwitches', 'nrErrorsOverTime', 'colMatchFound','Nedges','dendParams', 'pickRandomPuzzle', ...
        'randomPuzzleType', 'divAll', 'threshActive_divAll');
    
    if simTill==till
        %got stuck
        %    keyboard;
    end
    
    return;
    
end

%% plot state as a function of time 
figure(200);
imagesc(stateHistoryAll(:,1:simTill));
colorbar;
unitsPerMap = NrOfColors+1;

%== plot activity on each map
cols={'r','g','b','c','m','k','r--','g--','b--'};
nrNodesToPlot=25; %max nr nodes to plot

if nrNodesToPlot>nrNodes
    nrNodesToPlot=nrNodes;
end

nrColsToPlot=4;
if graphType==4
    nrColsToPlot=2;  % MIS
end

figure(201);
for k=1:nrNodesToPlot
    indsOfUnit = [1:unitsPerMap] + (k-1)*unitsPerMap;
    
    subplot(5,5,k);
    indTmp=indsOfUnit(1:unitsPerMap-1);
    for jj=1:nrColsToPlot
        if jj>1
            hold on;
        end
        plot( stateHistoryAll(  indTmp(jj) ,1:simTill)', cols{jj} );
     %plot( stateHistoryAll(  indsOfUnit(1:unitsPerMap-1) ,1:simTill)');
    end
    hold off
    hold on
        plot( stateHistoryAll( indsOfUnit(unitsPerMap),1:simTill)', 'm--');
    hold off
    title(['Map=' num2str(k)]);
end


%% divergence over time
threshActive=0.02;

divAll = calc_divOvertime( stateHistoryAll, nrNodes, NrOfColors+1, threshActive, simTill-200, alpha1 );
[mDiv, sDiv, seDiv, nDiv]=calcMeanSEOfSample(divAll);

figure(700);
subplot(2,2,1);
imagesc(divAll);
colorbar
ylabel('WTA #');
xlabel('time');
title([simGrpStr ' N=' num2str(nrNodes) ' thr=' num2str(threshActive)]);

subplot(2,2,2);
%plot( log(1:length(mDiv)), mDiv );
plot( log(1:length(mDiv)), mDiv );
xlabel('log time');
ylabel('average div');

subplot(2,2,3);
plot( mDiv );
ylabel('average div');
xlabel('time');

% subplot(2,2,4);
% tAll=1:length(mDiv);
% indsToPlot=tAll(end-7000):tAll(end);
% plot( tAll(indsToPlot), mDiv(indsToPlot));
% ylabel('average div');
% xlabel('time (zoom-in)');


%% plots to illustrate dendritic input, non-linearity g(z)
if WTAtype~=1
    plot_dendriticShunt_ofSinglerun(unitsPerMap, nrColsToPlot, simTill, dendParams, shuntHistoryAll,stateHistoryAll,nrNodesToPlot);
end

%% verify the max indp set
if ~massRun & graphType==4
    [colMatchFound,errPercentage, errorNodes] = checkColorSolutionOfSim(Conn, Maps, Mhistory,0,[],2 );
    if colMatchFound
        disp(['max independent set has errors : ' num2str(colMatchFound)]);
    else
        disp('max independent set is correct');
    end        
end        

%% sudoku result
if ~massRun & graphType==3
    
    figNr=121;
    figure(figNr);
    close;
    figure(figNr);
    
    % plot figure of solved sudoku
    graphPlot_sudokuResult(simTill, NrOfColors, Nnodes,Mhistory, paramStr, Maps, errorNodes, biasInp_forSudoku);    
end

%% export solution as jff/dot file and layout using graphviz ("dot")
if ~massRun & graphType~=3   %dont draw sudoku (graph is fixed)
    [~,nodeWinners] = getAmpOfWinner(Mhistory, Nnodes);
    exportSolvedGraphXML([basepathGraphs '/exportSolved_' num2str(simNrRand)  '.jff'], Nnodes,Nedges,stateTrans, nodeWinners);

    figure(3);
    plotChrobakDrawing_solvedGraph( stateTrans, Nedges, Nnodes, nodeWinners, errorNodes,density, colMatchFound );

    fnameDot=[basepathGraphs '/solvedGraph_N' num2str(Nnodes) '_s' num2str(simNrRand) '.dot'];
    fnameLayouted = [basepathGraphs '/layoutedGraph_N' num2str(Nnodes) '_s' num2str(simNrRand) '.pdf'];
    runMode=1;
    
    exportGraph_dotFile( stateTrans, Nedges, Nnodes, errorNodes, nodeWinners, runMode, fnameDot, fnameLayouted );
end


