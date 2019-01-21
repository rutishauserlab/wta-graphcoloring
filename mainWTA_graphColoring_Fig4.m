%
% mainWTA_graphColoring_fig4.m
%
% reproduces figure 4 of paper (4-coloring of random planar graph)
%

%% set parameters
simGrpStr = 'Fig4';   
runPrefix='Fig4';

graphType=2; % simulated random graph

% How many nodes should the random graph have. only numbers X^2 are allowed

%% execute (WTA_S version, 16 nodes); easier to play with
nrNodes = 16;  
WTAtype = 1;
mainWTA_graphColoring_run1(runPrefix,nrNodes,simGrpStr,graphType,WTAtype,basepathGraphs);

%% execute (WTA_S version, 49 nodes)
nrNodes = 49;   % like Fig 3A
WTAtype = 1;
mainWTA_graphColoring_run1(runPrefix,nrNodes,simGrpStr,graphType,WTAtype,basepathGraphs);

%% execute (WTA_E version)
WTAtype = 2;
mainWTA_graphColoring_run1(runPrefix,nrNodes,simGrpStr,graphType,WTAtype,basepathGraphs);
