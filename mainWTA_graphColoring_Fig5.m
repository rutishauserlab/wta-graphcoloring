%
% mainWTA_graphColoring_fig5.m
%
% reproduces figure 5 of paper (MIS graph coloring problem)
%

%% set parameters
simGrpStr = 'Fig5';   
runPrefix='Fig5';

graphType = 4; % simulated MIS graph

% How many nodes should the random graph have. only numbers 2^X are allowed
% to change to another fixed MIS, see flag % useFixedMIS_example in loadGenerateNewGraphDB.m
nrNodes = 8;  % 2^X for an arbitrary random MIS graph.  setting =8 leads to loading the special graph of Fig 4A.

%% execute (fig 5a graph)
nrNodes = 8;  % 2^X for an arbitrary random MIS graph.  setting =8 leads to loading the special graph of Fig 4A.
WTAtype = 1;
mainWTA_graphColoring_run1(runPrefix,nrNodes,simGrpStr,graphType,WTAtype,basepathGraphs);

%% execute (random graph)
nrNodes = 16;
WTAtype = 1;
mainWTA_graphColoring_run1(runPrefix,nrNodes,simGrpStr,graphType,WTAtype,basepathGraphs);