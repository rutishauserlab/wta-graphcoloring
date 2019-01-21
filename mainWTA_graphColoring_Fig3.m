%
% mainWTA_graphColoring_fig3.m
%
% reproduces figure 3 of paper (4-coloring of 4-node graph)
%

%% set parameters
simGrpStr = 'Fig3';   
runPrefix='Fig3';

graphType=1; % demo graph
nrNodes=4;
WTAtype=1;

%% execute
mainWTA_graphColoring_run1(runPrefix,nrNodes,simGrpStr,graphType,WTAtype);