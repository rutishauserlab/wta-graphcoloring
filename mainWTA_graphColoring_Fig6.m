%
% mainWTA_graphColoring_fig6.m
%
% reproduces figure 6 of paper (sudoku graph coloring problem)
%

%% set parameters
simGrpStr = 'Fig6';   
runPrefix='Fig6';

graphType = 3; % sudoku

nrNodes = 81;  % this is fixed for sudoku
WTAtype = 2;   % sudoku uses WTA_E version

%see pickRandomPuzzle and puzzleNr in mainWTA_graphColoring_run1.m to set the puzzle that is used.
%by default it is the one in Fig 6A.
%
%Also, note that the simulation execution time used by default might not be enough for all puzzles. increase till=.... to change this.

%% execute 
mainWTA_graphColoring_run1(runPrefix,nrNodes,simGrpStr,graphType,WTAtype, basepathGraphs);