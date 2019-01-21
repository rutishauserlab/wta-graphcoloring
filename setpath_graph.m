
%% change the following accordingly
basepath2='/home/urut/svnwork/DFA/matlab/release_sudoku/';

%% no changes needed below
basepathGraphs = [basepath2 '/graphFiles/'];

path(path,[basepath2 '3rdParty/']);
path(path,[basepath2 '3rdParty/matlab_bgl/']);
path(path,[basepath2 '3rdParty/GraphViz2Mat/']);

path(path,[basepath2 'collective/']);   % this is code from the collective computation paper (Rutishauser et al 2011).

path(path,[basepath2 'core/']);
path(path,[basepath2 'helpers_graph/']);
path(path,[basepath2 'helpers_general/']);