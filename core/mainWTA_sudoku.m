%
%
% mainWTA_sudoku
%
% this prepares a sudoku as a graph,so it can later be simulated in the usual graph coloring code simulator
%
%
%randomPuzzleType: 0 not random (use puzzleNr), >0 random subsypes
%
%
function [totNrBiases, biasInp_forSudoku,puzzleStr] = mainWTA_sudoku(puzzleNr,randomPuzzleType, basepathGraphs)

if nargin<3
    basepathGraphs = '~/svnwork/DFA/matlab/graph/graphFiles'; %path for jff and dot files
end

pSize=9;
boxsize=3;
simNr=1;

%pSize=4;
%boxsize=2;
%simNr=2;

if ~exist('puzzleNr')
    puzzleNr=3;
end

buildSudokuGraph_export(basepathGraphs, simNr, pSize, boxsize);

switch(randomPuzzleType)
    case 0   % not random
        switch( puzzleNr)
            case 1
                puzzleStr='.94...13..............76..2.8..1.....32.........2...6.....5.4.......8..7..63.4..8'; %hard
            case 2
                puzzleStr='.9..71.3.87.....49..3...1..7..3.9...1.......3...4.8..5..8...2..52.....17.3.12..9.'; %easy WORKS
            case 3
                puzzleStr  ='1..8..65....91..2..8..5.7.9.......9..53.4.17..4.......5.2.9..3..9..75....76..2..5'; %easy http://www.sudoku.ws/easy-3.htm  WORKS
            case 4
                puzzleStr = '..89....17..1.3.98..3....6.....481...17...64...521.....2....8..68.5.1..75....63..'; %easy http://www.sudoku.ws/easy-4.htm WORKS (easier)
            case 5
                puzzleStr = '.6...32...1.9..6.4....8...5...8..74.9..3.2..6.73..4...3...5....2.7..9.6...67...9.'; % standard http://www.sudoku.ws/standard-2.htm WORKS
                
            case 6
                puzzleStr  ='...2...633....54.1..1..398........9....538....3........263..5..5.37....847...1...'; %hard http://www.sudoku.ws/hard-1.htm WORKS
            case 7
                puzzleStr ='.1...4.....68.5..15.37.19..8.4..7...............3..6.9..15.82.46..4.31.....2...5.'; % http://www.sudoku.ws/hard-2.htm
            case 8
                
                puzzleStr ='.241..6...5..4...8.........51...83...3.9.5.8...63...97.........9...5..7...1..743.';%    hard11
            case 9
                puzzleStr ='.1.7.6...7.3.......2.54....34....25.5.......1.81....74....57.8.......3.9...81..4.';%    hard18
            case 10
                puzzleStr ='..6.....4...86.73..4.35...217.4..6...9.....8...8..6.172...81.4..67.43...8.....3..';%   http://www.sudoku.ws/expert-1.htm
            case 11
                puzzleStr ='7.9.....8...19...2.2.8...9..7.....432.4...5.998.....7..3...5.6.8...23...5.....3.7';  %http://www.sudoku.ws/expert-2.htm
            case 12
                puzzleStr ='..9748...7.........2.1.9.....7...24..64.1.59..98...3.....8.3.2.........6...2759..'; %http://www.sudoku.ws/extreme-1.htm
            case 13
                puzzleStr ='...3.8.7.3..71...46...4....1.....63.2.6...5.8.53.....7....8...17...64..5.1.2.7...'; %http://www.sudoku.ws/extreme-2.htm
                
                
            case 14
                puzzleStr ='85...24..72......9..4.........1.7..23.5...9...4...........8..7..17..........36.4.'; %hardest http://usatoday30.usatoday.com/news/offbeat/2006-11-06-sudoku_x.htm
                
            case 15
                puzzleStr ='8.5....3..3.9.....4.6.3....6...1.9...5.3.8.7...9.4...1....2.3.8.....9.2..7....5.4'; % hard sudoku in 2013 maas paper
                
            case 16
                puzzleStr ='..53.....8......2..7..1.5..4....53...1..7...6..32...8..6.5....9..4....3......97..'; %hardest by Arto Inkala http://www.mirror.co.uk/news/weird-news/worlds-hardest-sudoku-can-you-242294
                
                
            case 17
                puzzleStr='4.....8.5.3..........7......2.....6.....8.4......1.......6.3.7.5..2.....1.4......'; %top95,line 1
        end
        
        
    case 1
        D = textread([basepathGraphs 'top95.txt'],'%s','delimiter','\n','whitespace','');
        puzzleStr = D{puzzleNr};
    case 2
        D = textread([basepathGraphs 'hardest.txt'],'%s','delimiter','\n','whitespace','');
        puzzleStr = D{puzzleNr};
    case 3
        %https://projecteuler.net/problem=96
        %50 sudokus of variable difficulty
        D = textread([basepathGraphs 'p096_sudoku.txt'],'%s','delimiter','\n','whitespace','');
        puzzleStr = D{puzzleNr};
    otherwise
        error('unknown');
end

%= verify string is proper
if length(puzzleStr)~=81
    error(['Not a valid puzzleStr ' num2str(puzzleNr-100)]);
end
disp(['PuzzleStr=' puzzleStr]);

%= convert puzzleStr to biasInp vector
biasInp_forSudoku=ones(1,pSize^2)*-1;  %bias inp to each node
totNrBiases=0;
for k=1:length(puzzleStr)
    
    if puzzleStr(k) == '.'
        biasInp_forSudoku(k)=-1;
    else
        biasInp_forSudoku(k) = str2num(puzzleStr(k));
        totNrBiases=totNrBiases+1;
    end
end


