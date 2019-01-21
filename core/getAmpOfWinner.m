%
%get the amplitude of the winner on all maps at the end of the simulation
%t is the timepoint to use,if not given the last possible is used
%
%urut/jan12
function [maxAmpAll, winnerIDAll, ampInhib, ampInput] = getAmpOfWinner(Mhistory, Nnodes,t)
if nargin<3
    t = [];
end

if isempty(t)
        enableInputMean=1;
else
        enableInputMean=0;
end

maxAmpAll = zeros(1,Nnodes);  % amplitude of winner on this map
winnerIDAll = zeros(1,Nnodes);
ampInput=[];
for k=1:Nnodes
    
   S = Mhistory(k).S;
   if isempty(t)
       t = size(S,2);  %pick the last timepoint if not given
   end
   
   maxAmp = max(S(1:end-1,t));   %end-1 dont consider inhib unit
   maxAmpID = find( S(1:end-1,t) == maxAmp );
   winnerIDAll(k) = maxAmpID(1);
   maxAmpAll(k) = maxAmp;
   
   ampInhib(k) = S(end,t);   
end

%only run if no timepoint was given
if enableInputMean
    for k=1:Nnodes
        if isfield(Mhistory(k), 'inpHistory')
            ampInput(k) = mean( Mhistory(k).inpHistory(maxAmpID(1), :) );
        else
            ampInput(k)=0;
        end
    end
end