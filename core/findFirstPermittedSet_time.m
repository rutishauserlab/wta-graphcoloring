%
% finds the first point of time a permitted set was entered
%
%urut/CC2014
function firstWinnerInd = findFirstPermittedSet_time( DD_matchedSet, indsPermitted_sorted)
%find first time a permitted set was reached
isWinner = zeros(1,length(DD_matchedSet));
for kk=1:length(indsPermitted_sorted)
   isWinner(find(DD_matchedSet==indsPermitted_sorted(kk)))=1; 
end
indsWinners=find(isWinner==1);
firstWinnerInd = indsWinners(1);