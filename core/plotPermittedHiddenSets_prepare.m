%
%
% see plotPermittedHiddenSets.m for details
%
%
% DD_matchedSet: set # at every point of time; set nr is an index into allSets_sorted
% totNrSets: how many sets in total (permitted and forbidden together)
% trOfSet_sorted: trace (divergence) of every set, same order as allSets_sorted
%
%
function [DD_matchedSet, totNrSets, allSets_sorted, trOfSet_sorted, indsPermitted_sorted, indsForbidden_sorted, eigOfSet_sorted] = plotPermittedHiddenSets_prepare(  Wtotal, loadTerms, DD )

N = size(Wtotal,1);
Ne = N-1;


Npermitted=Ne;
Nforbidden=2^Ne-Npermitted;

totNrSets = 2^Ne;

%% enumerate them
allSets=[];
indsPermitted=[];
indsForbidden=[];

trOfSet = [];
eigOfSet=[];
for k=1:(totNrSets)
   
    binStr = dec2bin(k-1,Ne);
        
    for j=1:length(binStr)
        allSets(k,j) = str2num(binStr(j));
    end
    
    if sum(allSets(k,:))==1
        indsPermitted = [ indsPermitted k];
    else
        indsForbidden = [ indsForbidden k];
    end
    
    
    % For each, calculate the divergence
    
    Jtmp = ( diag( [allSets(k,:) 1] ) * Wtotal ) + loadTerms;   %1 because the inhib neuron is always active

    trOfSet(k) = trace(Jtmp);
    
    
    eigs = eig( getSymPart(Jtmp)  );    % symmetric part, theta transformed
    %eigs = eig( (Jtmp)  );    % symmetric part, theta transformed

    eigOfSet(k) = max(real(eigs));
end

[~,Itr]=sort(trOfSet, 'descend');  % rankorder

allSets_sorted = allSets(Itr,:);
trOfSet_sorted = trOfSet(Itr);
eigOfSet_sorted = eigOfSet(Itr);

indsPermitted_sorted =[];
indsForbidden_sorted = [];

for k=1:size(allSets,1)
    
    if sum(allSets_sorted(k,:))<=1
        indsPermitted_sorted = [ indsPermitted_sorted k];
    else
        indsForbidden_sorted = [ indsForbidden_sorted k];
    end
end

%% --- next: plot set # as funct of time to see if it follows the gradient
% DD_matchedSet=[];
% for i=1:size(DD,2)
%     
%     activeSet = (DD(:,i));
%     
%     % find which subset this is
%     matchedSet=0;
%     for k=1:size(allSets,1)
%     
%         if isequal(allSets_sorted(k,:), activeSet(1:end-1)' )
%             matchedSet=k;
%             break;
%         end
%     end
%     
%     
%     DD_matchedSet(i) = matchedSet;
% end

%===TODO
%=== allSets_sorted is used in the below for usual simulation
%=== but for big network simulation, needs to be allSets ?

DD_matchedSet = convert_DD_toSetNr(DD,allSets_sorted, 1:(size(Wtotal,1)-1));   % all are excitatory except last
