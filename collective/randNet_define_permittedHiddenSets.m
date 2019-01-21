%
%
% enumerate all possible sets for a network
%
%
function [allSets, indsPermitted, indsForbidden, eigOfSet, trOfSet, allSets_sorted, indsPermitted_sorted,indsForbidden_sorted,eigOfSet_sorted, trOfSet_sorted] = randNet_define_permittedHiddenSets( Wtotal, indsInhib, indsExcit, loadTerms )


Ne = length(indsExcit);
Ni = length(indsInhib);

totNrSets = 2^Ne;

if totNrSets>2^20
    error('too many sets, cannot enumerate');
end

%% enumerate them
allSets=zeros(totNrSets,Ne);
indsPermitted=[];
indsForbidden=[];

trOfSet = [];
eigOfSet=[];
for k=1:(totNrSets)
   
    binStr = dec2bin(k-1,Ne);
        
    for j=1:length(binStr)
        allSets(k,j) = str2num(binStr(j));
    end
    
    
    % For each, calculate the divergence
    
    diagTerms = ones(1,Ne+Ni);
    diagTerms(indsExcit(find( allSets(k,:)==0  ) )) = 0;   % those not part of this set are off
    
    Jtmp = ( diag( diagTerms ) * Wtotal ) + loadTerms;   %1 because the inhib neuron is always active

    trOfSet(k) = trace(Jtmp);

    % theta transform if these were the winners
    
    [V,D] = eig( Jtmp );  %eigenvectors
    ThetaInv = V;
    Theta = inv(V);
    
    Jtrans = Theta * Jtmp * ThetaInv;  %apply theta transform    
    
    eigs = eig( getSymPart(Jtrans)  );    % symmetric part, theta transform
    
    %eigs = eig( getSymPart(Jtmp)  );    % symmetric part, theta transformed
    %eigs = eig( (Jtmp)  );    % symmetric part, theta transformed
    eigOfSet(k) = max(real(eigs));
        
    % define whether set is forbidden or permitted based on its max eigenvalue
    if eigOfSet(k)<0
        indsPermitted = [ indsPermitted k];
    else
        indsForbidden = [ indsForbidden k];
    end    
end

%% reorder the sets according to permitted/forbidden
[~,Itr]=sort(trOfSet, 'descend');  % rankorder

allSets_sorted = allSets(Itr,:);
trOfSet_sorted = trOfSet(Itr);
eigOfSet_sorted = eigOfSet(Itr);

indsPermitted_sorted =[];
indsForbidden_sorted = [];
for k=1:size(allSets_sorted,1)    
    if eigOfSet_sorted(k)<0
        indsPermitted_sorted = [ indsPermitted_sorted k];
    else
        indsForbidden_sorted = [ indsForbidden_sorted k];
    end
end



