%
%
% for a collection of runs, generate set nr of all sets that were visited
%
function [allSets_generated, indsPermitted, indsForbidden, eigOfSet, trOfSet] = enumerateSets_dynamic( DD_all, nrRuns, indsExcit, Wtotal, loadTerms )

allSets_generated = [];


parfor k=1:nrRuns
    DD = DD_all(k).DDstats;
    
    toMatch = DD(indsExcit,:)';
    uniqueToMatch = unique( toMatch, 'rows');
    
    %if k==1
        %no sets exist yet
    %    allSets_generated = uniqueToMatch;
    %else
        allSets_generated = [allSets_generated; uniqueToMatch];
    %end
end

allSets_generated = unique( allSets_generated, 'rows');


indsPermitted=[];
indsForbidden=[];
eigOfSet=[];
trOfSet=[];

N=size(Wtotal,1);
%--- for the visited sets, evaluate whether they are hidden/forbidden
parfor k=1:size(allSets_generated,1)
   
    set = allSets_generated(k,:);
    
    
    diagTerms = ones(1,N);
    diagTerms(indsExcit(find( set==0  ) )) = 0;   % those not part of this set are off
    
    Jtmp = ( diag( diagTerms ) * Wtotal ) + loadTerms;   %1 because the inhib neuron is always active

    trOfSet(k) = trace(Jtmp);

    % theta transform if these were the winners
    
    %[V,D] = eig( Jtmp );  %eigenvectors
    %ThetaInv = V;
    %Theta = inv(V);
    
    %Jtrans = Theta * Jtmp * ThetaInv;  %apply theta transform    
    
    %eigs = eig( getSymPart(Jtrans)  );    % symmetric part, theta transform
    
    %eigs = eig( getSymPart(Jtmp)  );    % symmetric part, theta transformed
    %eigs = eig( (Jtmp)  );    % symmetric part, theta transformed
    %eigOfSet(k,:) = [max(real(eigs)) min(real(eigs))];
            
    % define whether set is forbidden or permitted based on its max eigenvalue
    %if eigOfSet(k)<0
    %    indsPermitted = [ indsPermitted k];
    %else
    %    indsForbidden = [ indsForbidden k];
    %end    

    
end