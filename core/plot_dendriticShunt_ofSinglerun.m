%
% plot properties and statistics of dendritic shunt g(z) of a single
% simulation run
%
%urut/aug15
function plot_dendriticShunt_ofSinglerun(unitsPerMap, nrColsToPlot, simTill, dendParams, shuntHistoryAll, stateHistoryAll,nrNodesToPlot)

plotModeShunt=2; % 1 all together, 2 winner/looser separately

%% == plot dendritic input to each map (to inspect g(z) term, which reflects overall dendritic excitatory 
%input from apical dendrites)
%
figure(202);
for k=1:nrNodesToPlot
    
    indsOfUnit = [1:unitsPerMap] + (k-1)*unitsPerMap;
    
    subplot(5,5,k);
    
    indTmp=indsOfUnit(1:unitsPerMap-1);
    for jj=1:nrColsToPlot
        if jj>1
            hold on;
        end
        plot( shuntHistoryAll(  indTmp(jj) ,1:simTill)', rotatingColorCode(jj) );
    
    end
    hold off
%    hold on
%    plot( shuntHistoryAll( indsOfUnit(unitsPerMap),1:simTill)', 'm--');
%    hold off
    title(['Shunt Map=' num2str(k)]);
    
    if k==25
        break;
    end
end


%% ==== histogram of g(z) values at certain periods of time, across all dendrites

timestepsToPlot=[200 500 2000 simTill-100];

figure(205);
for k=1:length(timestepsToPlot)
    allValsShunt_winner=[];
    allValsShunt_looser=[];
    
    timestep = timestepsToPlot(k);
    
    for kk=1:nrNodesToPlot
        indsOfUnit = [1:unitsPerMap] + (kk-1)*unitsPerMap;

        
        % which unit is the winner/looser on this map?
        
        excitatoryUnitInds = indsOfUnit(1:end-1);
        [~,maxInd] = max( stateHistoryAll( excitatoryUnitInds , timestep));
        winnerInd = excitatoryUnitInds(maxInd);
        looserInds = setdiff(excitatoryUnitInds, winnerInd);
        
        
        allValsShunt_winner = [ allValsShunt_winner shuntHistoryAll(  winnerInd ,timestep)];
        allValsShunt_looser = [ allValsShunt_looser shuntHistoryAll(  looserInds ,timestep)];        
    end
    
    
    subplot(2,2,k)
    
    if plotModeShunt==1
        %plot all together (jointly)
        allValsShunt = [ allValsShunt_winner(:); allValsShunt_looser(:) ];
        hist(allValsShunt(:),40);
    else
        % plot winners/loosers separately
        [n,x]=hist(allValsShunt_looser(:),20);
        bar(x,n,'r');
        hold on
        [n,x]=hist(allValsShunt_winner(:),20);
        bar(x,n,'g');
        hold off
    end
    
    xlim([0 1]);
    
    if plotModeShunt==2
        winners_nonZero = allValsShunt_winner(:);
        winners_nonZero = winners_nonZero(find(winners_nonZero>0.1));
        loosers_nonZero = allValsShunt_looser(:);
        loosers_nonZero = loosers_nonZero(find(loosers_nonZero>0.1));
        
        title(['slope=' num2str(dendParams(1)) ' t=' num2str(timestep) ' mW=' num2str(mean(winners_nonZero)) ' mL=' num2str(mean(loosers_nonZero))]);
    else
        
        
        percEliminated = length(find(allValsShunt<0.02))/length(allValsShunt);
        
        title(['slope=' num2str(dendParams(1)) ' t=' num2str(timestep) ' m=' num2str(mean(allValsShunt(:))) ' %out=' num2str(percEliminated) ]);
        
    end
    
    xlabel('dendritic sensitivity g(z)');
    ylabel('nr of branches');
end

