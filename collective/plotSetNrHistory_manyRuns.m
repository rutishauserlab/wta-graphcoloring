%
%
%
%urut/june14
function plotSetNrHistory_manyRuns(setNrAll, setNrsToUse, endPointToUse, subplotX,subplotY, plotPos, allSets_sorted, indsPermitted_sorted, indsForbidden_sorted)

t = 1:size(setNrAll,2);

totNrSets = size(allSets_sorted,1);

c=0;
for ii=1:length(setNrsToUse)
    
    selFor = setNrsToUse(ii);
    
    toPlot = selectRuns_byEndpoint_reached(setNrAll, endPointToUse, selFor);

    %===== history
    subplot(subplotX,subplotY, plotPos(1)+ii-1);
    
    %mark hidden/permitted sets with lines
    
    for setInd=1:length(indsPermitted_sorted)
        setNr = (indsPermitted_sorted(setInd));
        hLine = line([t(1) t(end)],[setNr setNr], 'color', 'g', 'LineStyle', '--');
    end
    for setInd=1:length(indsForbidden_sorted)
        setNr = (indsForbidden_sorted(setInd));
        hLine = line([t(1) t(end)],[setNr setNr], 'color', 'r', 'LineStyle', '--');
    end
    
   % allSets_sorted, indsPermitted_sorted, indsForbidden_sorted
    
    hold on
    
    if length(toPlot)>50
        %random subset
        R=randperm(length(toPlot));
        toPlotUse = toPlot(R(1:50));
    else
        toPlotUse=toPlot;
    end
    
    plot( t, setNrAll(toPlotUse,:)', 'k');
    hold off
    title(['endpoint set #' num2str(selFor)]);
    ylim([0 totNrSets+1]);
    ylabel('set #');
    xlabel('time');
    hold off
    
    
    
    %==== histogram
    subplot(subplotX,subplotY, plotPos(2)+ii-1);
    allVisited=[];
    for jj=1:length(toPlot)
        allVisited = [ allVisited unique(setNrAll(toPlot(jj),:))];
    end
    
    edges=[-0.5:1:totNrSets+1];
    [N]=histc( allVisited(:), edges );
    bar(edges,N,'histc');
    
    xlim([0 totNrSets+1]);
    
    title(['sets visited, endpoint #' num2str(selFor) ]);
    xlabel('set #');
    ylabel('# of visits');
end
