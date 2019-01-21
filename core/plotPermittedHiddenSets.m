% permitted/hidden sets analysis plot
% shows the set # as a function of time together with its divergence, ordered by divergence to illustrate the dynamics of computation in a WTA node
%
%input parameters: Wtotal: weight matrix (full), DD active subset at every point of time (switching matrix)
%
%plotMode: 0 no, 1 yes (all), 2 yes (only add additional state space trajectory, 3 only set nrs vs divergence
%
%written on LH450 june 2nd in-air, urut
%
function [allSets_sorted,trOfSet_sorted, firstWinnerInd] = plotPermittedHiddenSets( Wtotal, loadTerms, DD, plotMode, transitionPairs )
if nargin<4
    transitionPairs=[];
end

%%
[DD_matchedSet, totNrSets, allSets_sorted, trOfSet_sorted, indsPermitted_sorted, indsForbidden_sorted] = plotPermittedHiddenSets_prepare(  Wtotal, loadTerms, DD );
firstWinnerInd = findFirstPermittedSet_time( DD_matchedSet, indsPermitted_sorted);

if plotMode  % if ==0, no plotting
    
    x_DD=1:size(DD_matchedSet,2);
    indsUse = find(DD_matchedSet<=totNrSets);
    
    plotOrderInds = [indsForbidden_sorted indsPermitted_sorted]; %index into allSets_sorted_plot
    %relate new plotting pos to original set nr; DD_matchedSet refers to original order
    mapToNewOrderInds=[];
    for kk=1:size(allSets_sorted,1)
        mapToNewOrderInds(kk) = find( plotOrderInds==kk);
    end
    
    %figure(81);
    
    switch(plotMode)
        case 1
            
            subplot(2,2,1);
            
            toPlot=allSets_sorted( indsPermitted_sorted,:)';
            
            %plot boxes to illustrate the sets
            plotBoxed_hiddenPermitted( toPlot, {'w','g'} );
            
            %imagesc(  );
            
            title('permitted');
            
            set(gca,'XTick', [1:(size(Wtotal,1)-1)]-0.5  );
            set(gca,'XTickLabel', [1:(size(Wtotal,1)-1)] );
            
            set(gca,'YTick', [1:(size(Wtotal,1)-1)]-0.5  );
            set(gca,'YTickLabel', [1:(size(Wtotal,1)-1)] );
            
            subplot(2,2,2);
            
            allSets_sorted_plot = allSets_sorted;
            allSets_sorted_plot(indsPermitted_sorted,:) = allSets_sorted_plot(indsPermitted_sorted,:)*2;
            
            %imagesc( allSets_sorted_plot( [indsForbidden_sorted indsPermitted_sorted],:)' );
            
            toPlot=allSets_sorted_plot( plotOrderInds,:)';
            
            plotBoxed_hiddenPermitted( toPlot,{'w','r','g'}  );
            
            %colorbar;
            set(gca,'XTick', [1:size(allSets_sorted_plot,1)]-0.5);
            set(gca,'YTick', [1:(size(Wtotal,1)-1)]-0.5);
            set(gca,'XTickLabel', [1:size(allSets_sorted_plot,1)]);
            set(gca,'YTickLabel', [1:(size(Wtotal,1)-1)]);
            
            title('sets, rank-ordered by div (red=permitted,green=forbidden)');
            xlabel('set # (rank-ordered by tr)');
            ylabel('unit #');
            
            subplot(2,2,4);
            %hold on
            trVals_toPlot_forbidden = trOfSet_sorted( [indsForbidden_sorted ] );
            trVals_toPlot_permitted = trOfSet_sorted( [indsPermitted_sorted ] );
            
            trVals_Inds = plotOrderInds;  % mapping to set # plotted
            
            h1=plot( 1:length(trVals_toPlot_forbidden), trVals_toPlot_forbidden, 'or', 'MarkerSize', 10 , 'MarkerFaceColor', 'r');
            hold on
            xValsSetNr = [1:length(trVals_toPlot_permitted)]+length(trVals_toPlot_forbidden);
            h2=plot( xValsSetNr, trVals_toPlot_permitted, 'og', 'MarkerSize', 10 , 'MarkerFaceColor', 'g' );
            
            
            % ---- plot trajectory taken through state space
            
            xPos_setNr = mapToNewOrderInds(DD_matchedSet(indsUse));
            yPos_setNr = trOfSet_sorted((DD_matchedSet(indsUse)));
            
            %rOffset=rand(1,length(yPos_setNr))/2;  % rand offset to make distinguishable
            
            h3=plot( xPos_setNr, yPos_setNr, 'kx-');
            dXPos = sign(diff(xPos_setNr));
            dYPos = sign(diff(yPos_setNr))/4;
            scaleFact=4;
            h4=quiver( xPos_setNr(1:end-1), yPos_setNr(1:end-1), dXPos, dYPos,scaleFact,'k');
            
            
            set(h4,'LineWidth',3);
            hold off
            
            set(gca,'XTick',1:max(xValsSetNr));
            xlabel('set # (rank-ordered by tr)');
            ylabel('divergence');
            
            legend([h1 h2 h3],{'forbidden','permitted','trajectory'});
            
            
            
            subplot(2,2,2);
            hold on
            %plot( trVals_Inds(DD_matchedSet(indsUse)), 1, 'kx');
            
            hold off
            
            subplot(2,2,3);
            %plot (x_DD(indsUse), mapToNewOrderInds(DD_matchedSet(indsUse)) );
            
            ylabel('actiated subset # (re-sorted)');
        case 2
            
            
            
            %only add state space traj
            subplot(2,2,4)
            
            trVals_toPlot_forbidden = trOfSet_sorted( [indsForbidden_sorted ] );
            trVals_toPlot_permitted = trOfSet_sorted( [indsPermitted_sorted ] );
            
            trVals_Inds = [ indsForbidden_sorted indsPermitted_sorted ]; % mapping to set # plotted
            
            % ---- plot trajectory taken through state space
            xPos_setNr = mapToNewOrderInds(DD_matchedSet(indsUse));
            yPos_setNr = trOfSet_sorted((DD_matchedSet(indsUse)));
            
            %rOffset=rand(1,length(yPos_setNr))/2;  % rand offset to make distinguishable
            
            hold on
            h3=plot( xPos_setNr, yPos_setNr, 'kx-');
            dXPos = sign(diff(xPos_setNr));
            dYPos = sign(diff(yPos_setNr))/4;
            scaleFact=4;
            h4=quiver( xPos_setNr(1:end-1), yPos_setNr(1:end-1), dXPos, dYPos,scaleFact,'r');
            hold off
            
            set(h4,'LineWidth',3);
            
        case 3
            % only set numbers vs div

            trOfSet_sorted = round(trOfSet_sorted*1000)/1000;   % eliminate small numerical differences
            
            trVals_toPlot_forbidden = trOfSet_sorted( [indsForbidden_sorted ] );
            trVals_toPlot_permitted = trOfSet_sorted( [indsPermitted_sorted ] );
            
            trVals_Inds = plotOrderInds;  % mapping to set # plotted

            
            
            
            if ~isempty( transitionPairs )
                
                markerSizeOffset=2;
                markerSizeGain=100;

                lineWidthGain=10;
                lineWidthOffset=0;

                
                
                trValsUnique = sort(unique(trOfSet_sorted),'descend');

                %assign levels to each unique tr value
                
                useStepOffsets = 1;
                
                plotProbabilities = 1; %0 no, only plot lines; 1 yes modulate dot size/line thickness by probability
                
                
                
                setNrsRemappedOffset=zeros(1,length(trOfSet_sorted));
                if useStepOffsets
                    %for each set number 1...N, list an offset which depends on how many entries have a tr value bigger than the one of this set
                    for k=1:length(trOfSet_sorted)
                        setNrsRemappedOffset(k) = length(find( trOfSet_sorted>trOfSet_sorted(k)));
                    end
                end
                
                
                
                
                %===------ plot from/to pairs (connecting lines)
                startPos = unique(transitionPairs(:,1));
                
                for pos=1:length(startPos)
                    
                    indsFromThisSet = find(transitionPairs(:,1)==startPos(pos));
                    
                    for transNr=1:length(indsFromThisSet)
                    
                        transStats = transitionPairs(indsFromThisSet(transNr,:),:);
                        
                        if transStats(1)==0 || transStats(1)==max(startPos) || transStats(2)==max(startPos)     % dont plot dummy state (0) or zero state (last one)
                            continue;
                        end
                        
                        %transStats: from to count prob
                        %transTrVals: from to
                        transTrVals = [trOfSet_sorted(transStats(1)) trOfSet_sorted(transStats(2))];
                        
                        % make line thickness proportional to probability
                        %line(transStats(1:2), transTrVals, 'color', 'k', 'LineWidth', transStats(end)*10);
                        
                        %colToUse=[1 1 1];
                        
                        transProb = transStats(end);

                        if plotProbabilities
                            colToUse = [1 1 1].*(1-transProb);
                            lineSize = transProb*lineWidthGain+lineWidthOffset;
                        else
                            colToUse=[0 0 0];
                            lineSize = 2;
                        end
                        
                        
                        connectedSetNrs = transStats(1:2);

                        lineFromPos = connectedSetNrs-setNrsRemappedOffset(connectedSetNrs);
                        lineToPos = transTrVals;
                        %if plotProbabilities

                        
                        if lineToPos(2)<lineToPos(1)
                            % with gradient
                            lineFromPos = lineFromPos - 0.05;
                            
                            %colToUse = scaleColorForLines(colToUse, 3, 0.4);   %blue                           
                            %gray is down
                            line( lineFromPos , lineToPos, 'color', colToUse, 'LineWidth', lineSize);
                        
                        else
                            lineFromPos = lineFromPos + 0.05;
                            %against gradient
                            colToUse = scaleColorForLines(colToUse, 2, 0.4);   %green                           
                            line( lineFromPos , lineToPos, 'color', colToUse, 'LineWidth', lineSize);
                            
                        end
                        %else
                        %    arrow([lineFromPos(1) lineToPos(1)]+0.1, [lineFromPos(2) lineToPos(2)]+0.1, 'Length', 10);
                        %end
                        

                    end
                end

                %===------ how many times has a state been visited (to part) (plot dots)
                endPos = unique(transitionPairs(:,2));
                nrVisitStats=[]; %setNr nrVisits
                 for pos=1:length(endPos)
                    indsFromThisSet = find(transitionPairs(:,2) == endPos(pos));
                    nrVisits = sum(transitionPairs(indsFromThisSet,3));
                    nrVisitStats = [ nrVisitStats; [endPos(pos) nrVisits]];
                 end
                 
                 %list states that were never visited but that exist
                 setsNotVisited = setdiff( 1:length(trOfSet_sorted), endPos );
                 for pos=1:length(setsNotVisited)
                     nrVisitStats = [ nrVisitStats; [setsNotVisited(pos) 0]];
                 end
                 
                totVisits = sum(nrVisitStats(:,2));
                nrVisitStats(:,3)= nrVisitStats(:,2)./totVisits;
                for pos=1:size(nrVisitStats,1)
                    setNr = nrVisitStats(pos,1);
                    visitProb = nrVisitStats(pos,3);
                    
                    if setNr==0
                        continue;
                    end
                    
                    %is this one hidden or permitted
                    trOfSet = trOfSet_sorted(setNr);
                    
                    if ~isempty(find(trVals_toPlot_forbidden==trOfSet))
                        %is forbidden
                        colDot = 'r';
                                
                    else
                        %is permitted
                        colDot = 'g';
                    end
                    
                    hold on
                    
                    if plotProbabilities
                        %if isnan(visitProb)
                        %    dotSize=markerSizeOffset+markerSizeGain*0.1;
                        %else
                            dotSize=markerSizeOffset+markerSizeGain*visitProb;
                        %end
                    else
                        dotSize=20;
                    end
                                            
                    h1=plot( setNr - setNrsRemappedOffset(setNr), trOfSet, ['o' colDot], 'MarkerSize', dotSize , 'MarkerFaceColor', colDot);
                    
                    %label the sets
                    if useStepOffsets
                       
                        text( setNr - setNrsRemappedOffset(setNr), trOfSet, [num2str(setNr)]);
                    end
                            
                end
                
                if useStepOffsets
                    
                    setList = 1:length(trOfSet_sorted);
                    allSetNrsRemapped = setList - setNrsRemappedOffset;
                    legendPos = [ max(allSetNrsRemapped) min(trOfSet_sorted) ];
                else
                    legendPos = [ 1 min(trOfSet_sorted) ];
                end
                
                
                if plotProbabilities
                    %=== make a legend for lines
                    
                    %colMapGreen = makeColorMap([0 0.4 0],[0 1 0], 100);
                    
                    for transProb=0.1:0.1:1
                        colToUse = [1 1 1].*(1-transProb);
                        
                        
                        %colToUse = scaleColorForLines(colToUse, 2, 0.4);
                        lineSize = transProb*lineWidthGain+lineWidthOffset;
                        
                        line(legendPos(1)+[0.5 1], [legendPos(2) legendPos(2)]+transProb*1.5, 'color', colToUse, 'LineWidth', transProb*10);
                        text(legendPos(1)+ 1.1, legendPos(2)+transProb*1.5, ['p=' num2str(transProb) ]);
                    end
                    
                    %=== make a legend for dots
                    colDot='r';
                    for visitProb=0:0.05:0.3
                        h1=plot( legendPos(1)+2, legendPos(2)+visitProb*6+0.1, ['o' colDot], 'MarkerSize', markerSizeOffset+markerSizeGain*visitProb , 'MarkerFaceColor', colDot);
                        
                        text( legendPos(1)+2.5, legendPos(2)+visitProb*6+0.1, ['p=' num2str(visitProb) ]);
                    end
                    
                end
                
                
                
            end
            
            
            hold off
            xlabel('set # (rank-ordered by tr)');
            ylabel('divergence');     
            
            title('blue=with gradient, green=against gradient');
    end
    
end


%==
    function colToUse = scaleColorForLines(colToUse, whichCol, minVal)
        switch(whichCol)
            case 1
                %red
                %only one color
                colToUse(2)=0;
                colToUse(3)=0;
                colTmp = scaledata([0 colToUse(1) 1], minVal, 1);
                colToUse(1) = colTmp(2);
            case 2
                %green
                %only one color
                colToUse(1)=0;
                colToUse(3)=0;
                colTmp = scaledata([0 colToUse(2) 1], minVal, 1);
                colToUse(2) = colTmp(2);
            case 3
                %blue
                %only one color
                colToUse(1)=0;
                colToUse(2)=0;
                colTmp = scaledata([0 colToUse(3) 1], minVal, 1);
                colToUse(3) = colTmp(2);
        end
    
                        