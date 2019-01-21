%
% plot colored boxes to mark hidden/permitted sets
%
function plotBoxed_hiddenPermitted( toPlot, colOrder )


% 0 = white, 1=red, 2=green

for j=1:size(toPlot,2)
    for i=1:size(toPlot,1)
        
        switch toPlot(i,j)
            case 0
                %nothing
                col=colOrder{1};
            case 1
                col=colOrder{2};
            case 2
                col=colOrder{3};
        end
        if toPlot(i,j)>0
            rectangle( 'Position',[j-1, i-1, 1, 1], 'FaceColor',col, 'EdgeColor',col);
        end
    end
end
            
            
        
xlim([0 size(toPlot,2)]);
ylim([0 size(toPlot,1)]);
set(gca,'YDir','reverse')