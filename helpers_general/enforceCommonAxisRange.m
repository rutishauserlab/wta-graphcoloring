%
%axisList: empty (take all subplots of current figure) or list of axis
%whichAxis: 1=x, 2=y, 3=both
%
%urut/sept10/MPI
function enforceCommonAxisRange( axisList, whichAxis, presetX,presetY )
if nargin<3
    presetX=[];
end
if nargin<4
    presetY=[];
end

if isempty(axisList)
    axisList = findall(gcf,'Type','axes');
end

%% find the min/max of each subplot
Xlimits = [0 0];
Ylimits = [0 0];
for k=1:length(axisList)
   if strcmp(get(axisList(k),'Tag'),'legend')
       %dont modify if it is a legend
       continue;
   end
    xlimsNew = xlim( axisList(k) );
    ylimsNew = ylim( axisList(k) );
    
    [minValTmp,maxValTmp] = updateMinMaxVals(Xlimits(1),Xlimits(2), xlimsNew);
    Xlimits =  [minValTmp maxValTmp];
   
    [minValTmp,maxValTmp] = updateMinMaxVals(Ylimits(1),Ylimits(2), ylimsNew);
    Ylimits =  [minValTmp maxValTmp];
    
end

%% set common
if nargin>2
    Xlimits = presetX;
    Ylimits = presetY;
end

for k=1:length(axisList)
   if strcmp(get(axisList(k),'Tag'),'legend')
       %dont modify if it is a legend
       continue;
   end

    if whichAxis==1 || whichAxis==3
       xlim( axisList(k), Xlimits ); 
    end
    if whichAxis>=2
       ylim( axisList(k), Ylimits ); 
    end
end