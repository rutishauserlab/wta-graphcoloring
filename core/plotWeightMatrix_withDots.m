%
% plot a weight matrix using dots ; can be used as an alternative to using imagesc
%
%urut/feb16
function plotWeightMatrix_withDots(Woverall,colorMapWeightMatrix, dotSize)
if nargin<3
    dotSize=10;
end
hold on
%scale colormap

nColors = size(colorMapWeightMatrix,1);
minVal = min(Woverall(:));
maxVal = max(Woverall(:));
range=maxVal-minVal;


for k=1:size(Woverall,1)
    for j=1:size(Woverall,2)
        
        Wval = Woverall(k,j);
        
        if Wval~=0
           
            colInd = ((Wval/range)-minVal)*nColors;
            
            %scale color into colormap
            colInd = round(rescale_range(Wval, minVal, maxVal,1, nColors));
            if colInd==0
                colInd=1;
            end
            if colInd>nColors
                colInd=nColors;
            end
            
            plot( j, k, '.', 'color', colorMapWeightMatrix(round(colInd),:), 'MarkerSize', dotSize);
        end
        
    end
end
hold off;
set(gca,'YDir','reverse');

xlabel('unit nr [from]');
ylabel('unit nr [to]');

xlim([0 size(Woverall,2)+1]);
ylim([0 size(Woverall,1)+1]);