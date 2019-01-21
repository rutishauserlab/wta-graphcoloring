%
%dendritic non-linearity
%
%  params = [ s o ]   where s is slope, o is offset
%
function Y=fShunt(X, params)

%Y = 0.5.*(1-tanh((X./0.5)-8));  % DEFAULT, works for sudoku

%testing
offset=4; %where is maximal steepness
Y = 1 - (1/2 * ( tanh(params(1)*(X-params(2)))+1));    % sigmoid function is 1/2*(tanh(slope*(x-offset))+1)

%Y = 1-tanh((X.*slope));  % testing

%Y = 1-min(max(slope.*X,0),1);  % used for analyical arguments, also works for sudoku (a bit worse)

debugPlots=0;
if debugPlots
    %
    % %===== debug figures
    % slope=0.5;
    X=0:0.01:20;
    
    
    %Y2 = tanh((X.*slope) );
    
    % Y = 1-min(max(slope.*X,0),1);
    %Y=0.5.*(1-tanh((X./0.5)-8));
    figure(21);
    subplot(2,2,1);
    
    slopes = [ 0.1 0.2 0.3 0.5  1.1 2 4 ];
    hs=[];
    legendStrs=[];
    for k=1:length(slopes)
        slope=slopes(k);
        
        if slope~=1
        Y = (1-tanh((X.*slope) ));  % testing
        offset=0;
        %Y = 1 - (1/2 * ( tanh(slope*(X-offset))+1));    % sigmoid function is 1/2*(tanh(slope*(x-offset))+1)

        legendStrs{k} = ['s=' num2str(slope)];
        else
            
            offset=4; %where is maximal steepness
            slopeD=2;
            Y = 1 - (1/2 * ( tanh(slopeD*(X-offset))+1));
            %Y = 0.5*(1-tanh(( 2* X./0.5)-8));
            legendStrs{k} = ['s=' num2str(slope) ' -8'];
        end
        
        if k>1
           hold on;
        end
        hs(k)=plot(Y,X, 'color', rotatingColorCode(k));
        
    end
    hold off;
    legend(hs,legendStrs);
    ylabel('total inhibition');
    xlabel('dendritic sensitivity g(z)');
    %subplot(2,2,3);
    %plot(X,Y2);
    
end

 
