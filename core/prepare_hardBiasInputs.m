function [Maps,biasInps] = prepare_hardBiasInputs(Maps,Nnodes,graphType,biasInp_forSudoku,biasAmp,determineRandBiases)

%set the bias inps
for k=1:Nnodes
    biasInps = zeros(Maps(k).N,1);
    WofMap = Maps(k).W;
    if graphType==2
        
        
        %experiment with random bias to all
    
%        biasInps(1:Maps(k).N-1) = randn(1,Maps(k).N-1)/100+0.2;
        
% 
%         %only one fixed bias
%         if k==1 
% %           biasInps(1) = biasAmp+randn/100;  %so isnt exactly the same
%         end
%         if k==6
% %           biasInps(2) = biasAmp+randn/100;  %so isnt exactly the same
%         end
%         
% %         %if k==2
% %         if biasInp_forSudoku(k)>-1
% %             %biasInps(1) =biasAmp;    %want this color on this node
% %             
% %             %biasUnit=3;
% %             
% %             biasUnit=biasInp_forSudoku(k);
% %             
% %             for jj=1:NrOfColors   %loop through all to make sure old biases are removed
% %                 if jj==biasUnit
% %                     WofMap(jj,jj)=alpha1+0.2;
% %                     disp(['biasing map ' num2str(k) ' unit ' num2str(jj)]);
% %                 else
% %                     WofMap(jj,jj)=alpha1;
% %                 end
% %             end
% %         end
    end
    
    if graphType==3 || (graphType==2 & determineRandBiases==1)
        if exist('biasInp_forSudoku');
            if biasInp_forSudoku(k)>-1
                jj=biasInp_forSudoku(k);
                %WofMap(jj,jj)=alpha1+0.10;
                
                valToSet =  biasAmp+randn/100;   % default, low variance
                
                %valToSet =  biasAmp+ (rand-0.5)*15;   % testing, large variance
                
                % testing, inp is one of two values (randomly)
                %if rand<0.5
                %    valToSet = biasAmp/3+randn/100;
                %else
                %    valToSet = biasAmp+randn/100;
                %end
                                
                biasInps(biasInp_forSudoku(k)) = valToSet;  %so isnt exactly the same

                disp(['GraphType=3 biasing map ' num2str(k) ' unit ' num2str(jj) ' val=' num2str(valToSet) ]);
            end
        end
    end
    
    if graphType==4
        biasInps(1)=biasAmp;
        
        %WofMap(2,2)=1.5-0.3;   % lower gain of the other
        
    end
    
    Maps(k).biasInps = biasInps;
    Maps(k).W=WofMap;
end
