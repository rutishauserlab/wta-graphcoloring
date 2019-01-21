function graphPlot_sudokuResult(tillPlot, NrOfColors, Nnodes,Mhistory, paramStr, Maps, errorNodes, biasInp)


nodeMatrix=[];
for rowNr=1:NrOfColors  
    rowNodes=[1:NrOfColors] + (rowNr-1)*NrOfColors;    
    nodeMatrix(rowNr,:) = rowNodes;
end
nodeMatrixSolution=zeros(NrOfColors,NrOfColors);

[maxAmp,winnerID,ampInhib,ampInput] = getAmpOfWinner(Mhistory, Nnodes);

errBias=0;

subplot(2,2,[1 3]);
for k=1:Nnodes
    if k>1
        hold on
    end
   winnerOfNode = winnerID(k);
    
   if find(errorNodes==k)
       errStr='EE';
   else
       errStr='';
   end
   
   %errStr= [ errStr ' N' num2str(k)];  %to debug whether plotting is correct
   
   hasBias=0;
   biasStr='';
   if ~isempty(biasInp)
      if biasInp(k)>-1
         %is preset

         hasBias=1;
         
         %biasStr = [' S' num2str(biasInp(k))];
         
         if biasInp(k) ~= winnerOfNode
             errBias=errBias+1;
         end
      end
      
   end
   
   [rind,cind] = find(nodeMatrix==k);
   
   if hasBias
    text( cind, rind, [num2str(winnerOfNode) errStr ], 'FontWeight','bold','FontSize',24,'color','r'); 
   else
    text( cind, rind, [num2str(winnerOfNode) errStr ], 'FontWeight','bold','FontSize',24,'color','k');        
   end
   
   %if ~isempty(biasStr)
   %    text( cind+0.2, rind, [biasStr], 'FontWeight','bold','FontSize',12,'color','r');        
   %end
   
   nodeMatrixSolution(rind,cind) = winnerOfNode;
end


hold off
set(gca,'YDir','reverse');
xlim([0 NrOfColors+1]);
ylim([0 NrOfColors+1]);

set(gca,'XTick',1:9);
set(gca,'YTick',1:9);
ylabel('row number');
xlabel('column number');

subplot(2,2,2);
plot( 1:NrOfColors, sum(nodeMatrixSolution),'r', 1:NrOfColors, sum(nodeMatrixSolution'), 'b' );
legend('row sum','column sum');
xlabel('position');
ylabel('sum');
ylim([0 50]);


% check correctness of the solution
e=checkSudokuCorrectness(nodeMatrixSolution);
title(['hasErrors=' num2str(e) ' biasErrors=' num2str(errBias) ]);

if e==0
    display(['Sudoku solution is correct']);
else
    warning('Sudoku solution has errors');
end

