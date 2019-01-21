%
%count nr external biases that are set but not satisfied by the network
%
function [errBias,whichHaveError] = sudoku_nrBiasErrors(Mhistory,Nnodes, biasInp)
[maxAmp,winnerID,ampInhib,ampInput] = getAmpOfWinner(Mhistory, Nnodes);

errBias=0;
whichHaveError=[];
for k=1:Nnodes
    
    winnerOfNode = winnerID(k);
    
    if biasInp(k)>-1
        
        if biasInp(k) ~= winnerOfNode
            errBias=errBias+1;
            whichHaveError = [whichHaveError k];
        end        
    end
end