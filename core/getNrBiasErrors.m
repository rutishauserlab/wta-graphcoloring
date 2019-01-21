function errBias = getNrBiasErrors( Mhistory,Nnodes, biasInp,t)

if isempty(t)
     t= size( Mhistory(1).S, 2 );  % last timepoint   
end

[maxAmp,winnerID,ampInhib,ampInput] = getAmpOfWinner(Mhistory, Nnodes,t);

errBias=0;
for k=1:Nnodes
    if biasInp(k)>-1
        %is preset
        
        winnerOfNode = winnerID(k);
        
        
        if biasInp(k) ~= winnerOfNode
            errBias=errBias+1;
        end
    end
end