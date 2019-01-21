function W = graph_DBcells_normalizeWeights_column(W, targetWeight)

%normalize rowwise weights
for k=1:size(W,2)
    col = W(:,k);
    if sum(col)~=0
        col = col./sum(col);
        col = col.*targetWeight;
        W(:,k) = col;
        %sum( WDb_to_pyr(k,:) )
    end
end
