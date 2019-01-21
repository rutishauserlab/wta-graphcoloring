function W = graph_DBcells_normalizeWeights_row(W, targetWeight)

%normalize rowwise weights
for k=1:size(W,1)
    row = W(k,:);
    if sum(row)~=0
        row = row./sum(row);
        row = row.*targetWeight;
        W(k,:) = row;
        %sum( WDb_to_pyr(k,:) )
    end
end
