function W = graph_addProportionalWeightNoise( W, noiseStdProp )

indsToAdd = find( W  ~= 0);
for k=1:length(indsToAdd)
    W( indsToAdd(k) ) = W( indsToAdd(k) )+randn*W( indsToAdd(k) )*noiseStdProp;
end
