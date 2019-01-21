%
% all possible pairs of two vectors
% grp1 and grp2 need to be unique (each item)
%
%urut/sept11
function pairs = getAllPossiblePairs( grp1,grp2 )
pairs=[];
c=0;


for k=1:length(grp1)
    for j=2:length(grp2)
        from=grp1(k);
        to=grp2(j);
        
        if from~=to
            if isempty(pairs) | (isempty( find(pairs(:,1)==from & pairs(:,2)==to)) & isempty( find(pairs(:,2)==from & pairs(:,1)==to)) )
                c=c+1;
                pairs(c,:) = [ from to ];
            end
        end
    end
end