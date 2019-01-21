function plotPermittedHiddenSets_pairwise( setNr_map1, setNr_map2, indsPermitted_sorted, indsForbidden_sorted, totNrSets,trOfSet_sorted,eigOfSet_sorted)


pairwiseSets = zeros(totNrSets,totNrSets);

for k=1:totNrSets
    for j=1:totNrSets
%        pairwiseSets( k,j) = eigOfSet_sorted(k)+eigOfSet_sorted(j);
        pairwiseSets( k,j) = trOfSet_sorted(k)+trOfSet_sorted(j)-1;
    end
end


for k=1:length(indsPermitted_sorted)
    for j=1:length(indsPermitted_sorted)
        
        if k~=j
            pairwiseSets(indsPermitted_sorted(k),indsPermitted_sorted(j))=2;   % jointly permitted
    
        end
    end
end

%pairwiseSets( :,indsPermitted_sorted) = 1;

%%
figure(81);

subplot(2,2,1);
imagesc(pairwiseSets); colorbar;
hold on
plot( setNr_map1, setNr_map2, '-x' );
hold off

xlabel('set # (node 1) - rank-ordered by tr');
ylabel('set # (node 2) - rank-ordered by tr');

set(gca,'YDir','normal')
title('col: summed trace (if <0); >0 is permitted set combination ');

subplot(2,2,2);


%imagesc(pairwiseSets2); colorbar;
%hold on
plot( setNr_map1, setNr_map2, '-x' );
%hold off
title('trajectory in set space');
%set(gca,'YDir','normal')

    
xlim([0 totNrSets]);
ylim([0 totNrSets]);

xlabel('set # (node 1) - rank-ordered by tr');
ylabel('set # (node 2) - rank-ordered by tr');

subplot(2,2,3);

dDiv1 = setNr_map1;
dDiv2 = setNr_map2;

h=rectangle('Position',[0 indsPermitted_sorted(1)-0.5 length(dDiv1) indsPermitted_sorted(end)-indsPermitted_sorted(1)+1])
set(h,'FaceColor',[0.05 0.05 0.05]*15);

hold on
h1=plot(dDiv1, 'bx-');

h2=plot(dDiv2, 'rx-');

hold off
legend('Node 1','Node 2');
ylabel('set #');
xlabel('time');
ylim([0 10]);

title('set # as funct of time; grey=permitted');