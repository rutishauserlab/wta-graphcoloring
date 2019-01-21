%
%
% random network
%

%% setup
alpha1=1.2;
alpha2=0.5;

beta1=3;     % 
beta2=0.25;

beta1_DB=3;  % at least (3+beta2)/2
beta2_DB=0.25;  % max 0.5*(-4+alpha1+beta1+2*beta1_DB-beta2)

params=[];
params.dt = 0.01;
params.inpOnTime=[2000 16000; 0 0];
params.s=1;
params.maxCutoff = 0;
params.setConditionalBranchingEnabled=0;

T=0;
G_E=1;
G_I=1.5;

enumerateSets = 0;

%% design random network
gridSize = 10;

Poccupied = 0.4;
Pexcit = 0.8;  % Pinhib is 1-Pexcit
Pinhib2 = repmat(0.2, 1, 8);  % P that excitatory unit connects to a second / third inhib unit [Psecond Pthird ...]
PexcitCross = 0;

%generate a new network. keep as-is to re-use existing network
unitGrid=[];
unitList=[];
connectionList1=[];

connectionList2=[];


beta1_DB=1;  % at least (3+beta2)/2


[unitGrid,unitList,connectionList1, connectionList2,connectionListExcit, W,indsInhib,indsExcit] = build_rand_network(gridSize,Poccupied,Pexcit, Pinhib2, PexcitCross, alpha1,alpha2, beta1,beta2, beta1_DB, beta2_DB, unitGrid, unitList, connectionList1, connectionList2);
N=size(W,1);

figure(2);
subplot(3,3,1);
imagesc(unitGrid);
colorbar
title('location of units (red=inhib, green=excit, blue=none)');
xlabel('x position');
ylabel('y position');

subplot(3,3,2);
imagesc(W);
colorbar;
title('weight matrix');
xlabel('unit nr [from]');
ylabel('unit nr [to]');

%disp(['nrTot=' num2str(N) ' nrExcit=' num2str(length(indsExcit)) ' %' num2str(100*length(indsExcit)/N) ' nrInhibConn=' num2str(size(connectionList,1)) ' nrExcitConn=' num2str(size(connectionListExcit,1))])

paramsUsed.Poccupied=Poccupied;
paramsUsed.Pexcit=Pexcit;
paramsUsed.Pinhib2=Pinhib2;
paramsUsed.PexcitCross=PexcitCross;
paramsUsed.gridSize=gridSize;
paramsUsed.unitGrid=unitGrid;
paramsUsed.unitList=unitList;

%% 

Gnetwork = ones(N,1) * G_E;
Gnetwork(indsInhib) = G_I;

Tnetwork = ones(N,1) * T;
Tnetwork(indsInhib) = 0;    % none for inhib

inpDefault=zeros(N,1);  % no input

I=6.01; %external input 

inputAmp = randn(N,1)/5+I;
inputAmp(indsInhib)=0;   %no input to inhib

loadTerms=diag(-Gnetwork);

nrSteps=15000;

noiseEnabled=0;

%% define set numbers for this network
if enumerateSets
[allSets, indsPermitted, indsForbidden, eigOfSet, trOfSet, allSets_sorted, indsPermitted_sorted,indsForbidden_sorted,...
    eigOfSet_sorted, trOfSet_sorted] = randNet_define_permittedHiddenSets( W, indsInhib, indsExcit, loadTerms );
end

%% run

%run simulation
[x,inpHistory, DD] = collective_runSim(N, W, Gnetwork, Tnetwork, nrSteps, params, inputAmp, [], inpDefault, noiseEnabled );

%match this run against the set numbers
if enumerateSets
    DD_matchedSet = convert_DD_toSetNr(DD,allSets_sorted, indsExcit);
else
    DD_all_tmp=[];
    DD_all_tmp(1).DDstats  = DD;
    
    [allSets_generated, indsPermitted_sorted, indsForbidden_sorted, eigOfSet, trOfSet] = enumerateSets_dynamic( DD_all_tmp, 1, indsExcit, W, loadTerms );
    
    %generate setsVisited based on dynamic set numbers
    DD_matchedSet = convert_DD_toSetNr(DD,allSets_generated, indsExcit);
end


%% plot a single run
figure(2);
subplot(3,3,3);
imagesc(x);
colorbar
title('example run');
xlabel('time');
ylabel('unit nr');

subplot(3,3,4);
%if enumerateSets
    plot(DD_matchedSet);
%end

winnerSetNr = DD_matchedSet(end);

%when was the winner reached
indsReached = find( DD_matchedSet == winnerSetNr);
timeWinnerReached = indsReached(1)

if ~isempty( find( indsPermitted_sorted==winnerSetNr) )
    winnerIsPermitted=1;
else
    winnerIsPermitted=0;
end
title(['Example run; winnerID=' num2str(winnerSetNr) ' winnerReached=' num2str(timeWinnerReached) ]); %' isPermitted=' num2str(winnerIsPermitted)]);

%subplot(3,3,5);
%imagesc( allSets_sorted );

%allSets_sorted_plot = allSets_sorted;            
%allSets_sorted_plot(indsPermitted_sorted,:) = allSets_sorted_plot(indsPermitted_sorted,:)*2;
            
%plotBoxed_hiddenPermitted( allSets_sorted_plot, {'w','r','g'} )


%% many runs
tic

evalWinnerTime = 1;  % turn on to save in timeWinnerReachedAll the time the winning set was first entered
timeWinnerReachedAll=[];

nrRuns = 100;
setsVisited=[];
DD_all=[]; %zeros(runNr,nrSteps);
parfor runNr=1:nrRuns
    %generate new input
    
    inputAmp = randn(N,1)+I;
    inputAmp(indsInhib)=0;   %no input to inhib

    [x,inpHistory, DD] = collective_runSim(N, W, Gnetwork, Tnetwork, nrSteps, params, inputAmp, [], inpDefault, noiseEnabled );

    %match this run against the set numbers

    if enumerateSets
        DD_matchedSet = convert_DD_toSetNr(DD,allSets_sorted, indsExcit);
        setsVisited(runNr,:) = DD_matchedSet;
    end

    if evalWinnerTime
  
        DD_all_tmp=[];
        DD_all_tmp(1).DDstats  = DD;
          [allSets_generated, indsPermitted_sorted, indsForbidden_sorted, eigOfSet, trOfSet] = enumerateSets_dynamic( DD_all_tmp, 1, indsExcit, W, loadTerms );
        DD_matchedSet = convert_DD_toSetNr(DD,allSets_generated, indsExcit);

        winnerSetNr = DD_matchedSet(end);

        %when was the winner reached
        indsReached = find( DD_matchedSet == winnerSetNr);
        timeWinnerReachedAll(runNr) = indsReached(1);
    
    end
    
    
    DD_all(runNr).DDstats  = DD;

    disp(['run nr ' num2str(runNr) ' completed']); 
end
toc

%% winner time comparison
if ~exist('totStatsWinnertime')
totStatsWinnertime = [];
end
totStatsWinnertime=[totStatsWinnertime; [gridSize Pinhib2(1) mean(timeWinnerReachedAll) std(timeWinnerReachedAll) length(timeWinnerReachedAll)]];
timeWinnerReached_history{size(totStatsWinnertime,1)}=timeWinnerReachedAll;

%%
figure(40);

bar( 1:size(totStatsWinnertime,1), totStatsWinnertime(:,3)); %, 'bardwith',0.5);

hold on
errorbar( 1:size(totStatsWinnertime,1), totStatsWinnertime(:,3), totStatsWinnertime(:,4), '.');
hold off

set(gca,'XTick', 1:size(totStatsWinnertime,1));
set(gca,'XTickLabel', totStatsWinnertime(:,1));
ylim([0 10000]);

xlabel('Grid Size');
ylabel('time permitted set entered');

ttest2( timeWinnerReached_history{2}, timeWinnerReached_history{3} )
[h,p]=ttest2( timeWinnerReached_history{3}, timeWinnerReached_history{4} )


%% dynamically generate visited sets if too many to enumerate

if ~enumerateSets

    [allSets_generated, indsPermitted_sorted, indsForbidden_sorted, eigOfSet, trOfSet] = enumerateSets_dynamic( DD_all, nrRuns, indsExcit, W, loadTerms);
    
    %generate setsVisited based on dynamic set numbers
    parfor k=1:nrRuns
        setsVisited(k,:) = convert_DD_toSetNr(DD_all(k).DDstats,allSets_generated, indsExcit);
    end
        
end

winnerSetsAll = setsVisited(:,end);   % which permitted set was reached at the end

nrWinners = length(unique(winnerSetsAll));

nrWinners
size(allSets_generated)

%%
[H, nrSetsVisited] = calcEntropy_ofWTAstate(setsVisited);

[avTr, avSetNr] = calcSetStats_ofWTAstate( setsVisited, trOfSet );

xlimCommon=[0 7000];

figure(16);


t = 1:length(H);

subplot(2,2,1);

[ax,h1,h2] = plotyy(t, H, t, avTr);

set(h1(1),'color','k', 'linewidth', 1);
set(h2(1),'color','r', 'linewidth', 1);
set(ax(1),'YColor','k');
set(ax(2),'YColor','r');

ylabel(ax(1),'Entropy [bits]');
ylabel(ax(2),'Trace');

set(ax(2),'YTick', min(avTr)-2:5:0);
%ylabel('Entropy [bits]');

xlim(ax(1),xlimCommon);
xlim(ax(2),xlimCommon);



%plot(H);
ylabel('entropy [bits]');

subplot(2,2,2);
plot(nrSetsVisited);
ylabel('nr different sets visited')

subplot(2,2,3);
plot( avTr);
ylabel('trace');


%%
save('run3_constDens2_R300.mat', 'H','nrSetsVisited','avTr','avSetNr','W','paramsUsed');

%% summarize runs
runToLoad=[1:3];

figure(700);
for k=1:length(runToLoad)
    runNr = runToLoad(k);

    fName = ['run' num2str(runNr) '_constDens2_R300.mat'];
    
    load(fName); 
    if k>1
        hold on;
    end
    
    plot( t, H, rotatingColorCode(k));
    legendStrs{k}=[fName];
end
hold off
legend(legendStrs);

xlim([0 15000]);
xlabel('time [ms]');
ylabel('Entropy [bits]');

