
%
% main of collective computation simulations
%
% this file simulates the simple network (1-WTA with N units total)
%
%urut/june2014


%% setup

%% parameters
alpha1=1.2;
alpha2=0;

beta1=3;     % 
beta2=0.25;

beta1_DB=1.7;  % at least (3+beta2)/2
beta2_DB=0.3;  % max 0.5*(-4+alpha1+beta1+2*beta1_DB-beta2)

params=[];
params.dt = 0.01;
params.inpOnTime=[2000 14000; 0 0];
params.s=1;
params.maxCutoff = 0;

T=0;
%G_E=0.9;
G_E=1.1;
%G_I=1.0;
G_I=1.5;
%G_I=2.0;

%G_E=1;
%G_I=1.2;

gain = calcGain(alpha1, beta1, beta2, G_E, G_I)

% setup network
N=5;  % total nr units in the networkW=zeros(N,N);

W=zeros(N,N);

for k=1:N-1
    W(k,k)=alpha1;
end

W(N,1:end-1) = beta2;
W(1:end-1,N) = -beta1;

%W(N-1,1:end-2) = beta2;
%W(1:end-2,N-1) = -beta1;
%W(N-1,N-1)=0;
%W(N-1,N)=0;
%W(N,N-1)=0;

%W(2,1)=0.2;  %cond1
%W(2,4)=0.2;  %cond2

%W(2,4)=-0.2;  %cond2


Gnetwork = ones(N,1) * G_E;
Gnetwork(end) = G_I;

Tnetwork = ones(N,1) * T;
Tnetwork(end) = 0;    % none for inhib

inpDefault=zeros(N,1);  % no input

I=2.01; %external input 

inputAmp = randn(N,1)/5+I;
inputAmp(end)=0;   %no input to inhib

inputAmp2 = randn(N,1)+I;
inputAmp2(end)=0;   %no input to inhib

%inputAmp(end-1)=0;
loadTerms=diag(-Gnetwork);

nrSteps=13500;

noiseEnabled = 0;

params.setConditionalBranchingEnabled = 0;
params.setConditional = [0 0 1 1 1]';  % which set is the trigger
params.setConditionalInpNr =[1 6]; %which input is changed how much

%% run

%additional wiring
W(2,1)=0.2;  %cond1
W(2,4)=0.2;  %cond2

%no additional
%W(2,1)=0;  %cond1
%W(2,4)=0;  %cond2

%run simulation
[x,inpHistory, DD] = collective_runSim(N, W, Gnetwork, Tnetwork, nrSteps, params, inputAmp, inputAmp2, inpDefault, noiseEnabled );
% determine sets visited in this run
[DD_matchedSet1, totNrSets, allSets_sorted, trOfSet_sorted, indsPermitted_sorted, indsForbidden_sorted] = plotPermittedHiddenSets_prepare(  W, loadTerms, DD );


winnerTime=10000;
winnerID = find ( x(:,winnerTime) > 0.01 )


%% plot single sim
figNr=11;


[maxEig,minEig,trEig,condNr4] = calcContractStats_forRun(nrSteps, DD, N, W, loadTerms, winnerID );

figure(figNr);

unitsToPlot=1:5;
xlimCommon=[0 7000];
collective_plotOneRun( x, inpHistory, nrSteps, unitsToPlot, maxEig, minEig, trEig,xlimCommon );




%--- permitted/forbidden sets for this run
plotMode=1;
figure(81);
[allSets_sorted,trOfSet_sorted, firstWinnerInd] = plotPermittedHiddenSets( W, loadTerms, DD, plotMode, [] );


%===== plot set transitions for single run

[transitionPairs,dAllTr] = prepareSetSwitchStats( DD_matchedSet1', trOfSet_sorted );

figure(94+figNr);
%subplot(2,2,1);

plotMode=3;
[allSets_sorted,trOfSet_sorted, firstWinnerInd] = plotPermittedHiddenSets( W, loadTerms, DD, plotMode, transitionPairs );

ylim([-5.5 -0.5]);
set(gca,'YTick',-5:1-1);
xlim([0 7]);

%% many runs of the same network
tic
nrRuns = 500;

setsVisited=[];
DD_all=[]; %zeros(runNr,nrSteps);
parfor runNr=1:nrRuns
    %generate new input
    inputAmp = randn(N,1)+I;
    inputAmp(end)=0;   %no input to inhib
    
    %run simulation
    [x,inpHistory, DD] = collective_runSim(N, W, Gnetwork, Tnetwork, nrSteps, params, inputAmp, inputAmp2, inpDefault, noiseEnabled );
    
    DD_all(runNr).DDstats  = DD;
    %DD_all{runNr}=DD;
    
    % determine sets visited in this run
    [DD_matchedSet1, totNrSets, allSets_sorted, trOfSet_sorted, indsPermitted_sorted, indsForbidden_sorted] = plotPermittedHiddenSets_prepare(  W, loadTerms, DD );
    
%        DD_matchedSet = convert_DD_toSetNr(DD,allSets_sorted, indsExcit);
%        setsVisited(runNr,:) = DD_matchedSet;
        
    %store the set numbers that were visited
    setsVisited(runNr,:) = DD_matchedSet1;
   
    disp(['run nr ' num2str(runNr) ' completed']); 
end
toc

%% summarize transition history (transition probability) for a block of simulations

setsVisited_endPointConditional = [];

%endPointToUse=-1;  
%selFor=71;

endPointToUse=100;  
selFor=13;


toUseInds = selectRuns_byEndpoint_reached(setsVisited, endPointToUse, selFor);
setsVisited_endPointConditional = setsVisited(toUseInds,:);


[transitionPairs,dAllTr] = prepareSetSwitchStats( setsVisited, trOfSet_sorted );

%transitionPairs = prepareSetSwitchStats( setsVisited_endPointConditional, trOfSet_sorted );

figure(94);
%subplot(2,2,1);

plotMode=3;
[allSets_sorted,trOfSet_sorted, firstWinnerInd] = plotPermittedHiddenSets( W, loadTerms, DD, plotMode, transitionPairs );

xlim([0 9])

figure(95);
hist(dAllTr);
title('hist of trace jumps at transitions');
xlabel('difference divergence');
ylabel('nr transitions');

nrPos = length(find(dAllTr>0))
percPos = nrPos/length(dAllTr)

%% plot set number histograms, endpoint dependent

endPointToUse=100;
subplotX=3;
subplotY=3;

setNrsToUse=[12];
plotPos=[2 3];

figure(11);
plotSetNrHistory_manyRuns(setsVisited, setNrsToUse, endPointToUse, subplotX,subplotY, plotPos, allSets_sorted, indsPermitted_sorted, indsForbidden_sorted)

subplot(3,3,plotPos(1));
xlim([1000 10000]);

setNrsToUse=[13];
plotPos=[5 6];
plotSetNrHistory_manyRuns(setsVisited, setNrsToUse, endPointToUse, subplotX,subplotY, plotPos, allSets_sorted, indsPermitted_sorted, indsForbidden_sorted)

subplot(3,3,plotPos(1));
xlim([1000 10000]);


%% entropy calculation from a set of runs with different inputs
H = calcEntropy_ofWTAstate(setsVisited);

if ~exist('Hall')
    Hall=[];
end

%Hall(size(Hall,1)+1,:) = [alpha1 H];

Hall(size(Hall,1)+1,:) = [G_E H];

%% ---- compare entropy with other metrics of the network
[avTr, avSetNr] = calcSetStats_ofWTAstate( setsVisited, trOfSet_sorted );

avEigVals = calcEigStats_ofWTAstate(DD_all(1:30), N, W, loadTerms, nrSteps );


inpTime=zeros(1,nrSteps);

xlimCommon=[0 7000];

t=1:length(inpTime);
figure(31);
subplot(3,3,1)
%plot( H );

[ax,h1,h2] = plotyy(t, H, t, avTr);

set(h1(1),'color','k', 'linewidth', 1);
set(h2(1),'color','r', 'linewidth', 1);


ylabel(ax(1),'Entropy [bits]');
ylabel(ax(2),'Trace');

set(ax(2),'YTick', min(avTr)-2:1:0);
%ylabel('Entropy [bits]');

xlim(ax(1),xlimCommon);
xlim(ax(2),xlimCommon);

subplot(3,3,2);

[ax,h1,h2] = plotyy(t, H, [t; t]', avEigVals(1:2,:)' );

ylabel(ax(1),'Entropy [bits]');
ylabel(ax(2),'Eigenvalue');

set(h1(1),'color','k', 'linewidth', 1);
set(h2(1),'color','r', 'linewidth', 1);
set(h2(2),'color','b', 'linewidth', 1);


xlim(ax(1),xlimCommon);
xlim(ax(2),xlimCommon);


subplot(3,3,4);
[ax,h1,h2]=plotyy(t, H, t, avSetNr);
ylabel(ax(1),'average trace');
ylabel(ax(2),'average set nr');

set(h1(1),'color','k', 'linewidth', 1);
set(h2(1),'color','r', 'linewidth', 1);


xlim(ax(1),xlimCommon);
xlim(ax(2),xlimCommon);

subplot(3,3,7);
inpTime=zeros(1,nrSteps);
inpTime(params.inpOnTime(1,1):params.inpOnTime(1,2)) = 1;
plot( t, inpTime );
ylim([0 2]);
ylabel('input on/off');

xlim(xlimCommon);

%% compare different entropies

alphaVals = Hall(:,1);
HVals = Hall(:,2:end);

%strNotes={'noise on 1','noise off','noise on 2'};
%strNotes={'+0 Constr','+1 Constr','+2 Constr'};
%strNotes={'','','',''};
%strNotes={'G_I=1.0 gain=1.8','G_I=1.5, gain=3.3','G_I=2.0, gain=5.7'};
strNotes={'G_E=0.9, gain=5','G_E=1.1, gain=2.5'};

inpTime=zeros(1,nrSteps);
inpTime(params.inpOnTime(1,1):nrSteps) = 1;


t=1:1:size(HVals,2);

figure(101);
strs=[];
hs=[];
for k=1:length(alphaVals)

    if k>1
        hold on;
    end
    
    if k<length(alphaVals)
        hs(k) = plot( t,(HVals(k,:)), rotatingColorCode(k) );
    else
        
        [hAxes,hLine1,hLine2] = plotyy( t, (HVals(k,:)), t, inpTime); %, rotatingColorCode(k));
        hs(k)=hLine1;
        set(hLine1,'Color', rotatingColorCode(k));
        set(hLine2,'Color', 'k');
        set(hLine2,'LineWidth', 2);
        
        ylim(hAxes(2), [0 10]);
    end
        
    if length(strNotes)>=k
        %strs{k}=['\alpha=' num2str(alphaVals(k)) ' ' strNotes{k}];
        strs{k}=[strNotes{k}];
    else
        strs{k}='';
    end
    
end



hold off
legend(hs,strs);

ylabel(hAxes(1), 'Entropy [bits]');
ylabel(hAxes(2), 'Input On/Off');
xlabel('time');

xlim([1000 10000]);
xlim(hAxes(1),[1000 10000]);
xlim(hAxes(2),[1000 10000]);

