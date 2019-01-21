%simulate an arbitrary number of connected WTAs
%
%M and Conn are structures defining the network (see below for details)
%M contains the definition of each map
%Conn the connectivity between the M's
%
%enableInteract: yes/no to enable connectivity between the maps.
%enableNoiseInp: yes/no
%
%TOff is an offset for each unit: activation funct is f(x) = max(x+TOff-T,TOff)
%
%enableDendriticCompartment: use the dendritic segment yes or no. 1 is yes (threshold=0), 0 is no, threshold=-999
%economicMode: 0 no, 1 yes (dont keep history nor past state)
%
%debugMode: 0 no, 1 yes. returns a history of the state (slow)
%
%dendParams: [s o], s is slope, o is offset. s: >0 slope of dendritic non-linearity g(z). if ==-1, disabled
%
%
%this is an optimized version of runN_WTA_withDBs.m
%
%this version forked from runN_WTA_withDBs_optimized.m --- this version has the dendritic non-linearity
%
%
%initial: urut/march10
%
function [M,DB_state,i, paramsUsed,nrErrorsOverTime, nrSwitches, stateHistoryAll, shuntHistoryAll] = runN_WTA_withDBs_optimized_shunt( M, Mold, DB_state_old, Conn, dt, till, enableInteract, enableNoiseInp,noiseUpdateFreq,...
    noiseStd, noiseMean, extInpClamp, updateFunct, updateFreq, hGraph, TOff, WPyr_to_db, WDb_to_pyr,T_DB,terminateIfNoError,biasInp, keepHistory,...
    graphType, addDebugInfo, enableDendriticCompartment, economicMode, debugMode, dendParams, enableInitState, initStateVals )

N_DBs = size( WPyr_to_db,1);

%counterCorrectRequired=5;  % DEFAULT (5 rounds correct before terminate)
counterCorrectRequired=2;  % to save time for large simulations

%addDebugInfo = 1; % if enabled, stores effective weight matrices (slow)
stateHistoryAll=[];
shuntHistoryAll=[];
paramsUsed=[];

%allocate space for history on each map
for k=1:length(M)
    if ~economicMode
        M(k).S = zeros( M(k).N,till);
    else
        M(k).S = zeros( M(k).N,1);
    end
    if keepHistory
        M(k).inpHistory = zeros( M(k).N,till);
        M(k).DBInpHistory = zeros( length(M(k).overallInd), till);
        M(k).effWeightMatrix=[];
    end
end

%allocate space for history of the DBs
if ~economicMode
    DB_state = zeros( N_DBs, till);
else
    DB_state = zeros( N_DBs, 1);    
end

DB_stateMult = repmat( DB_state(:,1)', M(1).N, 1); 

%if previous sim is resumed, re-initialize by copying the state
if ~isempty(Mold)
    disp('Initialize state from old sim');
    for k=1:length(M)
        M(k).S(:,1) = Mold(k).S(:,end);
        if keepHistory
            M(k).inpHistory(:,1) = Mold(k).inpHistory(:,end);
        end
    end
end
if ~isempty(DB_state_old)
      disp('Initialize DB state from old sim');
  
      DB_state(:,1) = DB_state_old(:,end);
end


outFlag=0;
corrMechanismMsgGiven=0;

noiseMeanFrozen=[];
for k=1:length(M)
    noiseMeanFrozen(:,k) = abs(randn( M(k).N-1 ,1)*noiseStd+noiseMean) ;   %make sure it is never negative
end

TOff_perMapOrig = zeros(length(TOff), length(M) );
for k=1:length(M)
    TOff_perMapOrig(:,k) = TOff;
end
TOff_current =  TOff_perMapOrig;

%maxCrossInhib = -0.35*7 - 1 ;  % gamma*steadyAmp, plus some safety margin for variability
maxCrossInhib = -9999 ;  % unlimited
%updateFreqPerMap = round(rand(1,length(M))*5+10); %randomized
updateFreqPerMap = 10*ones(1,length(M)); %standard

inputOfMap = zeros(M(k).N,length(M));  % cache the external input in case noise is not updated this step

%gNonLinOffset_val = -100; % g(x)=x if negative and high; 0 makes g(x)=max(X,0)

if enableDendriticCompartment
    gNonLinOffset_val=0;
else
    gNonLinOffset_val=-999;
end



gNonLinOffset = repmat(gNonLinOffset_val, M(1).N,1);

paramsUsed.maxCrossInhib=maxCrossInhib;
paramsUsed.updateFreqPerMap=updateFreqPerMap;

nrUnits_perMap = M(1).N;
nrUnits = nrUnits_perMap*length(M);

if enableNoiseInp == 3
    savedPrecomputedNoise = 0;  % 0 none, 1 generate and save it, 2 load from file and use it   (for testing with same input, no new noise is generated)    
    shuffled1FNoise = 0; % 0 none, 1 yes; if yes, there are no autocorrelations in the noise
    if savedPrecomputedNoise==2
        warning('re-useing pre-computed noise,no new variability');
    end
    if shuffled1FNoise
        warning('shuffeld 1/f noise is enabled');
    end
    
    q_d = noiseStd/400;
    alphaForNoise = 2;
    
    maxNoiseLength=till/noiseUpdateFreq;
    
    %highPassFreq = [1 500];  %at 2000Hz
    highPassFreq = [2 ];  %at 2000Hz
    
    paramsUsed.q_d=q_d;
    paramsUsed.alphaForNoise=alphaForNoise;
    paramsUsed.maxNoiseLength=maxNoiseLength;
    paramsUsed.highPassFreq=highPassFreq;
    paramsUsed.shuffled1FNoise=shuffled1FNoise;
    paramsUsed.savedPrecomputedNoise=savedPrecomputedNoise;
       
    precomputedNoise=[];
    %fNamePrecompNoise = ['~/DFA/matlab/precompNoise/N1f_collapsed_' num2str(j) '_' num2str(highPassFreq(1)) '_' num2str(highPassFreq(2)) '_N' num2str(maxNoiseLength) '.mat'];
    
    %one for each unit in the network
    for kk=1:nrUnits        
        [filtNoise,rawNoise,b,a] = precompute1FNoise(maxNoiseLength, highPassFreq,q_d,alphaForNoise);
        if shuffled1FNoise
            precomputedNoise(kk,:) = precompute1FNoise_scramble(b,a,rawNoise);
        else
            precomputedNoise(kk,:) = filtNoise;
        end
    end
    disp(['sd of achieved 1/f noise is ' num2str(std(precomputedNoise(:))) ' requested is ' num2str(noiseStd) ' . mean is: ' num2str(mean(precomputedNoise(:)))]);
end

T_g = gNonLinOffset_val;
%gNonLinOffset(1:end-1)=0;    % no constant offset for the inhib
gNonLinOffset(1:end-1)=T_g;    % no constant offset for the inhib
gNonLinOffset(end)=0;    % const input to inhib

paramsUsed.gNonLinOffset=gNonLinOffset;
paramsUsed.T_g=T_g;

modeCompartment = 1; %1 single compartment, 2 several
paramsUsed.modeCompartment = modeCompartment;

%% ==== prepare overall weight matrix and collapsed state to allow speed optimization
Woverall = zeros(nrUnits,nrUnits);

% collapse all sub-connection matrices into large connectivity matrix
for k=1:length(M)
    Wlocal = M(k).W;
    inds =1+((k-1)*nrUnits_perMap): k*nrUnits_perMap;
    Woverall( inds,inds) = Wlocal;
end

%total state
Soverall = zeros(nrUnits,1);
inputOfMap_collapsed = zeros(nrUnits,1);

Toverall = reshape([M.Tall],1, length(M)*M(1).N)';
biasInpsOverall = reshape([M.biasInps],1, length(M)*M(1).N)';
overallInd_collapsed = reshape([M.overallInd],1, length(M)*M(1).N)';
gNonLinOffset_collapsed = repmat(gNonLinOffset_val, nrUnits,1);

Conn_collapsed=[];
for j=1:length(Conn)
    Conn_collapsed(j,:) = [  Conn(j).indFrom Conn(j).indTo ];
end

%% split weight matrices from constraint cells to pyramids in excitatory and inhibitory part
% this is later needed to deal with dendritic non-linearity,which is driven only by inhibitory input
%
WDb_to_pyr_inhib = WDb_to_pyr;
WDb_to_pyr_inhib(find(WDb_to_pyr_inhib>0))=0;
WDb_to_pyr_excit = WDb_to_pyr;
WDb_to_pyr_excit(find(WDb_to_pyr_excit<0))=0;

%% =====
inpToMapRep=[];
weightsMaskCached=[];
weightsNonZeroCached=[];
counterCorrect=0;
noiseUpdateCounter=0;
nrErrorsOverTime=[];
nrSwitches=[];
errorCounter=0;
winnerIDPrev=zeros(1,length(M));

if debugMode
    debugHistoryCounter=1;
    stateHistoryAll = zeros(nrUnits+N_DBs,till);
    shuntHistoryAll = zeros(nrUnits,till);
end

for i=2:till
    clMap=-1;
    %run each map

    if economicMode
        % dont keep state history
        currentStateInd=1;
        prevStateInd=1;
    else
        currentStateInd=i;
        prevStateInd=i-1;
    end
    
     
    %===== see if the noise needs to be updated
        if i-noiseUpdateFreq*floor(i/noiseUpdateFreq)==0 | i==2   % speed optimze mod
            if enableNoiseInp==3
                noiseUpdateCounter = noiseUpdateCounter+1;
                [inputOfMap_collapsed] = determineExtInput_optimized(enableNoiseInp,i, nrUnits, M(k).N, inputOfMap_collapsed, noiseMeanFrozen, noiseUpdateFreq,precomputedNoise(:,noiseUpdateCounter),noiseMean);
            else
                
                manual_bimodal_switch=0;   % 0 is off; 1 is on, for manual debugging of setting patterns of dendritic inputs
                
                if manual_bimodal_switch
                    noiseMean_all = ones(nrUnits,1)*noiseMean;
                    noiseStd_all  = ones(nrUnits,1)*noiseStd;
                    
                    noiseMean_all(1:2:end) = noiseMean_all(1:2:end)./2;
                    noiseStd_all(1:2:end) = noiseStd_all(1:2:end)./2;
                    
                    [inputOfMap_collapsed] = determineExtInput_optimized(enableNoiseInp,i, nrUnits, M(k).N, inputOfMap_collapsed, noiseMeanFrozen, noiseUpdateFreq,noiseStd_all,noiseMean_all);
                    
                else
                    [inputOfMap_collapsed] = determineExtInput_optimized(enableNoiseInp,i, nrUnits, M(k).N, inputOfMap_collapsed, noiseMeanFrozen, noiseUpdateFreq,noiseStd,noiseMean);
                end
            end
        end

    %inpToPyramids is negative (already inhibiting), thus multiply with -1 to make it positive
    
    %first step - calc input that the DBs provide to the pyramids

    % CLAMP if init state, fix some state values for a given amount of time. if disabled, start values are all zero    
    if i<=enableInitState
        %enableInitState: if >0, how long (how many iteration steps), initStateParams
       indsInit = find( initStateVals >= 0 );    % negative means disabled
       Soverall(indsInit) = initStateVals(indsInit);
    end
    
    if dendParams(1)==-1
        %if ==-1, dendritic non-linearity is disabled
        %original without dend nonlin
        inpToPyramids = WDb_to_pyr * DB_state(:,prevStateInd) ;
        Soverall = max( Soverall+ (-Soverall + Woverall*Soverall - Toverall + biasInpsOverall + max(inpToPyramids( overallInd_collapsed ) + ...
            inputOfMap_collapsed, gNonLinOffset_collapsed)  )*dt , 0 );
        
    else
        if debugMode
            inpToPyramids_fromInhib = WDb_to_pyr_inhib * DB_state(:,prevStateInd) ;  %neg constr cells
            inpToPyramids_fromExcit = WDb_to_pyr_excit * DB_state(:,prevStateInd) ;  %pos constr cells
            
            dendInput = fShunt(-1*inpToPyramids_fromInhib( overallInd_collapsed ), dendParams  );  % dendritic input
            Soverall = max( Soverall+ (-Soverall + Woverall*Soverall - Toverall + biasInpsOverall + ...
                (inputOfMap_collapsed+inpToPyramids_fromExcit) .* dendInput )*dt,0);
        else
            %collapse for speed reasons
            
            
            Soverall = max( Soverall+ (-Soverall + Woverall*Soverall - Toverall + biasInpsOverall + ...
                (inputOfMap_collapsed+WDb_to_pyr_excit * DB_state(:,prevStateInd)) .* fShunt(-1*WDb_to_pyr_fromInhib * DB_state(:,prevStateInd), dendParams  ) )*dt,0);
        end
        
    end
   
%     % CLAMP if init state, fix some state values for a given amount of time. if disabled, start values are all zero    
%     if i<=enableInitState
%         %enableInitState: if >0, how long (how many iteration steps), initStateParams
%        indsInit = find( initStateVals >= 0 );    % negative means disabled
%        Soverall(indsInit) = initStateVals(indsInit);
%     end

    %dB is pasted into below for speed reasons
    DB_state(:,currentStateInd) = max( DB_state(:,prevStateInd) + -DB_state(:,prevStateInd) + WPyr_to_db * Soverall - T_DB * dt, 0);   %apply max(x,0) rectification to the DBs
    
    %====end DB updating
    if debugMode
        debugHistoryCounter=debugHistoryCounter+1;
        stateHistoryAll(:,debugHistoryCounter) = [Soverall; DB_state(:,currentStateInd)];
        if dendParams(1)~=-1

                shuntHistoryAll(:,debugHistoryCounter) = [dendInput];
        end
    end
    
    %==== compute nr errors over time already here if in economicMode
    if economicMode
       if mod(i,100)==0  %do it every 100 steps

           %M=writeBack_state(M, Soverall);


           colMatchFoundLocal = checkColorSolutionOfSim_optimized(Conn_collapsed, Soverall, M, nrUnits_perMap );
           
           %[colMatchFoundLocal] = checkColorSolutionOfSim_optimized(Conn, Soverall );
           
           errorCounter= errorCounter+1;
           nrErrorsOverTime(errorCounter) = colMatchFoundLocal;
           
           % how many switches happened since last time
           %[~, winnerIDAll] = getAmpOfWinner(M, length(M), currentStateInd);

           %nrSwitches(errorCounter) = length( find( winnerIDAll ~= winnerIDPrev ));
           %winnerIDPrev = winnerIDAll;           
       end
    end
        
    %output debug messages and check for early termination
    if mod(i,1000)==0

        M=writeBack_state(M, Soverall);
        
        if graphType~=4
            [colMatchFound,errPercentage, errorNodes] = checkColorSolutionOfSim(Conn, M, M,0,[] );  % [] = use last timepoint available
        else
            [colMatchFound,errPercentage, errorNodes] = checkColorSolutionOfSim(Conn, M,M,0,[],2);
        end
        
        if isempty(biasInp)
            nrBiasErrors=0;
        else
            nrBiasErrors = getNrBiasErrors( M, length(M), biasInp,[]);
        end
        
        disp(['simulation step ' num2str(i) ' of ' num2str(till) ' mode=' num2str(enableNoiseInp) ' nrNodes=' num2str(length(M)) ' nrErrors=' num2str(colMatchFound) ' nrBiasErrors=' num2str(nrBiasErrors) ]);
    
        if colMatchFound
            counterCorrect=0;  %has to be correct twice consecutively
        end
        if terminateIfNoError==1 & ~colMatchFound & nrBiasErrors==0
            counterCorrect = counterCorrect+1;
            
            if counterCorrect >= counterCorrectRequired  %require X correct runs
                disp(['early terminate sim,no errors found anymore']);
                for k=1:length(M)
                    if ~economicMode
                        %M(k).S = M(k).S(:,1:i);
                        if keepHistory
                            M(k).inpHistory = M(k).inpHistory(:,1:i);
                            M(k).DBInpHistory = M(k).DBInpHistory(:,1:i);
                        end
                    end
                end
                
                % cut short the history appropriately
                stateHistoryAll = stateHistoryAll(:,1:debugHistoryCounter);   % remove zeros
                shuntHistoryAll = shuntHistoryAll(:,1:debugHistoryCounter);
                break;
            else
                disp(['all correct but not terminate yet']);
                
            end
        end
        
    end
    
end

M=writeBack_state(M, Soverall);
        
% determine the external input for this iteration
%k node number
%i integration step
%N nr units in each WTA
%
%
% enableNoiseInp: 1 frozen noise, 2 iid normal noise, 3 1/f noise
%
%
function [inp,inputOfMap, noiseChanged] = determineExtInput(enableNoiseInp,i, k, N, inputOfMap, noiseMeanFrozen,noiseUpdateFreq, noiseStd,noiseMean,precomputedNoiseTraces)
noiseChanged=0;
%== determine external input to be added

if enableNoiseInp
    switch(enableNoiseInp)

        case 1 %frozen noise only   
            if mod(i,noiseUpdateFreq)==0 || i==2      
                inputOfMap(1:end-1,k) = noiseMeanFrozen(:,k);   %constant input (frozen noise)
                noiseChanged=1;
            end
        case 2  %continously update noise
            if mod(i, noiseUpdateFreq)==0 || i==2    %initialize at very beginning (update moment is deterministic)
            %if mod(i, noiseUpdateFreq+round(rand*50) )==0 || i==2    %initialize at very beginning   (randomize exact update moment)


                %test with low-noise initial conditions
                %if i<500
                %    noiseStdReduced=noiseStd/100;
                %    inputOfMap(1:end-1,k) = randn(N-1,1)*noiseStdReduced+noiseMean;  %noisy external input   very small initial var
                %else
                
                
                %%==== testing, frozen noise initially
                
                
                % only update some subset, probabilisitc
                if i<100   % first one should update everything
                    inputOfMap(1:end-1,k) = randn(N-1,1)*noiseStd+noiseMean;   %noisy external input
                else
                    %indsUpdate = find( rand(1,N-1)<=1 );   % update probability (test,all)
                    %indsUpdate = find( rand(1,N-1)<=0.1 );   % update probability (only subset, randomly chosen); every 10th time (noiseUpdate=50)
                    indsUpdate = find( rand(1,N-1)<=0.2 );   % update probability (only subset, randomly chosen); every 5th time (noiseUpdate=100)
                    inputOfMap( indsUpdate, k) = randn( length(indsUpdate),1)*noiseStd + noiseMean;   %noisy external input
                end
                
                noiseChanged=1;
            end
            
        case 3  %continously update noise (1/f noise, pre-computed traces)
            if mod(i, noiseUpdateFreq)==0 || i==2    %initialize at very beginning (update moment is deterministic)

                    
                %if this mode is enabled, noiseStd is the precomputed noise trace for all units at this step
                %if i<100 
                 
                inputOfMap(1:end-1,k) = noiseStd + noiseMean;   %noisy external input
                
                %else
                %    indsUpdate = find( rand(1,N-1)<=0.2 );   % update probability (only subset, randomly chosen); every 5th time (noiseUpdate=100)
                %    inputOfMap( indsUpdate, k) = noiseStd(indsUpdate) + noiseMean;   %noisy external input
                %    
                %end
            end
            
        otherwise
            error('unsupported noise mode');
    end

    inp = inputOfMap(:,k);  %re-use the previous one if no new one was set
else
    inp=zeros(N, 1);
end

function [inputOfMap_collapsed] = determineExtInput_optimized(enableNoiseInp,i, nrUnitsTot, NperMap, inputOfMap_collapsed, noiseMeanFrozen, noiseUpdateFreq,noiseStd,noiseMean)

%noiseChanged=0;
%== determine external input to be added

if enableNoiseInp
    switch(enableNoiseInp)

        %case 1 %frozen noise only   
        %    if mod(i,noiseUpdateFreq)==0 || i==2      
        %        inputOfMap(1:end-1,k) = noiseMeanFrozen(:,k);   %constant input (frozen noise)
        %        noiseChanged=1;
        %    end
        
        case 2  %continously update noise
            %if mod(i, noiseUpdateFreq)==0 || i==2    %initialize at very beginning (update moment is deterministic)
            %if mod(i, noiseUpdateFreq+round(rand*50) )==0 || i==2    %initialize at very beginning   (randomize exact update moment)


                %test with low-noise initial conditions
                %if i<500
                %    noiseStdReduced=noiseStd/100;
                %    inputOfMap(1:end-1,k) = randn(N-1,1)*noiseStdReduced+noiseMean;  %noisy external input   very small initial var
                %else
                
                
                % testing with partially frozen noise
%                 if i<10000   % first one should update everything
%                 
%                     load('c:\temp\frozenVals3.mat');
%                     inputOfMap_collapsed=frozenVals;
%                     
%                     %inputOfMap_collapsed = randn(nrUnitsTot,1).*noiseStd+noiseMean;   %noisy external input
%                 
%                 else
%                     indsUpdate =1:nrUnitsTot;
%                     inputOfMap_collapsed( indsUpdate) = randn( length(indsUpdate),1).*noiseStd + noiseMean;   %noisy external input
%                 end

                %===DEFAULT(production)
                if i<100   % first one should update everything
                    inputOfMap_collapsed = randn(nrUnitsTot,1).*noiseStd+noiseMean;   %noisy external input
                else
                    indsUpdate =1:nrUnitsTot;
                    %indsUpdate = find( rand(1,nrUnitsTot)<=1 );   % update probability (test,all)
                    %indsUpdate = find( rand(1,N-1)<=0.1 );   % update probability (only subset, randomly chosen); every 10th time (noiseUpdate=50)
                    %indsUpdate = find( rand(1,nrUnitsTot)<=0.2 );   % update probability (only subset, randomly chosen); every 5th time (noiseUpdate=100)
                    inputOfMap_collapsed( indsUpdate) = randn( length(indsUpdate),1).*noiseStd + noiseMean;   %noisy external input
                end
                %=== end DEFAULT (production)
                
                
                % zero input to inhib units
                indsInhib = NperMap:NperMap:nrUnitsTot;
                inputOfMap_collapsed(indsInhib)=0;
                
                %noiseChanged=1;
            %end
            
         case 3  %continously update noise (1/f noise, pre-computed traces)
                  
                % inputOfMap(1:end-1,k) = noiseStd + noiseMean;   %noisy external input

                indsUpdate =1:nrUnitsTot;
                
                inputOfMap_collapsed( indsUpdate) = noiseStd + noiseMean;   %noisy external input
                 
                 
                % zero input to inhib units
                indsInhib = NperMap:NperMap:nrUnitsTot;
                inputOfMap_collapsed(indsInhib)=0;
                            
        otherwise
            error('unsupported noise mode');
    end

end


function M=writeBack_state(M, Soverall)
        
% write state back into maps
for k=1:length(M)
    inds = [1:M(k).N]+(k-1)*M(k).N;
    M(k).S = Soverall(inds);
end
