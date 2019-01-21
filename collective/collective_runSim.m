%
% run a simulation of a collective comp network
%
% nrSteps - how many time steps to simulate
% DD - is unit active/inactive
% x - state of the network
%
%urut/june14
function [x,inpHistory, DD] = collective_runSim(N, W, G, T, nrSteps, params, inputAmp, inputAmp2, inpDefault, noiseEnabled )

noiseUpdateFreq=100;

inpHistory = zeros(N,nrSteps);
x = zeros(N,nrSteps);
DD = zeros(N,nrSteps);
            
noiseInp = zeros(1,N);

inp = inpDefault;

SetTriggered=0;

for i=2:nrSteps    
    inpHistory(:,i) = inp;

    %euler integration (dt)
    dx = -x(:,i-1).*G + maxNonlin( W*x(:,i-1) + inp - T, params.s, params.maxCutoff);

    x(:,i) = x(:,i-1)+dx*params.dt;

    DD(:,i) = [x(:,i)>0.01];  %deriviative of the non-linearity
    
    if i>params.inpOnTime(1,1) & i<params.inpOnTime(1,2)
        inp = inputAmp;

        %if i>params.inpOnTime(1,1)+500
        %    inp(1)=inp(1)-4;
        %end
    
        % see if certain set is visited, then adjust input accordingly
        
        if params.setConditionalBranchingEnabled & ~SetTriggered
            if isequal(DD(:,i), params.setConditional )
                %set 6
                %inp(4)=inp(4)-4;
                %inp(1)=inp(1)+4;
                inputAmp(params.setConditionalInpNr(1)) = inputAmp(params.setConditionalInpNr(1)) + params.setConditionalInpNr(2);
                inp=inputAmp;
                SetTriggered=1;
            end
            
        end
        
        if noiseEnabled & mod(i,noiseUpdateFreq)==0      % update the noise only every X steps
            noiseInp = randn(1,N)./1;
            noiseInp(end) = 0; % no input to inhibitory neuron
        end
        
        inp = inp+noiseInp'; % add the noise every time
    else
        
        if i>params.inpOnTime(2,1) & i<params.inpOnTime(2,2)
            inp = inputAmp2;
            
        else
            inp = inpDefault;
            
        end
    end
    
end
