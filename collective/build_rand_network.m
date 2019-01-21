%generate a random collective computing network
%
%connectionList1: inhib connectivity, basic network (WTAs)
%connectionList2: inhib connectivity (constraints)
%connectionListExcit: excitatory constraints
%
%if unitGrid/ unitList/connectionList1/connectionList2 is pre-set will be reused and not overwritten. provide as [] to generate new
%
%
function [unitGrid,unitList,connectionList1, connectionList2,connectionListExcit,W,indsInhib,indsExcit] = build_rand_network(gridSize,Poccupied,Pexcit, Pinhib2, PexcitCross, alpha1, alpha2, beta1,beta2, beta1_DB, beta2_DB, unitGrid, unitList, connectionList1, connectionList2)


% if preset use the same units
if isempty(unitGrid)
    
    unitGrid=zeros(gridSize,gridSize);
    unitList=[]; %unitNr xPos yPos unitType
    c=0;
    for xPos=1:gridSize
        for yPos=1:gridSize
            unitType=0;
            if rand<=Poccupied
                % a unit is placed at this position
                if rand<=Pexcit
                    unitType=1; %excitatory
                else
                    unitType=2; %inhibitory
                end
            end
            
            unitGrid(yPos,xPos) = unitType;
            
            if unitType>0
                c=c+1;
                unitList(c,:) = [c xPos yPos unitType];
            end
        end
    end
    
end

%nr of units in the network
N = length(find(unitGrid(:)>0));   

cConns1=0;
cConns2=0;
%connectionList1=[]; % fromExcit toInhib ;    basic connectivity
connectionList2=[]; % fromExcit toInhib ;    between-map connectivity

cConnsExcit=0;
connectionListExcit=[];

%make connections in the network
indsInhib = find(unitList(:,4)==2);   % these are inhib units
indsExcit = find(unitList(:,4)==1);   % these are inhib units

nrInhib = length(indsInhib);
nrExcit = length(indsExcit);

if isempty(connectionList1)
    for k=1:size(unitList,1)
        %for each excitatory unit, pick which inhib units it connects to
        if unitList(k,4)==1
            %pick one inhib unit
            R = randperm(nrInhib);
            indInhibToConnect = indsInhib(R(1));
            cConns1 = cConns1+1;
            connectionList1(cConns1,:) = [ unitList(k,1) unitList(indInhibToConnect,1) ];
        end
    end
end


% find which inhib connections are already taken

usedList = connectionList1(:,2);
indInhibToConnectUsed=[];
for k=1:length(usedList)
    indInhibToConnectUsed = unitList(find( unitList(:,1)== usedList(k)));
end


if isempty(connectionList2)
    for k=1:size(unitList,1)
        %for each excitatory unit, pick which inhib units it connects to
        if unitList(k,4)==1
            
            % with some probability also connect to other inhib units
            for jj=1:length(Pinhib2)
                if rand<Pinhib2(jj)
                    indsInhibRemain = setdiff( indsInhib, indInhibToConnectUsed);
                    if length(indsInhibRemain)>0
                        R2 = randperm( length(indsInhibRemain) );
                        
                        indInhibToConnect2 = indsInhibRemain(R2(1));
                        cConns2 = cConns2+1;
                        connectionList2(cConns2,:) = [ unitList(k,1) unitList(indInhibToConnect2,1) ];
                        
                        indInhibToConnect = indInhibToConnect2;
                    end
                end
            end
        end
    end
end
for k=1:size(unitList,1)
    %for each excitatory unit, pick which inhib units it connects to
    if unitList(k,4)==1
        %for each excitatory unit, decide whether to insert a alpha2 connection to an other excitatory unit
        if rand<PexcitCross
           %pick one inhib unit
           R = randperm(nrExcit);
           indExcitToConnect = indsExcit(R(1));
           cConnsExcit = cConnsExcit+1;
           connectionListExcit(cConnsExcit,:) = [ unitList(k,1) unitList(indExcitToConnect,1) ];
        end
    end
end

%% convert connectivity into weight matrix
W=zeros(N,N);
for k=1:size(unitList,1)
    if unitList(k,4)==1
        W(k,k) = alpha1;
    else
        W(k,k) = 0; 
    end
end

for jj=1:2
    if jj==1
        connectionList = connectionList1;
        beta2_toUse=beta2;
        beta1_toUse=beta1;
    else
        connectionList = connectionList2;
        beta2_toUse=beta2_DB;
        beta1_toUse=beta1_DB;
    end
    
    for k=1:size(connectionList)
        fromInd = connectionList(k,1);
        toInd = connectionList(k,2);
        
        W(toInd,fromInd) = beta2_toUse; %excitatory input to inhib unit
        W(fromInd,toInd) = -beta1_toUse; %inhibition
    end
end

%excitatory cross connections
for k=1:size(connectionListExcit)
    fromInd = connectionListExcit(k,1);
    toInd = connectionListExcit(k,2);
    
    W(toInd,fromInd) = alpha2; 
end
% 
% usedInhibUnits = unique(connectionList(:,2));
% nrConsToInhibUnits=[];
% for k=1:length(usedInhibUnits)
%     nrConsToInhibUnits(k) = length(find( usedInhibUnits(k) == connectionList(:,2)));
% end
% nrConsToInhibUnits
