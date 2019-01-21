%
%
%new version that uses DB cells
function [Maps, W, WPyr_to_db, WDb_to_pyr] = graphSim_setupWTAs_network_DB(NrOfColors, alpha1,beta1,beta2, stateMapping,stateTrans, Tdefault,Tinhib,Nnodes,beta2_DB, beta1_DB, beta2P_DB, beta1P_DB, modeNormalizeWeights, modeWeightNoise,graphType, debugPlotsEnabled )
%the internal connectivity (same for all)
N = NrOfColors+1;
W = diag( repmat(alpha1,1,N) ); 
W(N,N)=0; %no self-inhibition
W(N , 1:N-1) = beta2; %input to inhibitory neuron
W(1:N-1, N ) = -beta1; %output of inhibitory neuron

Maps=[];

overallInd=0;
for k=1:Nnodes
    Maps(k).W=W;
    Maps(k).ind=k;
    Maps(k).id = stateMapping(k,1);
    Maps(k).Tall = repmat(Tdefault,N,1);    %normal graph coloring
    
    %Maps(k).Tall = [Tdefault Tdefault+0.3 0]'; % max independent set
    Maps(k).Tall(end)  = Tinhib;
    Maps(k).N = size(W,1);
    Maps(k).indFromList = []; % list of indices into Conn, connections that input into this map
    
    Maps(k).overallInd = overallInd+1:overallInd+N;
    overallInd=overallInd+N;
    
end

nrDBMax=200;

DBcells=[];
DBcellsUsed=0;

Wdb = zeros(nrDBMax, N*Nnodes);  %from pyramids to DBs

%prepare a new stateTrans table that has every 

statsAllRings=[];

switch (graphType)
    case {1,2}
        
        stateTrans(:,4)=0; % already covered. 0 is not covered,1 is covered.
        %stateTrans is From, Symbol, To;   From/To notation is zero-based
        
        for k=0:Nnodes-1
            % notation is:    From(Self) ---> RemoteMaster --> RemoteSecondary
            % each secondary node needs to be connected to the master node
            % and to all other secondary nodes to be included in the same ring
            %
            %find the transitions that belong to this state (either from or to)
            [indsTrans,targetMapIDs] = graph_findTransition_forNode( stateTrans, k, 0 );
            while ~isempty(indsTrans) % same node can have multiple rings
                if length(indsTrans)>0
                    
                    remoteNodeMaster = targetMapIDs(1);
                    
                    nodeIDsForRing = [];  %at least 2, starting with self
                    nodeIDsForRing(1) = k; %self
                    nodeIDsForRing(2) = remoteNodeMaster; % first one
                    
                    %test candidates for secondary
                    if length(indsTrans)==1
                        %only one-create new ring for this one
                        stateTrans( indsTrans(1), 4)=1;
                        
                        
                    else
                        
                        for remoteNodeNr=2:length(indsTrans)   %see if more can be added to this ring
                            remoteNodeSecondary = targetMapIDs(remoteNodeNr);
                            
                            %==test if this node is connected to the first remote node
                            %get list of all connections from this remote node
                            [indsForRemoteNode,targetsRemoteNodeMaster] = graph_findTransition_forNode( stateTrans, remoteNodeSecondary, 0 );
                            
                            %see if the list includes the first remote node
                            
                            indsMasterToSecondary = [];
                            %test whether this candidate remote node is connected to he master as well as to all other already selected secondaries
                            secondaryQualifies=1;
                            secondaryTransitionsIncluded=[];
                            for jj=2:length(nodeIDsForRing)
                                
                                indsCandidateTransitions = find( nodeIDsForRing(jj) == targetsRemoteNodeMaster );
                                
                                if length( indsCandidateTransitions ) ~= 1
                                    secondaryQualifies=0;
                                    break;
                                else
                                    secondaryTransitionsIncluded = [secondaryTransitionsIncluded indsForRemoteNode(indsCandidateTransitions) ];
                                end
                            end
                            
                            if secondaryQualifies
                                %this remote node is connected to the remoteNodeMaster node & self
                                nodeIDsForRing=[nodeIDsForRing remoteNodeSecondary];
                                %mark as used
                                stateTrans( secondaryTransitionsIncluded, 4)=1; % master to secondary
                                stateTrans(indsForRemoteNode(find(targetsRemoteNodeMaster==k )) , 4)=1; % self to secondary
                            end
                            
                            %mark self to master as used
                            stateTrans( indsTrans(1), 4)=1;
                        end
                    end
                    
                    %insert the new ring
                    [DBcellsUsed, Wdb,statsOfRing] = graph_insertNewDBRing(NrOfColors, DBcellsUsed, Wdb, Maps, nodeIDsForRing, beta2_DB);
                    
                    statsAllRings = [statsAllRings; statsOfRing];
                else
                    %none left for this node - nothing needs to be done
                end
                
                [indsTrans,targetMapIDs] = graph_findTransition_forNode( stateTrans, k, 0 ); %see if there are transitions left of this node that belong to an other ring
                
            end
        end
        
    case 3
        %special rings for sudoku
        
        %nodeIDsForRing-1 to pass to graph_insertNewDBRing because is zero based whereas here it is 1 based
        
        if NrOfColors==9
            boxsize=3;
        end
        if NrOfColors==4
            boxsize=2;
        end
        %rows / columns
        for jj=1:NrOfColors
            %row
            nodeIDsForRing=[1:NrOfColors] + (jj-1)*NrOfColors;
            [DBcellsUsed, Wdb,statsOfRing] = graph_insertNewDBRing(NrOfColors, DBcellsUsed, Wdb, Maps, nodeIDsForRing-1, beta2_DB);
            statsAllRings = [statsAllRings; statsOfRing];
            
            %column
            nodeIDsForRing = [1:NrOfColors:Nnodes-NrOfColors+1]' + (jj-1);
            [DBcellsUsed, Wdb,statsOfRing] = graph_insertNewDBRing(NrOfColors, DBcellsUsed, Wdb, Maps, nodeIDsForRing-1, beta2_DB);
            statsAllRings = [statsAllRings; statsOfRing];
        end
        
        
        nodeMatrix=[];
        for rowNr=1:NrOfColors
            rowNodes=[1:NrOfColors] + (rowNr-1)*NrOfColors;
            nodeMatrix(rowNr,:) = rowNodes;
        end
        %boxes
        nrBoxes=NrOfColors/boxsize;
        for boxNrRow=1:nrBoxes
            for boxNrColumn=1:nrBoxes
                columnInds = [1:boxsize]+(boxNrColumn-1)*boxsize;
                rowInds = [1:boxsize]+(boxNrRow-1)*boxsize;
                boxNodes = nodeMatrix(rowInds,columnInds);
                nodeIDsForRing = boxNodes(:);
                
                [DBcellsUsed, Wdb,statsOfRing] = graph_insertNewDBRing(NrOfColors, DBcellsUsed, Wdb, Maps, nodeIDsForRing-1, beta2_DB);
                statsAllRings = [statsAllRings; statsOfRing];
                
            end
        end
    case 4  % max independent sets
        
        %add one DB for each constraint
        
        
        stateTrans(:,4)=0; % already covered. 0 is not covered,1 is covered.
        %stateTrans is From, Symbol, To;   From/To notation is zero-based
        
        statsPositiveConstr=[];
        rewireIndsAll=[];
        for jj=1:size(stateTrans,1)
            nodeIDsForRing=[stateTrans(jj,1) stateTrans(jj,3)];
            
            %insert the new ring (negative)
            NrOfColors_subset=1; %only constrain the first color
            [DBcellsUsed, Wdb,statsOfRing] = graph_insertNewDBRing(NrOfColors_subset, DBcellsUsed, Wdb, Maps, nodeIDsForRing, beta2_DB, 1);
            statsAllRings = [statsAllRings; statsOfRing];

            %insert two positive rings
            colFrom=2; % the second color            
            [DBcellsUsed, Wdb,statsOfRingP1,PyrInd1] = graph_insertNewDBRing_positive(colFrom, DBcellsUsed, Wdb, Maps, nodeIDsForRing(1), beta2P_DB );
            statsAllRings = [statsAllRings; statsOfRingP1];
            DBcellsUsed1=DBcellsUsed;
            
            [DBcellsUsed, Wdb,statsOfRingP2,PyrInd2] = graph_insertNewDBRing_positive(colFrom, DBcellsUsed, Wdb, Maps, nodeIDsForRing(2), beta2P_DB );
            statsAllRings = [statsAllRings; statsOfRingP2];                        
            DBcellsUsed2=DBcellsUsed;

            %the positive cells will point from one map to an other,so save for later re-wiring
            rewireInds = [ DBcellsUsed1 PyrInd2-1; DBcellsUsed2 PyrInd1-1];   %-1 to go to the previous color; Ind2/Ind1 flipped to go to the other node
            
            rewireIndsAll =[ rewireIndsAll; rewireInds ];
            
            statsPositiveConstr=[statsPositiveConstr; statsOfRingP1; statsOfRingP2];
            
            
        end
        
end




% for k=0:Nnodes-1
%     
%     
%     nrToCover = length(indsTrans)
% 
%     if nrToCover>0
%         nrCBsNeeded = NrOfColors;
%         
%         DBindsFrom = DBcellsUsed+1;
%         DBindsTo   = DBcellsUsed+nrCBsNeeded;
%         
%         DBinds = DBindsFrom:DBindsTo;
%         
%         DBcellsUsed = DBindsTo;
%         %
%         
%         %TODO-modify ring selection here
%         %
%         %1) members of one ring need to have all:all connectivity
% 
%         %2) for one node several rings can be created,till all are part of a ring
%         
%         % loop over all nodes
%         %    loop over all connections that lead to this node that have not been marked yet
%         %      pick one connected node and test for all the others whether it is connected to the one currently being tested. 
%         %      all those belong to the same ring. assign and mark as used.
%         %      repeat if unmarked connected nodes still exist
%         %
%         for colInd=1:NrOfColors   %over all colors
%             for jj=1:nrToCover    %nr of constraints that this ring has
%                 
%                 stateEndpoints = stateTrans( indsTrans(jj), [1 3]);
%                 targetMapInd = setdiff( stateEndpoints, k); %setdiff to exclude itself
%                 
%                 Wdb( DBinds(colInd), Maps(targetMapInd+1).overallInd(colInd) ) = beta2_DB;
%             end
%             
%             %also,include itself
%             Wdb( DBinds(colInd), Maps(k+1).overallInd(colInd) ) = beta2_DB;
%             
%         end
%     end
%     stateTrans(indsTrans,4)=1;  %set these to covered
%     
% end

WPyr_to_db = Wdb(1:DBcellsUsed,:);
WDb_to_pyr = WPyr_to_db';
WDb_to_pyr(find(WDb_to_pyr==beta2_DB)) = -beta1_DB;

%if there are positive,switch the to_pyr connections to the first color (max indp set); the above inserts a symmetric pDB cell that feeds back onto its
%own color; instead,it should feedback to the other color (max indp set can only work with 2 colors)
if exist('statsPositiveConstr')
    pDBCells = unique(statsPositiveConstr(:,1));
    if ~isempty(pDBCells)
        for k=1:length(pDBCells)
            
            DBInd = pDBCells(k);
            
            indOrig = find( WDb_to_pyr(:,DBInd)==beta2P_DB );
            
            %indNew = indOrig-1;
            
            %WDb_to_pyr(indNew,DBInd) = beta1P_DB;   % set to the previous color
            WDb_to_pyr(indOrig,DBInd) = 0;          %remove from this color
            
            
        end
        
        for k=1:size(rewireIndsAll,1)
            
            WDb_to_pyr( rewireIndsAll(k,2), rewireIndsAll(k,1) ) = beta1P_DB;
        end
    end
end

%add noise to the weights
if modeWeightNoise
    noiseStdProp=0.02;
    WDb_to_pyr = graph_addProportionalWeightNoise( WDb_to_pyr, noiseStdProp );
end

if modeNormalizeWeights
    %from DBs to pyramids
    WDb_to_pyr = graph_DBcells_normalizeWeights_row( WDb_to_pyr, -beta1_DB);    
    %from pyramids to DBs
    %WPyr_to_db = graph_DBcells_normalizeWeights_row( WPyr_to_db, beta2_DB);
end

if debugPlotsEnabled
    %plot the Pyr to DB weight matrices
    figure(50);
    subplot(2,2,1);
    imagesc(WPyr_to_db); xlabel('from: Maps');ylabel('to:DBs');colorbar;
    title('pyramids to DBs');
    subplot(2,2,2);
    imagesc(WDb_to_pyr); xlabel('from: DBs');ylabel('to: Maps');colorbar;
    title('DBs to pyramids');
    
    subplot(2,2,3);
    plot( statsAllRings(:,1), statsAllRings(:,2), 'o','MarkerSize',5);
    xlabel('ring nr (1 ring per color, 4 == 1 ring)');
    ylabel('members [node#]');
    xlim([0 max(statsAllRings(:,1))+1]);
    ylim([-1 max(statsAllRings(:,2))+1]);
    
    subplot(2,2,4);
    allWeightsUsed1 = WDb_to_pyr(find(WDb_to_pyr~=0));
    allWeightsUsed2 = WPyr_to_db(find(WPyr_to_db~=0));
    hist( [allWeightsUsed1 allWeightsUsed2],30 );
    title('DB weight distribution');
    
end

%===== after: TODO;


% %-- standard method: cross-inhibition
% Wconn = diag( repmat(gamma,1,N) );
% Wconn(end,end)=0; %no inhib to inhib
% 
% % Wconn(end-1,end-1)=0; %for max independent set
% 
% %-- alternative: cross-excitation?
% %Wconn = repmat(1,N,N) - eye(N);
% %Wconn(end,:)=0;
% %Wconn(:,end)=0;
% %Wconn = Wconn*gamma*-1;
% 
% Conn=[];
% c=0;
% for k=1:Nedges
%     %each edge is bidirectional,add two entries for each
%     for j=1:2        
%         if j==1
%             idFrom = stateTrans(k,1);
%             idTo   = stateTrans(k,3);
%         else
%             idFrom = stateTrans(k,3);
%             idTo   = stateTrans(k,1);
%         end
%         
%         c=c+1;
%         Conn(c).idFrom = idFrom;  
%         Conn(c).idTo = idTo;
% 
%         % ids in conn are zero based. to convert to inds into Maps, add +1
%         Conn(c).indFrom=Conn(c).idFrom+1;
%         Conn(c).indTo=Conn(c).idTo+1;
% 
%         Conn(c).Wconn = Wconn;
%     end
% end
% 
% %add direct inds into Maps for speedup
% for k=1:Nnodes
% 
%     indToSearch = Maps(k).ind;
%     
%     indFromList=[];  %list of all Maps that give input to map k
%     for j=1:length(Conn)
%        
%         if Conn(j).indTo == indToSearch
%             indFromList = [ indFromList j]; %Conn(j).indFrom];
%         end
%     end
%     
%     Maps(k).indFromList = indFromList;
% end
% 
% %stats for this graph
% nrEdgesEach=[];
% for k=1:length(Maps)
%     nrEdgesEach(k) = length( Maps(k).indFromList );
% end
% 
% 
% %normalize
% if normConnsEnabled
%     for k=1:length(Maps)
%         
%         indFromList = Maps(k).indFromList;
%         
%         n = length(indFromList);
%         
%         WconnNorm = (Wconn./n)*normFact;  %this multiple of gamma is allowed
%         
%         for j=1:n
%             Conn( indFromList(j) ).Wconn = WconnNorm;
%         end
%     end
%     
% end