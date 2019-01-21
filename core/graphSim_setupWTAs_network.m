function [Maps,Conn,W,Wconn,nrEdgesEach] = graphSim_setupWTAs_network(NrOfColors, alpha1,beta1,beta2,gamma,stateMapping,stateTrans, Tdefault,Tinhib,Nnodes,Nedges, normConnsEnabled, normFact )
%the internal connectivity (same for all)
N = NrOfColors+1;
W = diag( repmat(alpha1,1,N) ); 
W(N,N)=0; %no self-inhibition
W(N , 1:N-1) = beta2; %input to inhibitory neuron
W(1:N-1, N ) = -beta1; %output of inhibitory neuron

Maps=[];
for k=1:Nnodes
    Maps(k).W=W;
    Maps(k).ind=k;
    Maps(k).id = stateMapping(k,1);
    Maps(k).Tall = repmat(Tdefault,N,1);    %normal graph coloring
    %Maps(k).Tall = [Tdefault Tdefault+0.3 0]'; % max independent set
    Maps(k).Tall(end)  = Tinhib;
    Maps(k).N = size(W,1);
    Maps(k).indFromList = []; % list of indices into Conn, connections that input into this map
end

%-- standard method: cross-inhibition
Wconn = diag( repmat(gamma,1,N) );
Wconn(end,end)=0; %no inhib to inhib

% Wconn(end-1,end-1)=0; %for max independent set

%-- alternative: cross-excitation?
%Wconn = repmat(1,N,N) - eye(N);
%Wconn(end,:)=0;
%Wconn(:,end)=0;
%Wconn = Wconn*gamma*-1;

Conn=[];
c=0;
for k=1:Nedges
    %each edge is bidirectional,add two entries for each
    for j=1:2        
        if j==1
            idFrom = stateTrans(k,1);
            idTo   = stateTrans(k,3);
        else
            idFrom = stateTrans(k,3);
            idTo   = stateTrans(k,1);
        end
        
        c=c+1;
        Conn(c).idFrom = idFrom;  
        Conn(c).idTo = idTo;

        % ids in conn are zero based. to convert to inds into Maps, add +1
        Conn(c).indFrom=Conn(c).idFrom+1;
        Conn(c).indTo=Conn(c).idTo+1;

        Conn(c).Wconn = Wconn;
    end
end

%add direct inds into Maps for speedup
for k=1:Nnodes

    indToSearch = Maps(k).ind;
    
    indFromList=[];  %list of all Maps that give input to map k
    for j=1:length(Conn)
       
        if Conn(j).indTo == indToSearch
            indFromList = [ indFromList j]; %Conn(j).indFrom];
        end
    end
    
    Maps(k).indFromList = indFromList;
end

%stats for this graph
nrEdgesEach=[];
for k=1:length(Maps)
    nrEdgesEach(k) = length( Maps(k).indFromList );
end


%normalize
if normConnsEnabled
    for k=1:length(Maps)
        
        indFromList = Maps(k).indFromList;
        
        n = length(indFromList);
        
        WconnNorm = (Wconn./n)*normFact;  %this multiple of gamma is allowed
        
        for j=1:n
            Conn( indFromList(j) ).Wconn = WconnNorm;
        end
    end
    
end