%
% simplified version of recurrent map circuit. 
% 2 WTAs, 2 winners each
% 2 DB cells that enforce "not same" between the two
%
% ueli rutishauser, urut@caltech.edu
%===================================

% == this is figure 2 in the sudoku paper.


%% setup
alpha1=1.2;
alpha2=0;

beta1=3;     % 
beta2=0.25;

beta1_DB=1.2;  % at least (3+beta2)/2
beta2_DB=0.3;  % max 0.5*(-4+alpha1+beta1+2*beta1_DB-beta2)

gamma=0;

G=1; %load
I=6.01; %external input 
s=1;

%== verification
T=0.0; %threshold for excitatory units
T_DB=0;   % high positive value disables
Tinhib=0;  % negative value is constant input

g=1/(1+beta1*beta2-alpha1);
gainNoDB =  g * ( 1+ beta1*Tinhib )  %gain of the entire system with no DBs

gDB = 1/(1+beta1*beta2 + beta1_DB*beta2_DB -  alpha1);
gainWithAllDB = gDB * ( 1 + beta1*Tinhib + beta1_DB*T_DB)

gainPredicted = gainWithAllDB;

%--calculate eigenvalues
W=[alpha1 alpha2 ;
   alpha2 alpha1];

J=[beta2 beta2 ];

WtoDB =[0        beta2_DB 0        beta2_DB;
        beta2_DB 0        beta2_DB 0];
nrDBs = size(WtoDB,1);

%onto themselfs
%WtoDB=[beta2_DB 0 0 0;
%       0        0 beta2_DB 0];
  
WfromDB=WtoDB';

WfromDB(find(WfromDB==beta2_DB))=beta1_DB;

%WtoDB = graph_DBcells_normalizeWeights_column( WtoDB, beta1_DB);   
%WfromDB = graph_DBcells_normalizeWeights_row( WfromDB, beta1_DB); 

disableDBInput = 0;
turnOffDBToLoosers=0; %after some time,turn of DB input to loosers
turnOffDBToWinners=0; %after some time,turn of DB input to winners
delayDBOnset=0; % if >0,only turn on DB input after this time

%== simulation
till=15000;
N=2; %how many excitatory

x1=ones(N,till)*0;
y1=zeros(1,till);

x2=ones(N,till)*0;
y2=zeros(1,till);

DB = zeros(nrDBs,till);

inpHistory=zeros(1,till);

%2*N is nr of target pyramids
fromDBHist=zeros(2*N,till);
initDBState=zeros(1,2*N);

inputAmp=[I;I*0.9];

%Fixed
inputAmp1=[I*0.6;I*0.7];
inputAmp2=[I*0.8;I*0.9];

%Randomized
%inputAmp1 = randn(N,1)+I;
%inputAmp2 = randn(N,1)+I;

inputAmp1orig = inputAmp1;
inputAmp2orig = inputAmp2;

%inputAmp1=[I*0.5;I*0.97;I*0.55];
%inputAmp2=[I*0.91;I*0.98;I*0.93];

%inputAmp2=inputAmp1*1.05;


%no competition
%inputAmp1=[I*1.01;I*0];
%inputAmp2=[I*1.02;I*0];

addNoise = 0; 
noiseUpdateFreq=100;
noiseStd=0.01;

inp=[0;0];

inpDefault=inp;
inp1=inp; 
inp2=inp;

maxCutoff=0;   %normal
%maxCutoff=-1000; %to remove the non-linearity

%perform the numerical integration

%to start with different initial values (next 2 lines)
%x1(:,1)=[0 0 0];
%y1(1)=0;

aI=0; %recurrent of inhibition (experimental)
dt=.01;


%% run simulation
DD=[];
inpHistory=[];
for i=2:till

    %inpHistory(i)=inp(2);
    
    inpHistory(i,:) = [inp1; inp2]';
    
    %uses this funct: x = maxNonlin( x, s, cutoff)

    S = [ x1(:,i-1); x2(:,i-1)];   %state, concat all excitatory units.

    if disableDBInput
        fromDB=initDBState';
    else
        fromDB=initDBState';
        if delayDBOnset==0 | (delayDBOnset & i>delayDBOnset)
            fromDB = WfromDB*DB(:,i-1);         
        end
        if i>7000 & turnOffDBToWinners
            fromDB(find(S>0.01))=0;
        end
        if i>7000 & turnOffDBToLoosers
            fromDB(find(S<0.01))=0;
        end
    end        
    
    fromDBHist(:,i)=fromDB;
    %map1
    dx1 = -x1(:,i-1)*G + maxNonlin( -fromDB(1:length(fromDB)/2) + W*x1(:,i-1) - beta1*y1(i-1) + inp1 - T, s, maxCutoff);
    dy1 = -y1(i-1)*G   + maxNonlin( y1(i-1)*aI + J*x1(:,i-1) - Tinhib, s, maxCutoff);
    
    
    %map2
    dx2 = -x2(:,i-1)*G + maxNonlin( -fromDB(length(fromDB)/2+1:end) + W*x2(:,i-1) - beta1*y2(i-1) +inp2 - T, s, maxCutoff);
    dy2 = -y2(i-1)*G   + maxNonlin( y2(i-1)*aI + J*x2(:,i-1) - Tinhib, s, maxCutoff);

    dDB = - DB(:,i-1) + maxNonlin( WtoDB * S - T_DB , s, maxCutoff);

    
    %euler integration (dt)
    dx2 = (dx2)*dt ;
    dy2 = (dy2)*dt ;
    dx1 = (dx1)*dt ;
    dy1 = (dy1)*dt ;
    dDB = (dDB)*dt;
    
    x1(:,i) = x1(:,i-1)+dx1;
    y1(i)   = y1(i-1)+dy1; 
    x2(:,i) = x2(:,i-1)+dx2;
    y2(i)   = y2(i-1)+dy2;
  
    DB(:,i) = DB(:,i-1)+dDB;
    
    %disable DBs for a while for testing
    %if i<6000
    %    DB(:,i)=0;
    %end
    
    %old
    %x1(:,i) = max( s * ( x1(:,i-1)+dx1 ), maxCutoff);
    %y1(i)   = max( s * ( y1(i-1)+dy1   ), maxCutoff);    
    %x2(:,i) = max( s * ( x2(:,i-1)+dx2 ), maxCutoff);
    %y2(i)   = max( s * ( y2(i-1)+dy2   ), maxCutoff);
    
    if i>2000 & i<14000
        inp1=inputAmp1;
        inp2=inputAmp2;
        
        if addNoise
                   if mod(i, noiseUpdateFreq)==0 || i==2
                   
                       inputAmp1 = inputAmp1orig + randn(N,1)*noiseStd;
                       inputAmp2 = inputAmp2orig + randn(N,1)*noiseStd;
                       
                   end 
        end
    else
%         if i>8000 & i<8500
%             inp1=inputAmp1;
%             inp1(2)=inp1(2)*1.8;
%         else
%             if i>8500 & i<9500
%                 inp1=inputAmp1;
%             else            
            
                inp1=inpDefault;
                inp2=inpDefault;
%            end
%        end
    end
    
    

    actVect=[x1(:,i); y1(i); x2(:,i); y2(i)];
    DD(:,i) = [actVect>0.01];  %deriviative of the non-linearity
    
end

%% permitted/hidden sets analysis

Wtotal = zeros(N+1,N+1);
Wtotal(1:N,1:N) = W;
Wtotal(end,1:N) = J;
Wtotal(1:N,end) = -beta1;

loadTerms=diag(repmat(-G,1,N+1));

DD_subsetToPlot1 = DD(1:N+1,:); % which Node to plot (one entire WTA)
plotMode=1;
[allSets_sorted,trOfSet_sorted] = plotPermittedHiddenSets( Wtotal, loadTerms, DD_subsetToPlot1,1, []  );


%[DD_matchedSet1, totNrSets, allSets_sorted, trOfSet_sorted, indsPermitted_sorted, indsForbidden_sorted] = plotPermittedHiddenSets_prepare(  Wtotal, loadTerms, DD_subsetToPlot1 );

%DD_subsetToPlot2 = DD(N+2:end,:); % which Node to plot (one entire WTA)
%[DD_matchedSet2, totNrSets, allSets_sorted, trOfSet_sorted, indsPermitted_sorted, indsForbidden_sorted, eigOfSet_sorted] = plotPermittedHiddenSets_prepare(  Wtotal, loadTerms, DD_subsetToPlot2 );

%timeRange=2010:14000;


%enable to see plots
%plotPermittedHiddenSets_pairwise( DD_matchedSet1(timeRange), DD_matchedSet2(timeRange), indsPermitted_sorted, indsForbidden_sorted, totNrSets, trOfSet_sorted,eigOfSet_sorted);


%setNr_map1 = DD_matchedSet1(timeRange);
%setNr_map2 = DD_matchedSet2(timeRange);

%return;

%% plot the overall weight matrix

Wall=zeros( (N+1)*2 + 2, (N+1)*2 + 2);
Wall(1:N,1:N)=W;
Wall(N+1,1:N) = J;
Wall(1:N,3) = -beta1;

Wall( 3+[1:N],3+[1:N])=W;
Wall(3+[N+1],3+[1:N]) = J;
Wall(3+[1:N],3+3) = -beta1;



%WtoDB =[0        beta2_DB 0        beta2_DB;
%        beta2_DB 0        beta2_DB 0];
%nrDBs = size(WtoDB,1);

%onto themselfs
%WtoDB=[beta2_DB 0 0 0;
%       0        0 beta2_DB 0];
  
%WfromDB=WtoDB';

%WfromDB(find(WfromDB==beta2_DB))=beta1_DB;

%WtoDB

Wall( 1, 2*(N+1)+1) = -beta1_DB;
Wall( 2, 2*(N+1)+2) = -beta1_DB;

Wall( 4, 2*(N+1)+1) = -beta1_DB;
Wall( 5, 2*(N+1)+2) = -beta1_DB;


Wall( 2*(N+1)+1,1) = beta2_DB;
Wall( 2*(N+1)+2,2) = beta2_DB;

Wall( 2*(N+1)+1,4) = beta2_DB;
Wall( 2*(N+1)+2,5) = beta2_DB;

%Wall(7:8,:) = Wall(:, 7:8)';

figure(20);
subplot(2,2,1);
imagesc(Wall);
colorbar;

xlabel('unit nr [from]');
ylabel('unit nr [to]');

subplot(2,2,2);
colorMapWeightMatrix=setmap();
plotWeightMatrix_withDots(Wall,colorMapWeightMatrix,20)
xlim([0.5 8.5]);
ylim([0.5 8.5]);

%% testing overall matrix, eigenvalues
% 
% W2=zeros( (N+1)*2 + 1, (N+1)*2 +1);
% W2(1:N+1,1:N+1)=Wtotal;
% W2(N+2:end-1,N+2:end-1)=Wtotal;
% 
% DB1_ind=9;
% W2(DB1_ind,1)=beta2_DB;
% W2(DB1_ind,5)=beta2_DB;
% W2(1,DB1_ind) = beta1_DB;
% W2(5,DB1_ind) = beta1_DB;
% 
% 
% switchMatrix=[0 1 0 1 1 0 0 1 1];
% 
% switchMatrix=[1 0 0 1 1 0 0 1 1];
% 
% Jtmp2 = diag(switchMatrix) * W2 - eye(size(W2,1));
% 
% [V,D] = eig( Jtmp2 );  %eigenvectors
% 
% ThetaInv = V;
% Theta = inv(V);
% 
% Jtrans = Theta * Jtmp2 * ThetaInv;  %apply theta transform    
% 
% eigs(getSymPart(Jtrans))

%%
figure(2);
simpleRecurrence_plotColumn_twoWinners(till, x1, y1, inputAmp1, [1 4 7], 'map1');
simpleRecurrence_plotColumn_twoWinners(till, x2, y2, inputAmp2, [2 5 8], 'map2');

hold on
line([0 till],[gainPredicted gainPredicted],'color','k');
hold off

subplot(3,3,3);
x=1:till;
plot( x, DB(1,:), 'r', x, DB(2,:), 'g');
legend('DB1','DB2');
title('DB cells');

subplot(3,3,6);
plot( x, fromDBHist);
title('input from DB to pyramids');


winnerAmp = max(x1(:,8000));
if ~exist('totHist')
    totHist=[];
end

totHist( size(totHist,1)+1,:) = [I winnerAmp T_DB Tinhib];

%inputs
subplot(3,3,9);

plot( x, inpHistory(:,1:4) );
title('inputs');
legend('M1-x1','M1-x2','M2-x1','M2-x2');


enforceCommonAxisRange( [], 1, [0 7000],[] )

%% plot gain
% noDB_TDB=5;
% indsDB1=find(totHist(:,3)==noDB_TDB); %no DB
% TinhibVals1 = unique( totHist(indsDB1,4) );
% 
% withDB_TDB=0;
% indsDB2=find(totHist(:,3)==withDB_TDB); %no DB
% TinhibVals2 = unique( totHist(indsDB2,4) );
% 
% figure(150);
% c=0;
% hs=[];
% for jj=1:2
% subplot(2,2,1);
% 
% if jj==1
%     TinhibVals_toLoop=TinhibVals1;
%     DBVal_toUse = noDB_TDB;
%     markerTyp='-';
% else
%     TinhibVals_toLoop=TinhibVals2;
%     DBVal_toUse = withDB_TDB;
%     markerTyp='--';
%     
% end
% 
% 
% 
% for k=1:length(TinhibVals_toLoop)
%     
%         inds=find(totHist(:,4)==TinhibVals(k) & totHist(:,3)==DBVal_toUse );
%     
%     if k>1 
%         hold on
%     end
%     
%     c=c+1;
%     hs(c)=plot( totHist(inds,1), totHist(inds,2), [rotatingColorCode(k) 'x-']);
%     
%     legendStrs{c}=['Tinhib=' num2str(TinhibVals_toLoop(k)) ' TDB=' num2str(DBVal_toUse)];
% end
% 
% 
% end
% hold off
% legend(hs,legendStrs);
% title('add constant to T_inhib');
% xlabel('input to winner');
% ylabel('output of winner');
%% == plot result
% figure(1);
% subplot(2,2,1)
% plot( [1:till]+50, x1(1,:), 'g', [1:till]+100, x1(2,:), 'b', [1:till]+150, x1(3,:), 'm', 1:till, x2(1,:), 'y', 1:till, x2(2,:), 'r', 1:till, inpHistory, 'k');
% title(['excitatory neurons. predictions: ampInp=' num2str(ampliRec2) ' ampMem=' num2str(memSteadyState) ]);
% 
% legend('x1(1)','x1(2)','x1(3)', 'y1(1)', 'y1(2)','external input to x1');
% xlabel('time [integration steps]');
% ylabel('activity');
% 
% text(6000,2,['\alpha=' num2str(alpha1) ' \beta_1=' num2str(beta1) ' \beta_2=' num2str(beta2) ' \gamma=' num2str(gamma) ' T=' num2str(T)]);
% 
% subplot(2,2,2)
% plot(1:till, y1, 'b', 1:till, y2, 'r');
% title('inhibitory neurons');
% xlabel('time [integration steps]');
% ylabel('activity');
% 
% subplot(2,2,3)
% plot(1:till-1, diff(x1(2,:)), 'b', 1:till-1, diff(x2(2,:)), 'r');
% title('derivate excitatory neurons');
% xlabel('time [integration steps]');
% ylabel('activity');
% 
% subplot(2,2,4)
% plot(1:till-1, diff(y1), 'b', 1:till-1, diff(y2), 'r');
% title('derivative inhibitory neurons');
% xlabel('time [integration steps]');
% ylabel('activity');
% 

%=====
%vector fields
%figure(10);
%plotDynamics(xout, fromPhase,toPhase, neuronToPlot, binsize, colCode)

