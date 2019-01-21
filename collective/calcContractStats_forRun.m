%for a run of a WTA network, compute various metrics describing the dynamics
%
%maxEig: max eigenval
%minEig: min eigenval
%trEig: sum of all eig vals (=trace)
%condNr4: condition nr
%
function [maxEig,minEig,trEig,condNr4] = calcContractStats_forRun(nrSteps, DD, N, Wtotal, loadTerms, winnerID )



winnerLocation = zeros(1,N);
winnerLocation(end)=1;
winnerLocation(winnerID)=1;

%winnerLocation(3)=1;

[V,D] = eig( (diag(winnerLocation)*Wtotal)+loadTerms );  %eigenvectors
ThetaInv = V;
Theta = inv(V);

condNr4=[];
trEig=[];
maxEig=[];
minEig=[];
eigVectmax=[];
eigVectmin=[];
eigVectSigns=[]; %count of the nr of eigenvectors with common sign
maxVectDist=[];
max2ndEig=[]; %second largest eig val
minAbsEig=[];

allEigs=[];
maxEigSkew=[];

for i=2:nrSteps
    
    
    Jtmp = ( diag(DD(:,i)) * Wtotal ) + loadTerms;
    %Jtmp(end,end)=-1;
    
    Jtrans = Theta * Jtmp * ThetaInv;  %apply theta transform    
    
    %Jtrans=Jtmp;

    % 

%    eigs_skew = eig( getSkewPart( Jtmp ) );
    
    eigs = eig( getSymPart(Jtrans)  );    % symmetric part, theta transformed
    %[V,D] = eig( getSymPart(Jtrans)  );    % symmetric part, theta transformed

    %eigs = eig( getSymPart(Jtmp)  );    % symmetric part, raw jacobian
    %[V,D] = eig( getSymPart(Jtmp)  );    % symmetric part, raw jacobian
% 
%     [V,D] = eig( (Jtmp)  );    % raw jacobian
%     
%     [maxVal,mInd]=max(diag(D));
%     [~,minInd]=min(diag(D));
% 
%     
%     eigenVectTrans = ThetaInv*V*Theta;
%     
%     maxVectDist(i) = ctranspose(eigenVectTrans(:,mInd))*eigenVectTrans(:,mInd);
% 
%     %maxVectDist(i) = ctranspose(V(:,mInd))*V(:,mInd);   % is, by definition, always one (norm)
%     
%     %eigVectmax(:,i) = eigenVectTrans(:,mInd);
%     
%     eigVectmax(:,i) = (V(:,mInd));
%     eigVectmin(:,i) = (V(:,minInd));
%     
%     % count common/mixed mode eigenvects
%     eigVectSigns(i,1:3) = [0 0 0];  % total, with positive eigenval, with negative eigenval
%     
%     eigValsForV=diag(D);
%     
%     for k=1:size(V,2)
%         
%         sOfV = sign(V(:,k));
%         
%         if abs(sum(sOfV))==size(V,2)
%             eigVectSigns(i,1) = eigVectSigns(i,1) +1;
%             
%             if eigValsForV(k)>0
%                 eigVectSigns(i,2) = eigVectSigns(i,2) +1;
%                 
%             else
%                 eigVectSigns(i,3) = eigVectSigns(i,3) +1;
%                 
%             end
%         end
%         
%     end
    
    %eigs = eig( (Jtmp)  );    % symmetric part, raw jacobian
    
    %[V,D] = eig( (Jtmp)  );    % symmetric part, raw jacobian
    
    %diag(D)
    
    %if i==2500
    %    keyboard;
    %end
    
    
    %V(:,1)
    
    %eigs=eig(Jtmp);  % take eigenvalues of the raw functional jacobian
    
    
    %apply theta, take hermitian part
    %eigs = eig( getSymPart(Jtrans)  );
    
    %[V,D] = eig( getSymPart(Jtrans)  );    % symmetric part, raw jacobian
    
    %diag(D)
    
    
    
    maxEig(i) = max(real(eigs));
    minEig(i) = min(real(eigs));

    condNr4(i) = max(abs(real(eigs)))./min(abs(real(eigs)));
    
    trEig(i) = sum(eigs);
    


    
    %second-max
    
    %max2ndEig(i) = max(setdiff( real(eigs), max(real(eigs))));

    %min abs
    %minAbsEig(i) = min(abs(real(eigs)));
    
    %allEigs(i,:) = sort(eigs,1,'descend');

    %eigs_skew_imag = imag(eigs_skew);
    
    %maxVals=max(eigs_skew_imag);
    %maxEigSkew(i,:) = [ maxVals(1) length(find(eigs_skew_imag~=0)) ];
end
