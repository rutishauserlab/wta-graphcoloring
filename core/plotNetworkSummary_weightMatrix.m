%
%plot weight matrix of the entire system (all WTAs, all DBs)
%
%urut/nov2014
function plotNetworkSummary_weightMatrix( Maps, WPyr_to_db, WDb_to_pyr, NrOfColors)

nrUnits_perMap = (NrOfColors+1);
nrDBcells = size(WPyr_to_db,1);
nrUnits = nrUnits_perMap*length(Maps)  + nrDBcells;

Woverall = zeros(nrUnits,nrUnits);

% collapse all sub-connection matrices into large connectivity matrix
for k=1:length(Maps)
    Wlocal = Maps(k).W;
    inds =1+((k-1)*nrUnits_perMap): k*nrUnits_perMap;
    Woverall( inds,inds) = Wlocal;
end

%add the DB cells
offset = nrUnits_perMap*length(Maps);
Woverall(1:offset, offset+1:end ) = WDb_to_pyr;
Woverall(offset+1:end, 1:offset ) = WPyr_to_db;

%load('colormap1.mat');

%=== to use imagesc to plot weight matrix
colorMapWeightMatrix=setmap;
colormap(colorMapWeightMatrix);

subplot(2,2,1);
imagesc(Woverall)
colorbar;
title('overall weight matrix');

subplot(2,2,2);
imagesc(WPyr_to_db);
title('pyr to DB (all pos)');
colorbar;

subplot(2,2,3);
imagesc(WDb_to_pyr);
title('DB to pyr (some neg, some pos)');
colorbar;

%% === to use dots to plot weight matrix
%subplot(2,2,4);

figure(15);

dotSize=4;
plotWeightMatrix_withDots(Woverall,colorMapWeightMatrix,dotSize);

disp('');


