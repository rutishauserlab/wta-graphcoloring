%
%plot the activity of all units, as a function of time, in the entire system (all WTAs, all DBs)
%
%urut/nov2014
function plotNetworkSummary_activity( Maps, WPyr_to_db, WDb_to_pyr, NrOfColors, Mhistory, DBhistory, tillPlot)

nrUnits_perMap = (NrOfColors+1);
nrDBcells = size(WPyr_to_db,1);
nrUnits = nrUnits_perMap*length(Maps)  + nrDBcells;

Soverall = zeros(nrUnits, tillPlot);

% collapse all sub-connection matrices into large connectivity matrix
for k=1:length(Mhistory)
    Wlocal = Maps(k).W;
    inds =1+((k-1)*nrUnits_perMap): k*nrUnits_perMap;
    Soverall( inds,:) = Mhistory(k).S(:,1:tillPlot);
end

%add the DB cells
offset = nrUnits_perMap*length(Maps);

Soverall(offset+1:end,:) = DBhistory(:,1:tillPlot);


imagesc(Soverall)
colorbar;

xlabel('time');
ylabel('unit nr');
title('activity');

