%
%
%
function collective_plotOneRun( x, inpHistory, nrSteps, unitsToPlot, maxEig, minEig, trEig, xlimCommon  )

t = 1:nrSteps;

%== plot result

%=========== activity on units
subplot(3,3,1)
strs=[];
hs=[];
for k=1:length(unitsToPlot)
    if k>1
        hold on
    end
    hs(k)=plot(t, x(unitsToPlot(k),1:nrSteps), rotatingColorCode(k) ); 
    strs{k}=['Unit ' num2str(unitsToPlot(k))];
end
hold off
legend(hs, strs);

title('units');


xlim(xlimCommon);

subplot(3,3,2);

plot( t, maxEig, 'k-', t, minEig, 'k--');
legend('\lambda_{max} (F_s)','\lambda_{min}(F_s)');
ylabel('eigenvalue');
xlim(xlimCommon);



subplot(3,3,3);

plot( t, trEig, 'k-');
ylabel('divergence');
xlim(xlimCommon);

%============ inputs to units
subplot(3,3,4);
strs=[];
hs=[];
for k=1:length(unitsToPlot)
    if k>1
        hold on
    end
    hs(k)=plot(t, inpHistory(unitsToPlot(k),1:nrSteps), rotatingColorCode(k) ); 
    strs{k}=['Unit ' num2str(unitsToPlot(k))];
end
hold off
legend(hs, strs);

%--- inputs - zoom-in
subplot(3,3,5);
strs=[];
hs=[];
for k=1:length(unitsToPlot)
    if k>1
        hold on
    end
    
    data = inpHistory(unitsToPlot(k),1:nrSteps);
    %indsOff=find(data==0);
    %indsOn=find(data>0);
    %data(indsOff) = data(indsOff)+1.5;
    %data(indsOn) = data(indsOn)-mean(data(indsOn))*0.8;
    
    hs(k)=plot(t, data, rotatingColorCode(k) ); 
    strs{k}=['Unit ' num2str(unitsToPlot(k))];
end
hold off
ylim([1.7 2.2]);


xlim(xlimCommon);
