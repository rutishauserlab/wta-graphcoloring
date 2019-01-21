
%%

I=1;

Ge=1;
Gi=1;
alpha1=1.2;

beta1=3;     % 
beta2=0.25;


Ge_var=1:0.2:10;


gain_GeVar = I./(Ge_var+beta1*beta2*(1/Gi) - alpha1)


Gi_var=0.5:0.1:3;

gain_GiVar = I./(Ge+beta1*beta2*(1./Gi_var) - alpha1)

maxGiVal = -beta1*beta2/(Ge-alpha1)

%
figure(20);
subplot(2,2,1);
plot(Ge_var,gain_GeVar,'x-');
xlabel('Ge, Gi=1');
ylabel('gain');

subplot(2,2,2);
plot(Gi_var,gain_GiVar, 'x-');
xlabel('Gi, Ge=1');
ylabel('gain');
title(['max Gi possible: ' num2str(maxGiVal) ]);
