%
%calculate the gain as a function of network parameters
%
%
function gain = calcGain(alpha1, beta1, beta2, Ge, Gi)

gain = 1./(Ge + (beta1*beta2)/Gi - alpha1);