
% mainWTA_graphColoring_assembleParams

paramsAll=[];

paramsAll = addFieldsToStruct( paramsAll, graphType, plotModeDynamic, plotModeNetwork, plotModeResult, alpha1, alpha2, beta1, beta2, beta1_DB, beta2_DB, beta1P_DB, beta2P_DB, gamma );

paramsAll = addFieldsToStruct( paramsAll,T_DB, Tdefault, Tinhib, dt, NrOfColors, stepSizeStates,till, noiseUpdateFreq, noiseStd, noiseMean);

paramsAll = addFieldsToStruct( paramsAll, enableInteract, enableNoiseInp, modeNormalizeWeights, modeWeightNoise, normConnsEnabled, normFact, fixedSeed, seedValTotal, massRun, density, NrOfColors, nrNodes);


paramsAll = addFieldsToStruct( paramsAll,biasAmp, keepHistory, terminateIfNoError, simNrRand, biasInp_forSudoku, Nnodes, Maps);

paramsAll = addFieldsToStruct( paramsAll, enableDendriticCompartment, economicMode);

if exist('dendSlope')
    paramsAll = addFieldsToStruct( paramsAll, dendSlope ); 
end