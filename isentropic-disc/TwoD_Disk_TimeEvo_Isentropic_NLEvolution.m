%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
function [pNLMatrix,pNLDecompMatrix,pNLZeroDecomp,vNLMatrix,vNLDecompMatrix,vNLZeroDecomp,uNLDecompMatrix] = TwoD_Disk_TimeEvo_Isentropic_NLEvolution( ...
    pIntMatrix,uIntDecompMatrix,rBasisSize,rGrid,rPoints,thetaBasisSize,thetaGrid,thetaPoints,timeGrid,timePoints,timeStepSize,pEigenBasis,pEigenFncZero,...
    pBasisWeight,vEigenBasis,vEigenFncZero,vBasisWeight,eigenFreqMatrix,sigmaBarConst,gamma,pertRIndex,pertThetaIndex,pBar,vBar)


%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%
%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%
%%\:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._\%%
%%\:._.:%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%~*:._\%%
%%\:._.:%_/▔\_/▔\_/▔\_/▔\_/▔\_/▔\_/▔\_/▔\_/▔\_/▔\_/▔\_ %~*:._\%%
%%\:._.:%|    NL TIME EVOLUTION - DISC - ISENTROPIC - MAIN     |%~*:._\%%
%%\:._.:%▔\_/▔\_/▔\_/▔\_/▔\_/▔\_/▔\_/▔\_/▔\_/▔\_/▔\_/▔%~*:._\%%
%%\:._.:%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%~*:._\%%
%%\:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._\%%
%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%
%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%


%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %



%NOTE 9/18/25: 
% This code is being written with a direct inner product decomposition in
% mind. More research is still needed on the possibility of some kind of
% FBT (fast Bessel transform)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%Initialization of matricies for values of NL variables p and v
% No such matrix is needed for the values of u as they are never actually
% used, only the decomposition coefficients u_nm. We do need the actual
% values of p and v, as the only way to go between them is:
%  1. recompose v by v = sum(v_nm * phi_nm)
%  2. apply inverse equation of state p = pBar * [(v/vBar)^gamma]
%  3. decompose p for the p_nm's
pNLMatrix = zeros(thetaPoints,rPoints,timePoints);
vNLMatrix = zeros(thetaPoints,rPoints,timePoints);

%Can make the u decomposition matrix with or without the 0-eigenfunction
%component, as psi_0 = 0 anyway
pNLDecompMatrix = zeros(thetaBasisSize,rBasisSize,timePoints);
vNLDecompMatrix = zeros(thetaBasisSize,rBasisSize,timePoints);
uNLDecompMatrix = zeros(thetaBasisSize,rBasisSize,timePoints);

pNLZeroDecomp = zeros(timePoints,1);
vNLZeroDecomp = zeros(timePoints,1);

matrixTempNL(:,:) = zeros(thetaPoints,rPoints);
matrixTempDecomp(:,:) = zeros(thetaBasisSize,rBasisSize);


%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%  %  %  %  %  % INITIALIZE P AND ITS DECOMPOSITION MATRIX  %  %  %  %  %  
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

%Initialize pNLMatrix at t=0 with provided matrix
pNLMatrix(:,:,1) = pIntMatrix(:,:);
matrixTempNL(:,:) = pNLMatrix(:,:,1);
[pNLDecompMatrix(:,:,1),pNLZeroDecomp(1)] = TwoD_Disk_TimeEvo_Isentropic_ScalarDecomp( ...
    matrixTempNL,pEigenBasis,pEigenFncZero,pBasisWeight,rBasisSize,thetaBasisSize, ...
    rGrid,rPoints,thetaGrid,thetaPoints,pertRIndex);


%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%  %  %  %  %  % INITIALIZE V AND ITS DECOMPOSITION MATRIX  %  %  %  %  %  
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

%Initialize vNLMatrix with equation of state for a gamma law gas
vNLMatrix(:,:,1) = vBar * (matrixTempNL.^((-1)*(1/gamma))) * (pBar^(1/gamma));
matrixTempNL(:,:) = vNLMatrix(:,:,1);
[vNLDecompMatrix(:,:,1),vNLZeroDecomp(1)] = TwoD_Disk_TimeEvo_Isentropic_ScalarDecomp( ...
    matrixTempNL,vEigenBasis,vEigenFncZero,vBasisWeight,rBasisSize,thetaBasisSize, ...
    rGrid,rPoints,thetaGrid,thetaPoints,pertRIndex);


%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%  %  %  %  %  %  % INITIALIZE U DECOMPOSITION MATRIX %  %  %  %  %  %  %  
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

uNLDecompMatrix(:,:,1) = uIntDecompMatrix;


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%NOTE 9/18/25:
% To use a leapfrog type scheme, we need a second piece of initial data to
% start the iteration. Need to redesign an initial data extension method by
% way of RK4. Currently, this code just assumes data is the same at second
% time step, as a rough approximation easy to get working.

%NOTE 9/18/25:
% To use a leapfrog type scheme, we need a second piece of initial data to
% start the process. Need to redesign an initial data extension method by
% way of RK4. Currently, this code just uses an/some Euler step(s) to 
% initialize the data at the second time step, as a rough approximation
% easy to get working.

numberEulerSteps = 5;
subStepSize = timeStepSize/(numberEulerSteps-1);

pNLSubSteps = zeros(thetaPoints,rPoints,numberEulerSteps);
vNLSubSteps = zeros(thetaPoints,rPoints,numberEulerSteps);

pNLDecompSubSteps = zeros(thetaBasisSize,rBasisSize,numberEulerSteps);
vNLDecompSubSteps = zeros(thetaBasisSize,rBasisSize,numberEulerSteps);
uNLDecompSubSteps = zeros(thetaBasisSize,rBasisSize,numberEulerSteps);

pNLZeroDecompSubSteps = zeros(numberEulerSteps,1);
vNLZeroDecompSubSteps = zeros(numberEulerSteps,1);

pNLSubSteps(:,:,1) = pNLMatrix(:,:,1);
vNLSubSteps(:,:,1) = vNLMatrix(:,:,1);
pNLDecompSubSteps(:,:,1) = pNLDecompMatrix(:,:,1);
vNLDecompSubSteps(:,:,1) = vNLDecompMatrix(:,:,1);
uNLDecompSubSteps(:,:,1) = uNLDecompMatrix(:,:,1);
pNLZeroDecompSubSteps(1) = pNLZeroDecomp(1);
vNLZeroDecompSubSteps(1) = vNLZeroDecomp(1);

pNLZeroDecompSubSteps = pNLZeroDecomp(1)*ones(numberEulerSteps,1);
vNLZeroDecompSubSteps = vNLZeroDecomp(1)*ones(numberEulerSteps,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for e1 = 2:numberEulerSteps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    for n = 1:thetaBasisSize
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for m = 1:rBasisSize
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            uNLDecompSubSteps(n,m,e1) = uNLDecompSubSteps(n,m,e1-1) + subStepSize*eigenFreqMatrix(n,m)*pNLDecompSubSteps(n,m,e1-1);
            vNLDecompSubSteps(n,m,e1) = vNLDecompSubSteps(n,m,e1-1) + subStepSize*eigenFreqMatrix(n,m)*uNLDecompSubSteps(n,m,e1-1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    matrixTempDecomp(:,:) = vNLDecompSubSteps(:,:,e1);

    [pNLSubSteps(:,:,e1),pNLDecompSubSteps(:,:,e1),pNLZeroDecompSubSteps(e1),vNLSubSteps(:,:,e1)] = TwoD_Disk_TimeEvo_Isentropic_VtoP( ...
    matrixTempDecomp,vNLZeroDecompSubSteps(e1),pEigenBasis,pEigenFncZero,pBasisWeight,vEigenBasis,vEigenFncZero, ...
    rBasisSize,thetaBasisSize,rGrid,rPoints,thetaGrid,thetaPoints,gamma,pertRIndex,pertThetaIndex,pBar,vBar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pNLMatrix(:,:,2) = pNLSubSteps(:,:,numberEulerSteps);
vNLMatrix(:,:,2) = vNLSubSteps(:,:,numberEulerSteps);
pNLDecompMatrix(:,:,2) = pNLDecompSubSteps(:,:,numberEulerSteps);
vNLDecompMatrix(:,:,2) = vNLDecompSubSteps(:,:,numberEulerSteps);
uNLDecompMatrix(:,:,2) = uNLDecompSubSteps(:,:,numberEulerSteps);
pNLZeroDecomp(2) = pNLZeroDecompSubSteps(numberEulerSteps);
vNLZeroDecomp(2) = vNLZeroDecompSubSteps(numberEulerSteps);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 3:timePoints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    if j > 3
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

        elapsed_NLEvolution = toc(tStart_EvolutionIteration)

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    vNLZeroDecomp(j) = vNLZeroDecomp(j-1);
    tStart_EvolutionIteration = tic;


    %perform leapfrog step on the decomposition components
    %from the weak* solution of the NL system this iterates a finite
    %difference/spectral version of
    % d_t u_j - lambda_j * p_j = 0
    % d_t v_j - lambda_j * u_j = 0


    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    for n = 1:thetaBasisSize
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for m = 1:rBasisSize
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            uNLDecompMatrix(n,m,j) = uNLDecompMatrix(n,m,j-2) + ...
                2*timeStepSize*eigenFreqMatrix(n,m)*pNLDecompMatrix(n,m,j-1);
            vNLDecompMatrix(n,m,j) = vNLDecompMatrix(n,m,j-2) + ...
                2*timeStepSize*eigenFreqMatrix(n,m)*uNLDecompMatrix(n,m,j-1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    matrixTempDecomp(:,:) = vNLDecompMatrix(:,:,j);

    [pNLMatrix(:,:,j),pNLDecompMatrix(:,:,j),pNLZeroDecomp(j),vNLMatrix(:,:,j)] = TwoD_Disk_TimeEvo_Isentropic_VtoP( ...
    matrixTempDecomp,vNLZeroDecomp(j),pEigenBasis,pEigenFncZero,pBasisWeight,vEigenBasis,vEigenFncZero, ...
    rBasisSize,thetaBasisSize,rGrid,rPoints,thetaGrid,thetaPoints,gamma,pertRIndex,pertThetaIndex,pBar,vBar);

end

elapsed_NLEvolution = toc(tStart_EvolutionIteration)


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




end
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %