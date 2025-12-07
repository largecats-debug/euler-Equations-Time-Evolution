%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
function [pNLMatrix,pNLDecompMatrix,pNLZeroDecomp,vNLMatrix,vNLDecompMatrix,vNLZeroDecomp,uNLDecompMatrix] = TwoD_Disk_TimeEvo_Isentropic_NLEvolution_Radial( ...
    pIntMatrix,uIntDecompMatrix,rBasisSize,rGrid,rPoints,timeGrid,timePoints,timeStepSize,pEigenBasis,pEigenFncZero,...
    pBasisWeight,vEigenBasis,vEigenFncZero,vBasisWeight,eigenFreqMatrix,sigmaBarConst,gamma,pertRIndex,pBar,vBar)


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


%NOTE 9/18/25: 
% This code is being written with a direct inner product decomposition in
% mind. More research is still needed on the possibility of some kind of
% FBT (fast Bessel transform)


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialization of matricies for values of NL variables p and v
% No such matrix is needed for the values of u as they are never actually
% used, only the decomposition coefficients u_nm. We do need the actual
% values of p and v, as the only way to go between them is:
%  1. recompose v by v = sum(v_nm * phi_nm)
%  2. apply inverse equation of state p = pBar * [(v/vBar)^gamma]
%  3. decompose p for the p_nm's
pNLMatrix = zeros(rPoints,timePoints);
vNLMatrix = zeros(rPoints,timePoints);

%Can make the u decomposition matrix with or without the 0-eigenfunction
%component, as psi_0 = 0 anyway
pNLDecompMatrix = zeros(rBasisSize,timePoints);
vNLDecompMatrix = zeros(rBasisSize,timePoints);
uNLDecompMatrix = zeros(rBasisSize,timePoints);


pNLZeroDecomp = zeros(timePoints,1);
vNLZeroDecomp = zeros(timePoints,1);

matrixTempNL = zeros(rPoints,1);
matrixTempDecomp = zeros(rBasisSize,1);

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%Initialize pNLMatrix at t=0 (index 1)
pNLMatrix(:,1) = pIntMatrix;

%Decompose pNLMatrix at t=0 to initialize pNLDecompMatrix
matrixTempNL(:) = pNLMatrix(:,1);
[pNLDecompMatrix(:,1),pNLZeroDecomp(1)] = TwoD_Disk_TimeEvo_Isentropic_InnerDecomp_Radial( ...
    matrixTempNL,pEigenBasis,pEigenFncZero,pBasisWeight,rBasisSize,rGrid,rPoints,pertRIndex);

%Calculate vNLMatrix and decompose for vNLDecompMatrix at t=0
vNLMatrix(:,1) = vBar * (matrixTempNL.^((-1)*(1/gamma))) * (pBar^(1/gamma));
matrixTempNL(:) = vNLMatrix(:,1);
[vNLDecompMatrix(:,1),vNLZeroDecomp(1)] = TwoD_Disk_TimeEvo_Isentropic_InnerDecomp_Radial( ...
    matrixTempNL,vEigenBasis,vEigenFncZero,vBasisWeight,rBasisSize,rGrid,rPoints,pertRIndex);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NOTE 9/18/25:
% To use a leapfrog type scheme, we need a second piece of initial data to
% start the process. Need to redesign an initial data extension method by
% way of RK4. Currently, this code just uses an/some Euler step(s) to 
% initialize the data at the second time step, as a rough approximation
% easy to get working.

numberEulerSteps = 4;
subStepSize = timeStepSize/numberEulerSteps;


pNLSubSteps = zeros(rPoints,numberEulerSteps);
vNLSubSteps = zeros(rPoints,numberEulerSteps);

pNLDecompSubSteps = zeros(rBasisSize,numberEulerSteps);
vNLDecompSubSteps = zeros(rBasisSize,numberEulerSteps);
uNLDecompSubSteps = zeros(rBasisSize,numberEulerSteps);

pNLZeroDecompSubSteps = zeros(numberEulerSteps,1);
vNLZeroDecompSubSteps = zeros(numberEulerSteps,1);


pNLSubSteps(:,1) = pNLMatrix(:,1);
vNLSubSteps(:,1) = vNLMatrix(:,1);

pNLDecompSubSteps(:,1) = pNLDecompMatrix(:,1);
vNLDecompSubSteps(:,1) = vNLDecompMatrix(:,1);
uNLDecompSubSteps(:,1) = uNLDecompMatrix(:,1);

pNLZeroDecompSubSteps(1) = pNLZeroDecomp(1);
vNLZeroDecompSubSteps(1) = vNLZeroDecomp(1);


for e1 = 1:numberEulerSteps
    for m = 1:rBasisSize
        uNLDecompSubSteps(m,e1+1) = uNLDecompSubSteps(m,e1) + timeStepSize*eigenFreqMatrix(m)*pNLDecompSubSteps(m,e1);
        vNLDecompSubSteps(m,e1+1) = vNLDecompSubSteps(m,e1) + timeStepSize*eigenFreqMatrix(m)*uNLDecompSubSteps(m,e1);
    end

    vNLZeroDecompSubSteps(e1+1) = vNLZeroDecompSubSteps(e1);
    matrixTempDecomp(:) = vNLDecompSubSteps(:,e1+1);

    [pNLSubSteps(:,e1+1),pNLDecompSubSteps(:,e1+1),pNLZeroDecompSubSteps(e1+1),vNLSubSteps(:,e1+1)] = TwoD_Disk_TimeEvo_Isentropic_VtoP_Radial( ...
    vNLDecompSubSteps(:,e1+1),vNLZeroDecompSubSteps(e1+1),pEigenBasis,pEigenFncZero,pBasisWeight,vEigenBasis,vEigenFncZero, ...
    rBasisSize,rGrid,rPoints,gamma,pertRIndex,pBar,vBar);
end

pNLMatrix(:,2) = pNLSubSteps(:,numberEulerSteps+1);
vNLMatrix(:,2) = vNLSubSteps(:,numberEulerSteps+1);
pNLDecompMatrix(:,2) = pNLDecompSubSteps(:,numberEulerSteps+1);
vNLDecompMatrix(:,2) = vNLDecompSubSteps(:,numberEulerSteps+1);
uNLDecompMatrix(:,2) = uNLDecompSubSteps(:,numberEulerSteps+1);


%vNLDecompMatrix(:,2) = vNLDecompMatrix(:,1);
%pNLDecompMatrix(:,2) = pNLDecompMatrix(:,1);
%pNLZeroDecomp(2) = pNLZeroDecomp(1);
%vNLZeroDecomp(2) = vNLZeroDecomp(1);
%pNLMatrix(2,:) = pNLMatrix(1,:);
%vNLMatrix(2,:) = vNLMatrix(1,:);

%for m = 1:rBasisSize
    %uNLDecompMatrix(m,2) = uNLDecompMatrix(m,1) + ...
        %timeStepSize*eigenFreqMatrix(m)*pNLDecompMatrix(m,1);
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 





% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 3:timePoints
    if j > 3
        elapsed_NLEvolution = toc(tStart_EvolutionIteration)
    end
    vNLZeroDecomp(j) = vNLZeroDecomp(j-1);
    tStart_EvolutionIteration = tic;
    %perform leapfrog step on the decomposition components
    %from the weak* solution of the NL system this iterates a finite
    %difference/spectral version of
    % d_t u_j - lambda_j * p_j = 0
    % d_t v_j - lambda_j * u_j = 0


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        for m = 1:rBasisSize

            uNLDecompMatrix(m,j) = uNLDecompMatrix(m,j-2) + ...
                2*timeStepSize*eigenFreqMatrix(m)*pNLDecompMatrix(m,j-1);
            vNLDecompMatrix(m,j) = vNLDecompMatrix(m,j-2) + ...
                2*timeStepSize*eigenFreqMatrix(m)*uNLDecompMatrix(m,j-1);

        end
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    matrixTempDecomp(:) = vNLDecompMatrix(:,j);

    [pNLMatrix(:,j),pNLDecompMatrix(:,j),pNLZeroDecomp(j),vNLMatrix(:,j)] = TwoD_Disk_TimeEvo_Isentropic_VtoP_Radial( ...
    matrixTempDecomp,vNLZeroDecomp(j),pEigenBasis,pEigenFncZero,pBasisWeight,vEigenBasis,vEigenFncZero, ...
    rBasisSize,rGrid,rPoints,gamma,pertRIndex,pBar,vBar);

end

elapsed_NLEvolution = toc(tStart_EvolutionIteration)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 




end
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %