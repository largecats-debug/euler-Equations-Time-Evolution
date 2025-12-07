%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%##\ .--. \##\ .-'. \##\ .--. \##\ .--. \##\ .--. \##\ .--. \##\ .`-. \##\.--%%%%
%%%%%:::::.\::::::::.\::::::::.\::::::::.\::::::::.\::::::::.\::::::::.\:::::::%%%%%
%%%%-' \##\ `--' \##\ `.-' \##\ `--' \##\ `--' \##\ `--' \##\ `-.' \##\ `--' \##%%%%
%%%%% \#\ /""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7 \#\ %%%%%
%%%%\  .--/""7__/""                                             /""7__/""7\ .--.%%%%
%%%%%:::::/""7__/"         2-D DISC - ISENTROPIC - MAIN          \"7__/""7:::::%%%%%
%%%%%:::::/""7__/"                RADIAL VARIANT                 \"7__/""7:::::%%%%%
%%%%-'  \%/""7__/""                                             /""7__/""7-' \#\%%%%
%%%%% \#\ /""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7 \#\.%%%%%
%%%%##\ .--. \##\ .-'. \##\ .--. \##\ .--. \##\ .--. \##\ .--. \##\ .`-. \##\.--%%%%
%%%%%:::::.\::::::::.\::::::::.\::::::::.\::::::::.\::::::::.\::::::::.\:::::::%%%%%
%%%%-' \##\ `--' \##\ `.-' \##\ `--' \##\ `--' \##\ `--' \##\ `-.' \##\ `--' \##%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Main iteration script to call once the user reaches the point in the main
% isentropic iteration script at which they are prompted for pertThetaIndex
% and set this value to 0


%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%Since no interaction terms can produce eigenfunction terms of any other
% theta-index, we set thetaBasisSize = 1 and prompt the user on if they
% would like to set a new rBasisSize

thetaBasisSize = 1;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
clear prompt
prompt = [newline,'The theta-index of the chosen k-mode was 0, so thetaBasisSize has been set to 1', newline, 'Set a new (or the same, if still desired) value for rBasisSize: '];
rBasisSize = input([prompt,newline]);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


%Prompt the user for the number of theta-points to use when reconstructing
% everything on the whole disk at the end of the problem

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
clear prompt
prompt = [newline,'This means that after the computations, the number of points in the discretization of [0,2pi]', newline, 'that surfaces are plotted on after the computations will be thetaPoints = ',int2str(thetaPoints), newline, 'Set a new (or the same) value for thetaPoints: '];
thetaPoints = input([prompt,newline]);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%We still define a thetaGrid so that after the computations are run, we
%   have a mesh covering the whole domain on which to recompose and plot
%   the matricies we have the values of on a single ray
thetaGridFull = linspace(0,2*pi,thetaPoints);



rPoints = 500;
rGrid = transpose(linspace(0,1,rPoints));


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


pBar = 101325;
pConstIntMatrix = pBar*ones(rPoints,1);

%Prompt user for a temperature with which to determine the baseline air
%   density to go along with the baseline air pressure pBar = 101325 Pa
for i = 1:radialEntropyLevels
    clear prompt
    prompt = ['Specify the temperature (in degrees Fahrenheit, between 50 and 80 F)', newline, 'to determine baseline air density (and thus specific volume): '];
    backgroundTemperature = input([prompt,newline]);
end

seaLevelAirCp = [1003.84,1003.87,1003.89,1003.92,1003.94,1003.97,1004.00,1004.02,1004.05,1004.08,1004.11,1004.13,1004.16,1004.19,1004.22,1004.25,1004.28,1004.30,1004.33,1004.36,1004.39,1004.42,1004.45,1004.48,1004.51,1004.55,1004.58,1004.61,1004.64,1004.67,1004.70];
seaLevelAirCv = [716.79,716.82,716.84,716.87,716.89,716.92,716.95,716.97,717.00,717.03,717.06,717.08,717.11,717.14,717.17,717.20,717.23,717.25,717.28,717.31,717.34,717.37,717.40,717.43,717.46,717.50,717.53,717.56,717.59,717.62,717.65];
seaLevelAirGamma = seaLevelAirCp./seaLevelAirCv;
seaLevelAirDensity = [1.246,1.244,1.241,1.239,1.237,1.234,1.232,1.229,1.227,1.225,1.222,1.22,1.218,1.215,1.213,1.211,1.208,1.206,1.204,1.202,1.199,1.197,1.195,1.193,1.19,1.188,1.186,1.184,1.182,1.179,1.177];

rhoBar = (1-mod(backgroundTemperature,1))*seaLevelAirDensity(floor(backgroundTemperature)-49) + ...
    mod(backgroundTemperature,1) * seaLevelAirDensity(floor(backgroundTemperature)-49 + 1);
Cp = (1-mod(backgroundTemperature,1))*seaLevelAirCp(floor(backgroundTemperature)-49) + ...
    mod(backgroundTemperature,1) * seaLevelAirCp(floor(backgroundTemperature)-49 + 1);
Cv = (1-mod(backgroundTemperature,1))*seaLevelAirCv(floor(backgroundTemperature)-49) + ...
    mod(backgroundTemperature,1) * seaLevelAirCv(floor(backgroundTemperature)-49 + 1);
gamma = Cp/Cv;
vBar = 1/rhoBar;


%Set sigmaBar = sigmaBar(pBar,s)  ( the inverse of the wave speed c )
%   by sigmaBar = sqrt(-dv/dp)|_{p=pBar} with equation of state
%   v(p,s) = vBar * [ (p/pBar)^(-1/gamma) ] * [e^(-s/Cp)]
%   where the last term can be ignored, as this is isentropic, s(x) = 0
sigmaBarConst = (1/sqrt(gamma))*sqrt(vBar)*(1/sqrt(pBar));


%Now generate the eigenfunction basis for p and v (which are the same
% because the v-basis is sigmaBar^2*pBasis, but since this is the
% isentropic case, sigmaBar^2=const does nothing after normalizing)
[pEigenBasis,pEigenFncZero,pBasisWeight,vEigenBasis,vEigenFncZero,vBasisWeight,eigenFreqMatrix,JZeroDerivLocations] = ...
    TwoD_Disk_TimeEvo_Isentropic_Basis_Radial(rBasisSize,sigmaBarConst,domainRadius,rGrid,rPoints);


pEigenFncPertTemp = zeros(rPoints,1);
pEigenFncPertTemp(:) = pEigenBasis(:,pertRIndex);
%pEigenFncPertTemp = zeros(rPoints,1);
%pEigenFncPertTemp(:,1) = pEigenBasis(pertRIndex,:);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%Set percentage of pBar chosen k-mode eigenfunction should perturb the
%   background state by:

%-*-%-*-%-*-%-*-%-*
   alpha = 0.04;  %
%-*-%-*-%-*-%-*-%-*

perturbationCoeff = (alpha/abs((min(min(pEigenFncPertTemp)))))*pBar;

%Set the time period so that the eigenfunction we are perturbing the
%   background pressure by is a solution of the linearization around pBar.
%   We consider our solution as being 2T-periodic so set T = pi/lambda_k.
timePeriod = pi/eigenFreqMatrix(pertRIndex);
timeGrid = linspace(0,timePeriod,timePoints);
timeStepSize = timeGrid(2) - timeGrid(1);

%Set initial data for u, which will just be u(0,x) = 0 -> u_k(0) = 0 if we
% want to run an iteration and find a solution
% (although we can still run the NL evolution with non zero u-data)
uIntDecompMatrixEvo = zeros(rBasisSize,1);

%The unstored variables in the NL evolution are the matrix and decomp 
%   variables for v in case thsoe are needed for something later:
%   [pNLMatrix,pNLDecompMatrix,pNLZeroDecomp,vNLMatrix,vNLDecompMatrix,vNLZeroDecomp,uNLDecompMatrix]


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if performBasicIteration == true
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    if performInitialGuessCorrection == true
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

        pNLMatrix = zeros(rPoints,timePoints,numberIterations+2);
        pNLDecompMatrix = zeros(rBasisSize,timePoints,numberIterations+2);
        pNLZeroDecomp = zeros(timePoints,numberIterations+2);
        vNLMatrix = zeros(rPoints,timePoints,numberIterations+2);
        vNLDecompMatrix = zeros(rBasisSize,timePoints,numberIterations+2);
        vNLZeroDecomp = zeros(timePoints,numberIterations+2);
        uNLDecompMatrix = zeros(rBasisSize,timePoints,numberIterations+2);
        pIntMatrixIteration = zeros(rPoints,numberIterations+2);

        pIntMatrixIteration(:,1) = pConstIntMatrix + perturbationCoeff*pEigenFncPertTemp;

        [pIntMatrixIteration(:,2),pNLDecompMatrix(:,1,2),pNLZeroDecomp(1,2)] = ...
            TwoD_Disk_TimeEvo_Isentropic_InitialCorrection_Radial(pIntMatrixIteration(:,1),pEigenBasis, ...
            pEigenFncZero,pBasisWeight,rBasisSize,rGrid,rPoints, ...
            pEigenFncPertTemp,eigenFreqMatrix,perturbationCoeff,pertRIndex, ...
            pBar,vBar,sigmaBarConst,gamma,timePeriod);

        tStart_NLEvolution = tic;
        [pNLMatrix(:,:,1),pNLDecompMatrix(:,:,1),pNLZeroDecomp(:,1),vNLMatrix(:,:,1),vNLDecompMatrix(:,:,1),vNLZeroDecomp(:,1),uNLDecompMatrix(:,:,1)] = ...
            TwoD_Disk_TimeEvo_Isentropic_NLEvolution_Radial(pIntMatrixIteration(:,1),uIntDecompMatrixEvo,rBasisSize,rGrid,rPoints, ...
            timeGrid,timePoints,timeStepSize,pEigenBasis,pEigenFncZero,pBasisWeight,vEigenBasis,vEigenFncZero,vBasisWeight,eigenFreqMatrix, ...
            sigmaBarConst,gamma,pertRIndex,pBar,vBar);
        elapsed_NLEvolution = toc(tStart_NLEvolution)

        tStart_NLEvolution = tic;
        [pNLMatrix(:,:,2),pNLDecompMatrix(:,:,2),pNLZeroDecomp(:,2),vNLMatrix(:,:,2),vNLDecompMatrix(:,:,2),vNLZeroDecomp(:,2),uNLDecompMatrix(:,:,2)] = ...
            TwoD_Disk_TimeEvo_Isentropic_NLEvolution_Radial(pIntMatrixIteration(:,2),uIntDecompMatrixEvo,rBasisSize,rGrid,rPoints, ...
            timeGrid,timePoints,timeStepSize,pEigenBasis,pEigenFncZero,pBasisWeight,vEigenBasis,vEigenFncZero,vBasisWeight,eigenFreqMatrix, ...
            sigmaBarConst,gamma,pertRIndex,pBar,vBar);
        elapsed_NLEvolution = toc(tStart_NLEvolution)

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    if performInitialGuessCorrection == false
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

        pNLMatrix = zeros(rPoints,timePoints,numberIterations+1);
        pNLDecompMatrix = zeros(rBasisSize,timePoints,numberIterations+1);
        pNLZeroDecomp = zeros(timePoints,numberIterations+1);
        vNLMatrix = zeros(rPoints,timePoints,numberIterations+1);
        vNLDecompMatrix = zeros(rBasisSize,timePoints,numberIterations+1);
        vNLZeroDecomp = zeros(timePoints,numberIterations+1);
        uNLDecompMatrix = zeros(rBasisSize,timePoints,numberIterations+1);
        pIntMatrixIteration = zeros(rPoints,numberIterations+1);

        pIntMatrixIteration(:,1) = pConstIntMatrix + perturbationCoeff*pEigenFncPertTemp;

        tStart_NLEvolution = tic;
        [pNLMatrix(:,:,1),pNLDecompMatrix(:,:,1),pNLZeroDecomp(:,1),vNLMatrix(:,:,1),vNLDecompMatrix(:,:,1),vNLZeroDecomp(:,1),uNLDecompMatrix(:,:,1)] = ...
            TwoD_Disk_TimeEvo_Isentropic_NLEvolution_Radial(pIntMatrixIteration(:,1),uIntDecompMatrixEvo,rBasisSize,rGrid,rPoints, ...
            timeGrid,timePoints,timeStepSize,pEigenBasis,pEigenFncZero,pBasisWeight,vEigenBasis,vEigenFncZero,vBasisWeight,eigenFreqMatrix, ...
            sigmaBarConst,gamma,pertRIndex,pBar,vBar);
        elapsed_NLEvolution = toc(tStart_NLEvolution)

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if performBasicIteration == false
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pNLMatrix = zeros(rPoints,timePoints);
    pNLDecompMatrix = zeros(rBasisSize,timePoints);
    pNLZeroDecomp = zeros(timePoints,1);
    uNLDecompMatrix = zeros(rBasisSize,timePoints);
    pIntMatrix = zeros(rPoints,1);

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    if performInitialGuessCorrection == true
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

        pNLMatrixCorrected = zeros(rPoints,timePoints);
        pNLDecompMatrixCorrected = zeros(rBasisSize,timePoints);
        pNLZeroDecompCorrected = zeros(timePoints,1);
        uNLDecompMatrixCorrected = zeros(rBasisSize,timePoints);
        pIntMatrix = pConstIntMatrix + perturbationCoeff*pEigenFncPertTemp;
        pIntMatrixCorrected = zeros(rPoints,1);

        [pIntMatrixCorrected(:),pNLDecompMatrixCorrected(:,1),pNLZeroDecompCorrected(1)] = ...
            TwoD_Disk_TimeEvo_Isentropic_InitialCorrection_Radial(pIntMatrix,pEigenBasis, ...
            pEigenFncZero,pBasisWeight,rBasisSize,rGrid,rPoints, ...
            pEigenFncPertTemp,eigenFreqMatrix,perturbationCoeff, ...
            pertRIndex,pBar,vBar,sigmaBarConst,gamma,timePeriod);

        tStart_NLEvolution = tic;
        [pNLMatrix,pNLDecompMatrix,pNLZeroDecomp,~,~,~,uNLDecompMatrix] = ...
            TwoD_Disk_TimeEvo_Isentropic_NLEvolution_Radial(pIntMatrix,uIntDecompMatrixEvo, ...
            rBasisSize,rGrid,rPoints,timeGrid,timePoints,timeStepSize,pEigenBasis, ...
            pEigenFncZero,pBasisWeight,vEigenBasis,vEigenFncZero,vBasisWeight, ...
            eigenFreqMatrix,sigmaBarConst,gamma,pertRIndex,pBar,vBar);
        elapsed_NLEvolution = toc(tStart_NLEvolution)

        tStart_NLEvolution = tic;
        [pNLMatrixCorrected,pNLDecompMatrixCorrected,pNLZeroDecompCorrected,~,~,~,uNLDecompMatrixCorrected] = ...
            TwoD_Disk_TimeEvo_Isentropic_NLEvolution_Radial(pIntMatrixEvoCorrection, ...
            uIntDecompMatrixEvo,rBasisSize,rGrid,rPoints,timeGrid,timePoints,timeStepSize, ...
            pEigenBasis,pEigenFncZero,pBasisWeight,vEigenBasis,vEigenFncZero,vBasisWeight, ...
            eigenFreqMatrix,sigmaBarConst,gamma,pertRIndex,pBar,vBar);
        elapsed_NLEvolution = toc(tStart_NLEvolution)

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    if performInitialGuessCorrection == false
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

        pIntMatrix = pConstIntMatrix + perturbationCoeff*pEigenFncPertTemp;

        tStart_NLEvolution = tic;
        [pNLMatrix,pNLDecompMatrix,pNLZeroDecomp,~,~,~,uNLDecompMatrix] = ...
            TwoD_Disk_TimeEvo_Isentropic_NLEvolution_Radial(pIntMatrix,uIntDecompMatrixEvo, ...
            rBasisSize,rGrid,rPoints,timeGrid,timePoints,timeStepSize,pEigenBasis, ...
            pEigenFncZero,pBasisWeight,vEigenBasis,vEigenFncZero,vBasisWeight, ...
            eigenFreqMatrix,sigmaBarConst,gamma,pertRIndex,pBar,vBar);
        elapsed_NLEvolution = toc(tStart_NLEvolution)

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 






%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
if performBasicIteration == true
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    %Our iteration by which we refine a solution to the Euler equations is:
    %   P^-1*f(y^(k)_0) = P^-1*Df(y^(k)_0)[r^(k)]
    %   where f(y^(k)_0) is a column vector representation of u^(k)(T,x),
    %   Df(y^(k)_0)[r^k] is a matrix mapping P^(k)(0,x) to u^(k)(T,x),
    %   and P is a suitable preconditioning matrix

    %NOTE 9/18/25:
    %   In the case of the 2D disc (at least for isentropic?) the only resonant
    %   eigenfunction to any given eigenfunction is that with the same theta-index
    %   and r-index but opposite trig function sin/cos. Would this be considered
    %   resonance (or some kind of it)? What is their iteraction term?
    %What we can say for now is this:
    %   this means we never need to consider how to separately handle resonant
    %   eigenfunctions as every eigenfunction we are using will be nonresonant
    %   to our chosen starting k-mode (at least in isentropic, research needed)

    smallDivisorsBasicNoZero = zeros(rBasisSize,1);

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    for m = 1:rBasisSize
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

        smallDivisorsBasicNoZero(m) = sin(eigenFreqMatrix(m)*timePeriod);

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    %For every j-mode for j ~= k (our chosen fixed mode), we are able to solve
    %   the j-mode equation with the j-mode.
    %For the k-mode equation however, we do not have a k-mode variable to solve
    %   it with. Instead this becomes a bifurcation problem, and we use the
    %   0-mode to solve the k-mode equation.
    %   From Duhamel, we find that the k-mode effect due to adding
    %   1*phi_(00) to our initial guess pBar + c*phi_k is
    %   f(pBar + c*phi_k * phi_0) = f(pBar) + Df(pBar)[c*phi_k + phi_0]
    %   + c^2/2*D^2f(pBar)[phi_k,phi_k] + 1/2 D^2f(pBar)[phi_0,phi_0]
    %   + c D^2f(pBar)[phi_0,phi_k]

    %Only have alpha*D^2f(pBar)[phi_0,phi_k] which was found to be
    %   alpha*D^2f(pBar)[phi_0,phi_k] = -alpha*(pi/2)*|phi_0|*v''(pBar)
    %   = -(pi/2)*[(1+gamma)/(gamma*pBar)]*sigmaBar^2*alpha
    %   Here we are assuming isentropic so phi_0 is a scalar, so |phi_0| is just
    %   meant to mean the value of the zero eigenfunction at any point
    
    %NOTE APRIL 25: pBar vs pBar in calculation of these second derivative
    %   constants/small divisors. Why are we linearizing around pBar instead of
    %   the pBar we are actually using? Seems like we need to use pBar though
    %   not important until we consider non-isentropic since in isentropic pBar =
    %   pBar though
    %We store delta_(0,0) as just <D^2f(pBar)[phi_(00),phi_k],psi_k>
    %   as this is the same value as <D^2f(pBar)[phi_(00),phi_j],psi_j> for phi_j
    %   that are resonant with phi_k (lambda_k = lambda_j), so the alpha will be
    %   added on where needed

    %smallDivisorZeroMode = (-1) * ...
        %(((1+gamma)*vBar*eigenFreqMatrix(2,pertRIndex)*timePeriod)/(2*(gamma^2)*(pBar^2))) * ...
        %pEigenFncZero(1,1);

    smallDivisorZeroMode = (-1) * ...
        ( ( (1+gamma)*vBar*eigenFreqMatrix(pertRIndex)*timePeriod ) / ( 2*(gamma^2)*(pBar) ) ) * ...
        ( 1/sigmaBarConst );


    %NOTE APRIL 25: Change this to loop over arbitrary nubmer of iterations
    %   after testing single basic iteration. Otherwise need to set up function
    %   for linearization evolution around nonlienar state

    zeroCorrectionBasic = zeros(1,numberIterations);
    nonResCorrectionBasic = zeros(rBasisSize,numberIterations);
    smallDivisorsPBarPlusZ = zeros(rBasisSize,numberIterations);
    nonResCorrectionPBarPlusZ = zeros(rBasisSize,numberIterations);
    pEigenFncIterationLoop = zeros(rPoints,1);



    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~performInitialGuessCorrection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    
    
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        for i = 1:numberIterations
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

            tStart_SolutionIteration = tic;

            pIntMatrixIteration(:,i+1) = pIntMatrixIteration(:,i);
    
            %We start by finding the corrections for the psi_k term, 
            %   which will be solved with the zero eigenfunction (where phi_k is our
            %   perturbed eigenfunction)
            %The phi_k equation for which we pick the z in +z*phi_(00) to solve is
            %   z*alpha * D^2f(pBar)[phi_(00),phi_k] = -<f(p(0,x)),psi_k>
    
            zeroCorrectionBasic(i) = (-1)*uNLDecompMatrix(pertRIndex,timePoints,i) /...
                (perturbationCoeff * smallDivisorZeroMode);


            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
            for m = 1:rBasisSize
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
       
                pEigenFncIterationLoop(:) = pEigenBasis(:,m);

                %For the nonresonant modes, the basic iteration will just
                %   attempt to correct the psi_j terms by choosing beta_j such that
                %   beta_j*<Df(pBar)[phi_j],psi_j> = -<f(p(0,x)),psi_j>
                %   where p(0,x) is the previous guess for the initial data that is
                %   being improved in this iteration

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if m ~= pertRIndex
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    nonResCorrectionBasic(m,i)=(-1)*uNLDecompMatrix(m,timePoints,i)/ ...
                        smallDivisorsBasicNoZero(m);

                    pIntMatrixIteration(:,i+1) = pIntMatrixIteration(:,i+1) + ...
                        nonResCorrectionBasic(m,i) * pEigenFncIterationLoop(:);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                nonResCorrectionBasic(m,i)=0;

                pIntMatrixIteration(:,i+1) = pIntMatrixIteration(:,i+1) + (zeroCorrectionBasic(i)*...
                    pEigenFncZero(:));

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
            end
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


            [pNLMatrix(:,:,i+1),pNLDecompMatrix(:,:,i+1),pNLZeroDecomp(:,i+1),vNLMatrix(:,:,i+1),vNLDecompMatrix(:,:,i+1),vNLZeroDecomp(:,i+1),uNLDecompMatrix(:,:,i+1)] = ...
                TwoD_Disk_TimeEvo_Isentropic_NLEvolution_Radial(pIntMatrixIteration(:,i+1),uIntDecompMatrixEvo,rBasisSize,rGrid,rPoints, ...
                timeGrid,timePoints,timeStepSize,pEigenBasis,pEigenFncZero,pBasisWeight,vEigenBasis,vEigenFncZero,vBasisWeight, ...
                eigenFreqMatrix,sigmaBarConst,gamma,pertRIndex,pBar,vBar);


            elapsed_SolutionIteration = toc(tStart_SolutionIteration)

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        end
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 




    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 





    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if performInitialGuessCorrection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        for i = 2:numberIterations+1
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

            tStart_SolutionIteration = tic;

            pIntMatrixIteration(:,i+1) = pIntMatrixIteration(:,i);

            %We start by finding the corrections for the psi_k term, 
            %   which will be solved with the zero eigenfunction (where phi_k is our
            %   perturbed eigenfunction)
            %The phi_k equation for which we pick the z in +z*phi_(00) to solve is
            %   z*alpha * D^2f(pBar)[phi_(00),phi_k] = -<f(p(0,x)),psi_k>

            zeroCorrectionBasic(i) = (-1)*uNLDecompMatrix(pertRIndex,timePoints,i) /...
                (perturbationCoeff * smallDivisorZeroMode);


                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
                for m = 1:rBasisSize
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
       
                    pEigenFncIterationLoop(:) = pEigenBasis(:,m);

                    %For the nonresonant modes, the basic iteration will just
                    %   attempt to correct the psi_j terms by choosing beta_j such that
                    %   beta_j*<Df(pBar)[phi_j],psi_j> = -<f(p(0,x)),psi_j>
                    %   where p(0,x) is the previous guess for the initial data that is
                    %   being improved in this iteration

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if m ~= pertRIndex
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        nonResCorrectionBasic(m,i)=(-1)*uNLDecompMatrix(m,timePoints,i)/ ...
                            smallDivisorsBasicNoZero(m);

                        pIntMatrixIteration(:,i+1) = pIntMatrixIteration(:,i+1) + ...
                            nonResCorrectionBasic(m,i) * pEigenFncIterationLoop(:);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    else
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        nonResCorrectionBasic(m,i)=0;

                        pIntMatrixIteration(:,i+1) = pIntMatrixIteration(:,i+1) + (zeroCorrectionBasic(i)*...
                            pEigenFncZero(:,1));

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
                end
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


            [pNLMatrix(:,:,i+1),pNLDecompMatrix(:,:,i+1),pNLZeroDecomp(:,i+1),vNLMatrix(:,:,i+1),vNLDecompMatrix(:,:,i+1),vNLZeroDecomp(:,i+1),uNLDecompMatrix(:,:,i+1)] = ...
            TwoD_Disk_TimeEvo_Isentropic_NLEvolution_Radial(pIntMatrixIteration(:,i+1),uIntDecompMatrixEvo,rBasisSize,rGrid,rPoints, ...
            timeGrid,timePoints,timeStepSize,pEigenBasis,pEigenFncZero,pBasisWeight,vEigenBasis,vEigenFncZero,vBasisWeight,eigenFreqMatrix, ...
            sigmaBarConst,gamma,pertRIndex,pBar,vBar);

            elapsed_SolutionIteration = toc(tStart_SolutionIteration)

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        end
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%for n = 1:thetaPoints

%    p231Ray(n,:) = pEigenBasis(3,3,1,n,:);
%    p241Ray(n,:) = pEigenBasis(3,4,1,n,:);
%    p131Ray(n,:) = pEigenBasis(2,3,1,n,:);
%    p232Ray(n,:) = pEigenBasis(3,3,2,n,:);

%    pFncDiffRRayInt(n) = trapz(rGrid,(p231Ray(n,:).*p241Ray(n,:)).*rGrid);
%    pFncDiffThetaRayInt(n) = trapz(rGrid,(p231Ray(n,:).*p131Ray(n,:)).*rGrid);
%    pFncDiffTrigRayInt(n) = trapz(rGrid,(p231Ray(n,:).*p232Ray(n,:)).*rGrid);
%    pFncSameRayInt(n) = trapz(rGrid,(p231Ray(n,:).*p231Ray(n,:)).*rGrid);

%end

%pFncDiffRDoubleInt = trapz(thetaGrid,pFncDiffRRayInt);
%pFncDiffRInnerProd = trapz(thetaGrid,pFncDiffRRayInt);
%pFncDiffTrigInnerProd = trapz(thetaGrid,pFncDiffTrigRayInt);
%pFncSameInnerProd = trapz(thetaGrid,pFncSameRayInt);


%for n = 1:thetaPoints
%J11CosCosVectorRayInt(n) = (1/pi)*cos(thetaGrid(n))*cos(thetaGrid(n))*trapz(rGrid,(J11VectorTemp.*J11VectorTemp).*rGrid);
%J11SinSinVectorRayInt(n) = (1/pi)*sin(thetaGrid(n))*sin(thetaGrid(n))*trapz(rGrid,(J11VectorTemp.*J11VectorTemp).*rGrid);
%J11CossSinVectorRayInt(n) = (1/pi)*sin(thetaGrid(n))*cos(thetaGrid(n))*trapz(rGrid,(J11VectorTemp.*J11VectorTemp).*rGrid);
%end


%J11CosCosInnerProd = trapz(thetaGrid,J11CosCosVectorRayInt)
%J11SinSinInnerProd = trapz(thetaGrid,J11SinSinVectorRayInt)
%J11CosSinInnerProd = trapz(thetaGrid,J11CosSinVectorRayInt)


%clear pNLMatrix pNLDecompMatrix vNLMatrix vNLDecompMatrix



%for n = 1:thetaBasisSize
%    for m = 1:rBasisSize
%        JScaledVectorLoop(:) = JScaledValues(n,m,:);

%        JScaledValuesNorm(n,m) = sqrt(trapz(rGrid,(JScaledVectorLoop.*JScaledVectorLoop).*rGrid));
%        JScaledValuesNormalized(n,m,:) = JScaledValues(n,m,:)/JScaledValuesNorm(n,m);
%    end
%end

%JInnerProdChoice(:) = JScaledValuesNormalized(2,3,:);
%JInnerProdCheck = trapz(rGrid,(JInnerProdChoice.*JInnerProdChoice).*rGrid)






