%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
function [pIntMatrixEvoCorrection,pIntCorrectedDecomp,pIntCorrectedDecompZero] = TwoD_Disk_TimeEvo_Isentropic_InitialCorrection(pIntMatrixEvo, ...
    pEigenBasis,pEigenFncZero,pBasisWeight,rBasisSize,thetaBasisSize,rGrid,rPoints,thetaGrid, ...
    thetaPoints,pEigenFncPertTemp,eigenFreqMatrix,perturbationCoeff,pertThetaIndex,pertRIndex,pertTrigIndex,pBar,vBar,sigmaBarConst,gamma,timePeriod)


%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%
%%%\:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:%%\%
%%%\%._.:*~\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%~*:._.%%%\
%%%\%._.:*~\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%~*:._.%%%\
%%%\:._.:*~*%\%\    2-D TIME EVOLUTION - DISC - ISENTROPIC    %\%*~*:._.:%%\%
%%%\:._.:*~*%\%\        INITIAL O(alpha)^2 CORRECTION         %\%*~*:._.:%%\%
%%%\:._.:*~*%\%\        SELECTIVE EIGENFUNCTION BASES         %\%*~*:._.:%%\%
%%%\%._.:*~\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%~*:._.%%%\
%%%\%._.:*~\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%~*:._.%%%\
%%%\:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:%%\%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%


%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %


%this function is to improve the initial guess of p(0,x) = pBar +
% alpha*phi_(k1,k2,l) from having error O(alpha^2) to O(alpha^3)

%For a chosen k-mode to be perturbed to a solution, we have some
% corresponding indicies thetaPert=j and rPert. The theta-indicies still
% range from 0 to N-1, but since we use "selective" eigenfunction bases,
% these indicies 0, 1, ..., N-1 actually correspond to:
% cos(0*j*theta), cos(1*j*theta), ..., cos((N-1)*j*theta)

%The leading order term in f(pBar + alpha*phi_k) that we are trying to
% cancel out to improve the starting solution is:
% (1/2)alpha^2*D^2f(pBar)[phi_k,phi_k]

%The interaction of the k-mode with itself has a term with
% cos(j*theta)*cos(j*theta) = (1/2)*[ cos(0*j*theta) - cos(2*j*theta) ]
% In the selective basis, these two trig terms making up the interaction
% term correspond to thetaIndex = 1 (0-mode) and 3 (2j-mode)

%The full interaction term (1/2)alpha^2*D^2f(pBar)[phi_k,phi_k]
% includes the square of the Bessel function in the k-mode. This must be
% decomposed back into our basis to find the corrections.

%This initial correct depends on the calculation:
% <D^2f(pBar)[phi_k,phi_k],psi_j> (for the appropriate j's)
% = sigmaBar^2*lambda_k*((1+gamma)/(2*gamma*pBar))*
%   * [-1/(lambda_j + 2*lambda_k) + 1/(lambda_j - 2*lambda_k)]
%   * sin(lambda_j*T) * <phi_k^2,phi_j>

%The effect of adding a_j phi_j to the initial guess looks like
% <f(pBar + alpha*phi_k + a_j*phi_j),psi_j>
%   = (...) + a_j*delta_j + (1/2)*alpha^2*<D^2f(pBar)[phi_k,phi_k],psi_j>
% where delta_j is the small divisor sin(lambda_j*T)



%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


% 1. calculate the decomposition of phi_k^2
% 2. compute <D^2f(pBar)[phi_k,phi_k],psi_j> for each psi_j corresponding
%    to thetaIndex = 1 (0-mode) or 3 (2j-mode)
% 3. set a_j = -(1/2)*alpha^2*<D^2f(pBar)[phi_k,phi_k],psi_j>/delta_j
% -> a_j = -sigmaBar^2*lambda_k*((1+gamma)/(4*gamma*pBar))*
%     * [-1/(lambda_j + 2*lambda_k) + 1/(lambda_j - 2*lambda_k)]
%     * <phi_k^2,phi_j>

pEigenFncPertSquared = pEigenFncPertTemp.^2;

%unstored value is zero eigenfunction decomposition
%[pEigenFncPertSquaredDecomp,~] = TwoD_Disk_TimeEvo_IsentropicScalarDecomposition( ...
%    pEigenFncPertSquared,pEigenBasis,pEigenFncZero,rBasisSize,thetaBasisSize, ...
%    rGrid,rPoints,thetaGrid,thetaPoints,pertRIndex,pertThetaIndex);
[pEigenFncPertSquaredDecomp,pEigenFncPertSquaredZeroDecomp] = TwoD_Disk_TimeEvo_Isentropic_ScalarDecomp( ...
    pEigenFncPertSquared,pEigenBasis,pEigenFncZero,pBasisWeight,rBasisSize,thetaBasisSize, ...
    rGrid,rPoints,thetaGrid,thetaPoints,pertRIndex);


pIntMatrixEvoCorrection = pIntMatrixEvo;
pEigenFncCorrectionTemp = zeros(thetaPoints,rPoints);
pIntCorrectedDecomp = zeros(thetaBasisSize,rBasisSize);
pIntCorrectedDecomp(2,pertRIndex) = perturbationCoeff;
pIntCorrectedDecompZero = pBar/pEigenFncZero(1,1);


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
for m = 1:rBasisSize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pEigenFncCorrectionTemp(:,:) = pEigenBasis(1,m,:,:);

    correctionCoeffTemp = pEigenFncPertSquaredDecomp(1,m) * vBar * ...
        ((1+gamma)/((gamma^2)*(pBar^2))) * ...
        ((2*(eigenFreqMatrix(2,pertRIndex)^2))/((eigenFreqMatrix(1,m)^2)-(4*(eigenFreqMatrix(2,pertRIndex)^2)))) * ...
        sin(eigenFreqMatrix(1,m)*timePeriod);

    pIntMatrixEvoCorrection = pIntMatrixEvoCorrection - (1/2) * (perturbationCoeff^2) * ...
        correctionCoeffTemp * ( sin(eigenFreqMatrix(1,m)*timePeriod) ) * pEigenFncCorrectionTemp;

    pEigenFncCorrectionTemp(:,:) = pEigenBasis(3,m,:,:);

    correctionCoeffTemp = pEigenFncPertSquaredDecomp(3,m) * vBar * ...
        ((1+gamma)/((gamma^2)*(pBar^2))) * ...
        ((2*(eigenFreqMatrix(2,pertRIndex)^2))/((eigenFreqMatrix(3,m)^2)-(4*(eigenFreqMatrix(2,pertRIndex)^2)))) * ...
        sin(eigenFreqMatrix(3,m)*timePeriod);

    pIntMatrixEvoCorrection = pIntMatrixEvoCorrection - (1/2) * (perturbationCoeff^2) * ...
        correctionCoeffTemp * ( sin(eigenFreqMatrix(3,m)*timePeriod) ) * pEigenFncCorrectionTemp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end

