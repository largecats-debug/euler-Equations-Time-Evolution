%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
function [pIntMatrixEvoCorrection,pIntCorrectedDecomp,pIntCorrectedDecompZero] = TwoD_Disk_TimeEvo_Isentropic_InitialCorrection_Radial(pIntMatrixEvo, ...
    pEigenBasis,pEigenFncZero,pBasisWeight,rBasisSize,rGrid,rPoints,pEigenFncPertTemp, ...
    eigenFreqMatrix,perturbationCoeff,pertRIndex,pBar,vBar,sigmaBarConst,gamma,timePeriod)


%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%
%%%\:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:%%\%
%%%\%._.:*~\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%~*:._.%%%\
%%%\%._.:*~\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%~*:._.%%%\
%%%\:._.:*~*%\%\    2-D TIME EVOLUTION - DISC - ISENTROPIC    %\%*~*:._.:%%\%
%%%\:._.:*~*%\%\        INITIAL O(alpha)^2 CORRECTION         %\%*~*:._.:%%\%
%%%\:._.:*~*%\%\         RADIAL EIGENFUNCTION BASES           %\%*~*:._.:%%\%
%%%\%._.:*~\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%~*:._.%%%\
%%%\%._.:*~\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%~*:._.%%%\
%%%\:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:%%\%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%



%this function is to improve the initial guess of p(0,x) = pBar +
% alpha*phi_(0,l) from having error O(alpha^2) to O(alpha^3)
% (so this is specifically for the case that our k-mode radial only)

%The leading order term in f(pBar + alpha*phi_k) that we are trying to
% cancel out to improve the starting solution is:
% (1/2)alpha^2*D^2f(pBar)[phi_k,phi_k]

%Unlike for any other choice of fixed k-mode, in the radial case, the
% corrections will be done with eigenfunctions of the same theta-index,
% which is the entirety of the rest of our basis anyway



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
[pEigenFncPertSquaredDecomp,pEigenFncPertSquaredZeroDecomp] = TwoD_Disk_TimeEvo_Isentropic_InnerDecomp_Radial( ...
    pEigenFncPertSquared,pEigenBasis,pEigenFncZero,pBasisWeight,rBasisSize,rGrid,rPoints,pertRIndex);

pIntMatrixEvoCorrection = pIntMatrixEvo;
pEigenFncCorrectionTemp = zeros(rPoints,1);
pIntCorrectedDecomp = zeros(rBasisSize,1);
pIntCorrectedDecomp(pertRIndex,1) = perturbationCoeff;
pIntCorrectedDecompZero = 2*sqrt(pi)*pBar;


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m = 1:rBasisSize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if m == pertRIndex
        continue
    end

    pEigenFncCorrectionTemp(:) = pEigenBasis(:,m);

    correctionCoeffTemp = pEigenFncPertSquaredDecomp(m) * vBar * ...
        ((1+gamma)/((gamma^2)*(pBar^2))) * ...
        ((2*(eigenFreqMatrix(pertRIndex)^2))/((eigenFreqMatrix(m)^2)-(4*(eigenFreqMatrix(pertRIndex)^2)))) * ...
        sin(eigenFreqMatrix(m)*timePeriod);

    pIntMatrixEvoCorrection = pIntMatrixEvoCorrection - (1/2) * (perturbationCoeff^2) * ...
        correctionCoeffTemp * ( sin(eigenFreqMatrix(m)*timePeriod) ) * pEigenFncCorrectionTemp(:);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

