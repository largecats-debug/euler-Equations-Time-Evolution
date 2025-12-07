
function [pMatrixTemp,pMatrixDecompTemp,pZeroDecompTemp,vMatrixTemp] = TwoD_Disk_TimeEvo_Isentropic_VtoP_GPU( ...
    decompMatrix,zeroFncDecomp,pEigenBasis,pEigenFncZero,pBasisWeight,vEigenBasis,vEigenFncZero, ...
    rBasisSize,thetaBasisSize,rGrid,rPoints,thetaGrid,thetaPoints,gamma,pertRIndex,pertThetaIndex,pBar,vBar)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%\:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:%%\%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% 2-D time evolution - iserntropic disk %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%      v(t,x) to p(t,x) calculator      %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%\:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:%%\%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tStart_GPUtimer = tic;

%This is the 2D version of a function to take in the eigenfunction
%decomposition of v and return both the matrix of actual values of p and
%the eigenfunction decomposition of p

%this function works with vMatrixDecomp being an N x M x 2 matrix
%   N = thetaBasisSize, M = rBasisSize
%   The (n,m) of the first two indicies specify the theta/r-index

%Construct the matrix of values of v on the disc by recomposing it as a
%linear combination of its basis eigenfunctions
vMatrixTemp(:,:) = zeros(thetaPoints,rPoints);
vMatrixTemp = TwoD_Disk_TimeEvo_Isentropic_ScalarRecomp(decompMatrix,zeroFncDecomp,...
    vEigenBasis,vEigenFncZero,rBasisSize,thetaBasisSize,rPoints,thetaPoints);

pMatrixTemp(:,:) = zeros(thetaPoints,rPoints);
pMatrixTemp = pBar*(vMatrixTemp.^((-1)*gamma))*(vBar^(gamma));

%[pMatrixDecompTemp,pZeroDecompTemp] = TwoD_Disk_TimeEvo_IsentropicScalarDecomposition( ...
%    pMatrixTemp,pEigenBasis,pEigenFncZero,pBasisWeight,rBasisSize,thetaBasisSize, ...
%    rGrid,rPoints,thetaGrid,thetaPoints,pertRIndex,pertThetaIndex);

%[pMatrixDecompTemp,pZeroDecompTemp] = eigenFncDecompInnerProdDiscSelective2( ...
%    pMatrixTemp,pEigenBasis,pEigenFncZero,rBasisSize,thetaBasisSize, ...
%    rGrid,rPoints,thetaGrid,thetaPoints,pertRIndex);
[pMatrixDecompTemp,pZeroDecompTemp] = TwoD_Disk_TimeEvo_Isentropic_ScalarDecomp_GPU( ...
    pMatrixTemp,pEigenBasis,pEigenFncZero,pBasisWeight,rBasisSize,thetaBasisSize, ...
    rGrid,rPoints,thetaGrid,thetaPoints,pertRIndex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elapsed_CPUexecution = toc(tStart_GPUtimer)

end

