function [pMatrixTemp,pMatrixDecompTemp,pZeroDecompTemp,vMatrixTemp] = TwoD_Disk_TimeEvo_Isentropic_VtoP_Radial( ...
    decompMatrix,zeroFncDecomp,pEigenBasis,pEigenFncZero,pBasisWeight,vEigenBasis,vEigenFncZero, ...
    rBasisSize,rGrid,rPoints,gamma,pertRIndex,pBar,vBar)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%\:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:%%\%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% 2-D time evolution - iserntropic disk %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%      v(t,x) to p(t,x) calculator      %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                RADIAL VARIANT                 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%\:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:%%\%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%This is the 2D version of a function to take in the eigenfunction
%decomposition of v and return both the matrix of actual values of p and
%the eigenfunction decomposition of p. This version of the function is
%specifically for the case where the variables can be treated as functions
%of the radius only (when pertThetaIndex = 0)

%this function works with vMatrixDecomp being an N x M x 2 matrix
%   N = thetaBasisSize, M = rBasisSize
%   The (n,m) of the first two indicies specify the theta/r-index

%Construct the matrix of values of v on the disc by recomposing it as a
%linear combination of its basis eigenfunctions
vMatrixTemp = zeros(rPoints,1);
vMatrixTemp = TwoD_Disk_TimeEvo_Isentropic_ScalarRecomp_Radial(decompMatrix,zeroFncDecomp,...
    vEigenBasis,vEigenFncZero,rBasisSize,rPoints);


pMatrixTemp = zeros(rPoints,1);
pMatrixTemp(:) = pBar*(vMatrixTemp(:).^((-1)*gamma))*(vBar^(gamma));


[pMatrixDecompTemp,pZeroDecompTemp] = TwoD_Disk_TimeEvo_Isentropic_ScalarDecomp_Radial( ...
    pMatrixTemp,pEigenBasis,pEigenFncZero,pBasisWeight,rBasisSize,rGrid,rPoints,pertRIndex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

