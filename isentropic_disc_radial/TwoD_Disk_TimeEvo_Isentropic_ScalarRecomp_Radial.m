function valueMatrix = TwoD_Disk_TimeEvo_Isentropic_ScalarRecomp_Radial(valueMatrixDecomp,zeroFncDecomp,eigenBasis,eigenFncZero, ...
    rBasisSize,rPoints)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% 2D time evolution - RADIAL - iserntropic disk %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Scalar eigenfunction recomposition   %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this function is to construct a matrix of values of a given variable
%provided its eigenfunction decomposition in the case of the 2D disc domain

%this function works with valueMatrixDecomp being an N x M matrix
%N = thetaBasisSize, M = rBasisSize
%So the indicies (n,m) specify the theta/r-index

valueMatrix = zeros(rPoints,1);
eigenFncMatrixTemp = zeros(rPoints,1);

for m = 1:rBasisSize

    eigenFncMatrixTemp(:) = eigenBasis(:,m);
    valueMatrix(:) = valueMatrix(:) + valueMatrixDecomp(m)*eigenFncMatrixTemp(:);

end

valueMatrix(:) = valueMatrix(:) + zeroFncDecomp*eigenFncZero(:);

end

