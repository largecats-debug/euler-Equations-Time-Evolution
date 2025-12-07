function valueMatrix = TwoD_Disk_TimeEvo_Isentropic_ScalarRecomp(valueMatrixDecomp,zeroFncDecomp,eigenBasis,eigenFncZero, ...
    rBasisSize,thetaBasisSize,rPoints,thetaPoints)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% 2-D time evolution - iserntropic disk %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Scalar eigenfunction recomposition   %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this function is to construct a matrix of values of a given variable
%provided its eigenfunction decomposition in the case of the 2D disc domain

%this function works with valueMatrixDecomp being an N x M matrix
%N = thetaBasisSize, M = rBasisSize
%So the indicies (n,m) specify the theta/r-index

valueMatrix = zeros(thetaPoints,rPoints);
for n = 1:thetaBasisSize

    for m = 1:rBasisSize
        eigenFncMatrixTemp(:,:) = eigenBasis(n,m,:,:);
        valueMatrix = valueMatrix + valueMatrixDecomp(n,m)*eigenFncMatrixTemp;


    end

end

valueMatrix = valueMatrix + zeroFncDecomp*eigenFncZero;

end

