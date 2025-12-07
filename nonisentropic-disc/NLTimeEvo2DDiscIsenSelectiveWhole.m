function [pNLMatrix,pNLDecompMatrix,pNLZeroDecomp,vNLMatrix,vNLDecompMatrix,vNLZeroDecomp,uNLDecompMatrix] = NLTimeEvo2DDiscIsenSelectiveWhole( ...
    pIntMatrix,uIntDecompMatrix,rBasisSize,rGrid,rPoints,thetaBasisSize,thetaGrid,thetaPoints,timeGrid,timePoints,timeStepSize,pEigenBasis, ...
    vEigenBasis,pEigenFncZero,vEigenFncZero,eigenFreqMatrix,entropyConst,sigmaBarConst,gamma,cP,cV,pertRIndex)

%this is the higher level function to compute the nonlinear evolution (in
%time) of the Euler equations with piecewise constant entropy. This
%function is provided:

%NOTE MAY 1:
% this is the version of the nonlinear evolution meant for the simplified
% version for the version of the problem that does not require the sin and
% cos variant of every eigenfunction. So pEigenBasis for instance has gone
% from being a 5D matrix with dimensions 
% thetaSize x rSize x 2 x thetaPoints x rPoints to a 4D matrix with
% dimensions thetaSize x rSize x 2 x thetaPoints x rPoints
% similarly any varialbe that had a column that was just 1 or 2 for the
% cos/sin variant of the eigenfunctions is one dimension smaller

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initial data: pIntArray = p(x,0) and uIntArray = u(x,0) (but this code
%is being written at the moment assuming u_0(x) = 0) as doubles
%the entropyLevels and the already calculated corresponding sigmaBar's
%the structure of the entropy profile (the widths/locations of jumps)
%the timeGrid over which to perform the NL evolution
%the xGrid we are using/the int data is specified on
%the eigenfunction basis for p,v (phi) and u (psi)

%for now, this code will be written with a direct inner product
%decomposition into the eigenfunctions in mind, at least until I can
%possibly figure out howt to properly use the nonuniform FFT

pNLMatrix = zeros(timePoints,thetaPoints,rPoints);
%uNLMatrix = zeros(timePoints,rPoints,thetaPoints);
vNLMatrix = zeros(timePoints,thetaPoints,rPoints);

%the NL decomposition matricies will be set up so that each column of the
%matrix is the vertical vector of the eigenfunction coefficents at a fixed
%time step. That is, the n_th column will be (f_0,f_1,...,f_basisSize)^T
%can make the u decomposition matrix with or without the 0-eigenfunction
%component, as again, psi_0 = 0, so it would just be a row of zeros
pNLDecompMatrix = zeros(thetaBasisSize,rBasisSize,timePoints);
vNLDecompMatrix = zeros(thetaBasisSize,rBasisSize,timePoints);
uNLDecompMatrix = zeros(thetaBasisSize,rBasisSize,timePoints);

%pNLZeroDecomp = zeros(1,timePoints);
%vNLZeroDecomp = zeros(1,timePoints);
pNLZeroDecomp = zeros(timePoints,1);
vNLZeroDecomp = zeros(timePoints,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for the things within this subsection, these can probably be moved to the
%higher order program that will actually be running the iteration/that
%will be calling the evolution functions. As the eigenfunction matricies,
%the jumpLocationOnGrid array, the entropyXGridLevel array, the innerProdWeight
%will be used every time the evolution functions are called. So we can avoid 
%computingb the same thing multiple times by calculating these outside this
%function and feeding them as yet another input


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pNLMatrix(1,:,:) = pIntMatrix;
vNLMatrix(1,:,:) = (pNLMatrix(1,:,:).^((-1)*(1/gamma)))*(exp(entropyConst/cP));

matrixTemp(:,:) = pNLMatrix(1,:,:);
%[pNLDecompMatrix(:,:,1),pNLZeroDecomp(1)] = eigenFncDecompInnerProdDiscSelective(matrixTemp,pEigenBasis, ...
 %   pEigenFncZero,rBasisSize,thetaBasisSize,rGrid,rPoints,thetaGrid,thetaPoints);
[pNLDecompMatrix(:,:,1),pNLZeroDecomp(1)] = eigenFncDecompInnerProdDiscSelective2(matrixTemp,pEigenBasis, ...
    pEigenFncZero,rBasisSize,thetaBasisSize,rGrid,rPoints,thetaGrid,thetaPoints,pertRIndex);
matrixTemp(:,:) = vNLMatrix(1,:,:);
%[vNLDecompMatrix(:,:,1),vNLZeroDecomp(1)] = eigenFncDecompInnerProdDiscSelective(matrixTemp,vEigenBasis, ...
 %   vEigenFncZero,rBasisSize,thetaBasisSize,rGrid,rPoints,thetaGrid,thetaPoints);
[vNLDecompMatrix(:,:,1),vNLZeroDecomp(1)] = eigenFncDecompInnerProdDiscSelective2(matrixTemp,vEigenBasis, ...
    vEigenFncZero,rBasisSize,thetaBasisSize,rGrid,rPoints,thetaGrid,thetaPoints,pertRIndex);
uNLDecompMatrix(:,:,1) = zeros(thetaBasisSize,rBasisSize);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%to use a leapfrog scheme, need a second piece of initial data
%call function that takes in u,p,v arrays at t=0 and uses RK4 type method
%to initialize u,p,v at t=dt
%NOTE: need only give RK4 initial data extension function the 1,...,N
%eigenfunctions of the p-basis as the 0-eigenfunction component of p will
%contribute nothing when calculating d_x p which is how the RK4 code
%currently works

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%[uNLMatrix(2,:),pNLMatrix(2,:),vNLMatrix(2,:)] = NlTimeEvo1DRK4Ext(uNLMatrix(1,:),pNLMatrix(1,:),vNLMatrix(1,:), ...
%    pEigenBasis,pEigenBasisDx,uEigenBasis,uEigenBasisDx,xGrid,xPoints,timeStepSize, ...
%    basisSize,pInnerProdWeight,uInnerProdWeight,entropyXGridLevels,gamma,cP);

%pNLDecompMatrix(:,2) = eigenFncDecompInnerProd1D(pNLMatrix(2,:),pEigenBasis,basisSize+1, ...
%    pInnerProdWeight,xGrid);
%vNLDecompMatrix(:,2) = eigenFncDecompInnerProd1D(vNLMatrix(2,:),vEigenBasis,basisSize+1, ...
%    vInnerProdWeight,xGrid);
%uNLDecompMatrix(:,2) = eigenFncDecompInnerProd1D(uNLMatrix(2,:),uEigenBasis,basisSize, ...
%    uInnerProdWeight,xGrid);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%APRIL 10:
%Looking for error in code causing complex outputs
%presumably p/v going negative at some point
%for now: implent single euler step initial data extension

%vNLDecompMatrix(1,2) = vNLDecompMatrix(1,1);
%for k = 1:basisSize
%    uNLDecompMatrix(k,2) = uNLDecompMatrix(k,1) + timeStepSize*eigenFreqArray(k)*pNLDecompMatrix(k+1,1);
%    vNLDecompMatrix(k+1,2) = vNLDecompMatrix(k+1,1) + timeStepSize*eigenFreqArray(k)*uNLDecompMatrix(k,1);
%end

%APRIL 16:
%for extension of initial data will initially hard code assumption that
%this will only encounter u_0 = 0 so that every decomp coefficent of u at
%t=0 is zero
vNLDecompMatrix(:,:,2) = vNLDecompMatrix(:,:,1);
pNLDecompMatrix(:,:,2) = pNLDecompMatrix(:,:,1);
pNLZeroDecomp(2) = pNLZeroDecomp(1);
vNLZeroDecomp(2) = vNLZeroDecomp(1);
pNLMatrix(2,:,:) = pNLMatrix(1,:,:);
vNLMatrix(2,:,:) = vNLMatrix(1,:,:);



for n = 1:thetaBasisSize

    for m = 1:rBasisSize

        %for q = 1:2

            %uNLDecompMatrix(n,m,2) = uNLDecompMatrix(n,m,1) + ...
            %    timeStepSize*(sigmaBarConst^2)*eigenFreqMatrix(n,m)*pNLDecompMatrix(n,m,1);

            %NOTE MAY 12:
            % trying to move both sigma terms from u eqn to both
            % equations, not sure of effect between different EVPs
            uNLDecompMatrix(n,m,2) = uNLDecompMatrix(n,m,1) + ...
                timeStepSize*(sigmaBarConst)*eigenFreqMatrix(n,m)*pNLDecompMatrix(n,m,1);

        %end
        
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%[pNLMatrix(2,:),pNLDecompMatrix(:,2),vNLMatrix(2,:)] = timeEvoVtoP(vNLDecompMatrix(:,2),pEigenBasis,vEigenBasis, ...
   %pInnerProdWeight,basisSize+1,xGrid,xPoints,entropyXGridLevels,gamma,cP);

%the function that takes in v's decomposition array and returns p's
%decomposition array also returns the actual values of p and v at the
%current time step, so to store of the NL variable matricies we need
%only compute u's here:

%for k = 1:basisSize
%    uNLMatrix(2,:) = uNLMatrix(2,:) + uNLDecompMatrix(k,2)*uEigenBasis(k,:);
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 3:timePoints

    vNLZeroDecomp(j) = vNLZeroDecomp(j-1);

    %perform actual leapfrog step on the decomposition components
    %from the weak* solution of the NL system this iterates a finite
    %difference/spectral version of
    % d_t u_j - lambda_j * p_j = 0
    % d_t v_j - lambda_j * u_j = 0

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    
    for n = 1:thetaBasisSize

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

        for m = 1:rBasisSize

            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

            %for q = 1:2

                %uNLDecompMatrix(n,m,q,j) = uNLDecompMatrix(n,m,q,j-2) + ...
                %2*timeStepSize*(sigmaBarConst^2)*eigenFreqMatrix(n,m)*pNLDecompMatrix(n,m,q,j-1);
                %
                %vNLDecompMatrix(n,m,q,j) = vNLDecompMatrix(n,m,q,j-2) + ...
                %2*timeStepSize*eigenFreqMatrix(n,m)*uNLDecompMatrix(n,m,q,j-1);

                %uNLDecompMatrix(n,m,j) = uNLDecompMatrix(n,m,j-2) + ...
                %2*timeStepSize*(sigmaBarConst^2)*eigenFreqMatrix(n,m)*pNLDecompMatrix(n,m,j-1);
                %%%%%%%%%%%%%%%
                %vNLDecompMatrix(n,m,j) = vNLDecompMatrix(n,m,j-2) + ...
                %2*timeStepSize*eigenFreqMatrix(n,m)*uNLDecompMatrix(n,m,j-1);

                %NOTE MAY 12:
                % trying to move both sigma terms from u eqn to both
                % equations, not sure of effect between different EVPs
                uNLDecompMatrix(n,m,j) = uNLDecompMatrix(n,m,j-2) + ...
                2*timeStepSize*(sigmaBarConst)*eigenFreqMatrix(n,m)*pNLDecompMatrix(n,m,j-1);
                vNLDecompMatrix(n,m,j) = vNLDecompMatrix(n,m,j-2) + ...
                2*timeStepSize*eigenFreqMatrix(n,m)*sigmaBarConst*uNLDecompMatrix(n,m,j-1);

            %end

            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

        end

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    end

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    [pNLMatrix(j,:,:),pNLDecompMatrix(:,:,j),pNLZeroDecomp(j),vNLMatrix(j,:,:)] = ...
        timeEvo2DDiscIsenVtoPSelective(vNLDecompMatrix(:,:,j),vNLZeroDecomp(j),pEigenBasis, ...
        vEigenBasis,pEigenFncZero,vEigenFncZero,rBasisSize,thetaBasisSize,rGrid,rPoints, ...
        thetaGrid,thetaPoints,entropyConst,gamma,cP,pertRIndex);




end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end

