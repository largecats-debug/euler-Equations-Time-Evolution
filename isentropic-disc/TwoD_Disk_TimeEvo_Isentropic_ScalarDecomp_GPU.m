
function [matrixDecomp,zeroDecomp] = TwoD_Disk_TimeEvo_Isentropic_ScalarDecomp_GPU( ...
    valueMatrix,eigenBasis,eigenFncZero,basisWeight,rBasisSize,thetaBasisSize, ...
    rGrid,rPoints,thetaGrid,thetaPoints,pertRIndex)

%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%
%%%\:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*%\%%
%%%\:._.:*~\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%*:._.:*%\%%
%%%\:._.:*~\%\%\                                              %\%\%*:._.:*%\%%
%%%\:._.:*~\%\%\    2-D TIME EVOLUTION - DISC - ISENTROPIC    %\%\%*:._.:*%\%%
%%%\:._.:*~\%\%\      SCALAR EIGENFUNCTION DECOMPOSITION      %\%\%*:._.:*%\%%
%%%\:._.:*~\%\%\                                              %\%\%*:._.:*%\%%
%%%\:._.:*~\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%*:._.:*%\%%
%%%\:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*%\%%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%


%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %


%This function takes in:
% 1. a matrix of values on the 2D disk
% 2. a basis of eigenfunctions on the disk in some insentropic case

%This function returns:
% an N x M matrix, where the indicies designate the theta (n) index and 
% r (m) index, respectively.

%NOTE 9/17/25:
% need separate scalar eigenfunction decomposition function for the case of
% a non-isentropic disk. The weight which makes the basis into an
% orthogonal basis depends on sigmaBar (and thus s)



%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


matrixDecomp = zeros(thetaBasisSize,rBasisSize);

%To increase accuracy, as we find the decomposition coefficients for the
% valueMatrix corresponding to 1. phi_0 2. whatever eigenfunction phi_nm
% is fixed and subtract that part from the valueMatrix so that the inner 
% product of the valueMatrix with the next eigenfunction isn't affected
% by it. (these two decomposition coefficients will be much larger
% magnitude than the rest)


rGridGPU = gpuArray(rGrid);
rGridGPUWeighted = gpuArray(rGrid*basisWeight);
thetaGridGPU = gpuArray(thetaGrid);


%First, find the decomposition coefficent for the (0,0)-eigenfunction (the
% constant eigenfunction on the disc) separatly
eigenBasisTemp = gpuArray(zeros(thetaPoints,rPoints));
valueMatrixTemp = gpuArray(zeros(thetaPoints,rPoints));
valueBasisIntegral = gpuArray(zeros(thetaPoints,1));


eigenBasisTemp(:,:) = gpuArray(eigenFncZero(:,:));
valueMatrixTemp(:,:) = gpuArray(valueMatrix(:,:));
valueBasisIntegral(:) = arrayfun(@trapz2,rGridGPU,(eigenBasisTemp.*valueMatrixTemp).*rGridGPUWeighted,2);


zeroDecomp = arrayfun(@trapz2,thetaGridGPU,valueBasisIntegral);
%Subtract portion of given matrix in the phi_0 direction to reduce
% interference of this larger magnitude term from smaller ones
valueMatrixTemp = sum(valueMatrixTemp, (-1)*zeroDecomp*eigenBasisTemp);


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


%Second, find the decomposition coefficent for the (0,0)-eigenfunction (the
% constant eigenfunction on the disc) separatly
eigenBasisTemp(:,:) = gpuArray(eigenBasis(2,pertRIndex,:,:));
valueBasisIntegral(:) = arrayfun(@trapz2,rGridGPU,(eigenBasisTemp.*valueMatrixTemp).*rGridGPUWeighted,2);
matrixDecomp(2,pertRIndex) = arrayfun(@trapz2,thetaGridGPU,valueBasisIntegral);

%Again, subtract the portion of the given matrix in the direction of phi_k
% to reduce interference with smaller magnitude terms
valueMatrixTemp = sum(valueMatrixTemp, (-1)*matrixDecomp(2,pertRIndex)*eigenBasisTemp);


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


%Iterate over each r-index first. We put this one on the
% outside because the eigenfunction with larger r-indicies behave in a more
% stable manner if so.
for m = 1:rBasisSize

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    %iterate over each theta-index second
    for n = 1:thetaBasisSize

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

        if n == thetaBasisSize-1 && m == pertRIndex
            continue
        end

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

        %Second, find the decomposition coefficent for the (0,0)-eigenfunction (the
        % constant eigenfunction on the disc) separatly
        eigenBasisTemp(:,:) = gpuArray(eigenBasis(n,m,:,:));
        valueBasisIntegral(:) = arrayfun(@trapz2,rGridGPU,(eigenBasisTemp.*valueMatrixTemp).*rGridGPUWeighted,2);
        matrixDecomp(n,m) = arrayfun(@trapz2,thetaGridGPU,valueBasisIntegral);

        %Again, subtract the portion of the given matrix in the direction of phi_k
        % to reduce interference with smaller magnitude terms
        valueMatrixTemp = sum(valueMatrixTemp, (-1)*matrixDecomp(n,m)*eigenBasisTemp);

    end

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





end
