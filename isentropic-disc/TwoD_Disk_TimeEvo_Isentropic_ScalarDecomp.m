function [matrixDecomp,zeroDecomp] = TwoD_Disk_TimeEvo_Isentropic_ScalarDecomp( ...
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



%First, find the decomposition coefficent for the (0,0)-eigenfunction (the
% constant eigenfunction on the disc) separatly

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for l = 1:thetaPoints       %% ZERO MODE DECOMP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    eigenBasisRayTemp(:) = eigenFncZero(l,:);
    valueMatrixRayTemp(:) = valueMatrix(l,:);

    valueBasisRayInt(l) = simps2DCopy(rGrid,(eigenBasisRayTemp.*valueMatrixRayTemp).*rGrid*basisWeight);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%zeroDecomp = (pertThetaIndex-1)*simps2DCopy(thetaGrid,valueBasisRayInt);
zeroDecomp = simps2DCopy(thetaGrid,valueBasisRayInt);

%Subtract portion of given matrix in the phi_0 direction to reduce
% interference of this larger magnitude term from smaller ones 
valueMatrix = valueMatrix - zeroDecomp*eigenFncZero;


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

eigenFncTemp(:,:) = eigenBasis(2,pertRIndex,:,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for l = 1:thetaPoints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    eigenBasisRayTemp(:) = eigenFncTemp(l,:);
    valueMatrixRayTemp(:) = valueMatrix(l,:);

    %valueBasisRayInt(l) = trapz(rGrid,(eigenBasisRayTemp.*valueMatrixRayTemp).*rGrid);
    valueBasisRayInt(l) = simps2DCopy(rGrid,(eigenBasisRayTemp.*valueMatrixRayTemp).*rGrid*basisWeight);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrixDecomp(2,pertRIndex) = simps2DCopy(thetaGrid,valueBasisRayInt);

%Again, subtract the portion of the given matrix in the direction of phi_k
% to reduce interference with smaller magnitude terms
%valueMatrix = valueMatrix - matrixDecomp(2,pertRIndex)*eigenFncTemp;
valueMatrix = valueMatrix - matrixDecomp(2,pertRIndex)*eigenFncTemp;


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

        %iterate over each thetaGrid element (over each ray)
        for l = 1:thetaPoints

            %eigenBasisRayTemp(:) = eigenBasis(n,m,l,:);
            %valueMatrixRayTemp(:) = valueMatrix(l,:);
            %valueBasisRayInt(l) = trapz(rGrid,(eigenBasisRayTemp.*valueMatrixRayTemp).*rGrid);

            eigenBasisRayTemp(:) = eigenBasis(thetaBasisSize-n+1,m,l,:);
            valueMatrixRayTemp(:) = valueMatrix(l,:);
            valueBasisRayInt(l) = simps2DCopy(rGrid,(eigenBasisRayTemp.*valueMatrixRayTemp).*rGrid*basisWeight);

        end

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

        %matrixDecomp(n,m) = trapz(thetaGrid,valueBasisRayInt);
        matrixDecomp(thetaBasisSize-n+1,m) = simps2DCopy(thetaGrid,valueBasisRayInt);

        eigenFncTemp(:,:) = eigenBasis(thetaBasisSize-n+1,m,:,:);
        %valueMatrix = valueMatrix - matrixDecomp(thetaBasisSize-n+1,m)*eigenFncTemp;
        valueMatrix = valueMatrix - matrixDecomp(thetaBasisSize-n+1,m)*eigenFncTemp;

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    end

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





end
