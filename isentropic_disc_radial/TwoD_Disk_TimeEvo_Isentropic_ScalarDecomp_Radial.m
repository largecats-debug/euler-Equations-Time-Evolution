function [matrixDecomp,zeroDecomp] = TwoD_Disk_TimeEvo_Isentropic_ScalarDecomp_Radial( ...
    valueMatrix,eigenBasis,eigenFncZero,basisWeight,rBasisSize,rGrid,rPoints,pertRIndex)


%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%
%%%\:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*%\%%
%%%\:._.:*~\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%*:._.:*%\%%
%%%\:._.:*~\%\%\                                              %\%\%*:._.:*%\%%
%%%\:._.:*~\%\%\    2-D TIME EVOLUTION - DISC - ISENTROPIC    %\%\%*:._.:*%\%%
%%%\:._.:*~\%\%\      SCALAR EIGENFUNCTION DECOMPOSITION      %\%\%*:._.:*%\%%
%%%\:._.:*~\%\%\                RADIAL VARIANT                %\%\%*:._.:*%\%%
%%%\:._.:*~\%\%\                                              %\%\%*:._.:*%\%%
%%%\:._.:*~\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%*:._.:*%\%%
%%%\:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*~*:._.:*%\%%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%


%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %



%NOTE MAY 1:
% this is the version of the eigenfunction decomposition via inner product
% function for the version of the problem that does not require a cos/sin
% variant of every eigenfunction.

%this function is to return the eigenfunction decomposition in the case of
%a 2D disc domain. The returned matrix will be a N x M x 2 matrix, where
%the first two indicies designate the theta (n) index and r (m) index and
%the third component is a 1 (for the cos variant) or 2 (for the sin
%variant)

%this version of the 2D eigenfunction decomposition also assumes that we
%are working in an isentropic case, s = s_0 constant. Can perform such an
%inner product decomposition for nonconstant entropy, however will take
%some work to specify how to use the inner product weight with which the
%eigenbasis is orthogonal.
%NOTE: perhaps can specify inner product weight separately for the
%r-integral and the theta-integral. Seems like it would be much easier to
%have the entropy be s = s(r,theta) in this regard instead of s=s(x,y)

matrixDecomp = zeros(rBasisSize,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%find the decomposition coefficent for the (0,0)-eigenfunction (the
%constant eigenfunction on the disc) separatly
zeroDecomp = simps2DCopy(rGrid,(eigenFncZero.*valueMatrix).*rGrid*basisWeight);

%subtracting the 0-eigenfunction component to minimize interference
%valueMatrix = valueMatrix - zeroDecomp*eigenFncZero;
valueMatrix = valueMatrix - zeroDecomp*eigenFncZero;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%find the decomposition coefficent for the k-eigenfunction (the
%perturbed eigenfunction) separatly as well
eigenFncTemp = zeros(rPoints,1);
eigenFncTemp(:) = eigenBasis(:,pertRIndex);
valueBasisRayInt = simps2DCopy(rGrid,(eigenFncTemp.*valueMatrix).*rGrid*basisWeight);


%again subtract component from valueMatrix to minimize interference
matrixDecomp(pertRIndex) = valueBasisRayInt;
%valueMatrix = valueMatrix - matrixDecomp(pertRIndex,1)*eigenFncTemp;
valueMatrix = valueMatrix - matrixDecomp(pertRIndex)*eigenFncTemp;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%iterate over each r-eigenfunction index

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m = 1:rBasisSize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %iterate over the cos and sin variants
        %for q = 1:2

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        if m == pertRIndex
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

            continue

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        end
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

        eigenFncTemp(:) = eigenBasis(:,m);
        valueBasisRayInt = simps2DCopy(rGrid,(eigenFncTemp.*valueMatrix).*rGrid*basisWeight);

        matrixDecomp(m) = valueBasisRayInt;
        %valueMatrix = valueMatrix - matrixDecomp(m)*eigenFncTemp;
        valueMatrix = valueMatrix - matrixDecomp(m)*eigenFncTemp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

