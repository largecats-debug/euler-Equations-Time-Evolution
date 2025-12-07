function [pEigenBasis,pEigenFncZero,pBasisWeight,vEigenBasis,vEigenFncZero,vBasisWeight,eigenFreqMatrix] = ...
    TwoD_Disk_TimeEvo_NonIsentropicBasis(rBasisSize,thetaBasisSize,sigmaBarProfile,domainRadius,rGridBase, ...
    rPoints,thetaGrid,thetaPoints,pertThetaIndex,pertTrigIndex,radialEntropyLevels,entropyLevelStartingIndex,entropyLevelWidths)


%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%
%\%%%\%%.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   %%\%%%\%
%\%%%\%% \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ /%%\%%%\%
%\%%%\%%  `-/-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   -`-'%%\%%%\%
%\%%%\%%-.___ %\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\% ("""%%\%%%\%
%\%%%\%%   __)%\%\      2-D TIME EVOLUTION - DISC - ISENTROPIC       \%\%) "")%%\%%%\%
%\%%%\%%  (__)%\%\        p,v EIGENFUNCTION BASES GENERATION         \%\%)("")%%\%%%\%
%\%%%\%%  (__,%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%,-"" %%\%%%\%
%\%%%\%%.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   %%\%%%\%
%\%%%\%% \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ /%%\%%%\%
%\%%%\%%  `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   -`-'%%\%%%\%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%

%This function is to create the eigenfunction bases for p and v, 
% in the case of a non-isentropic 2-D disc, D = B_r^2(0), with boundary
% conditions u * n = 0, dp/dn = 0, and a radially piecewise constant
% entropy profile: s(r) = s_i for r_(i-1) < r < r_i (r_0 = 1, r_N = R)

%These eigenfunctions are of the form:
% phi_nm(theta,r) = O_n(theta) * R_(nm)(r), where
%   O_n(theta) = cos(n*theta) or sin(n*theta)
%   R_nm(r) = beta_i * J_n(lambda_nm*sigmaBarProfile_i*[r - gamma_i]) for r_(i-1) < r < r_i

%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %

%The eigenfunctions are essentially described by the Helmholtz equation,
%   and thus the p and v bases would actually be infinite-dimensional. Since
%   the problem is 2-D, the eigenfunctions are indexed by 2 variables, an
%   r-index and a theta-index. The main script sets how many different
%   theta and r indicies to use, N and M. But the theta indicies actually
%   used are the N values 0, k, 2*k, ..., (N-1)*k where k is the
%   theta-index of our chosen k-mode
thetaEigenFncIndicies = (pertThetaIndex-1)*(0:1:thetaBasisSize-1);
rEigenFncIndicies = 1:1:rBasisSize;

%initialize matrix for the (n,m)th eigenvalues
eigenFreqMatrix = zeros(thetaBasisSize,rBasisSize);
JScaledValues = zeros(thetaBasisSize,rBasisSize,rPoints);
JScaledValuesNorm = zeros(thetaBasisSize,rBasisSize);
rEigenFunctionMatrix = zeros(thetaBasisSize,rBasisSize,rPoints);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load pre-calculated data tables for the zeros of J_n(r) and J_n'(r)
temp1 = load('besselFunctionPrecalculatedZeros.mat');
temp2 = fieldnames(temp1);
besselZeroDataTable = temp1.(temp2{1});
clear temp1 temp2

temp1 = load('besselFunctionPrecalculatedDerivativeZeros.mat');
temp2 = fieldnames(temp1);
besselDerivativeZeroDataTable = temp1.(temp2{1});
clear temp1 temp2

%The non-isentropic disc is a case with its own dedicated function for
%   calculating the eigenvalues. Call this function and feed it all the
%   info about the problem and the pre-calculated data tables
[eigenFreqMatrix,gammaMatrix,betaMatrix] = TwoD_Disk_TimeEvo_NonIsentropicEigenvalueCalcBackwards(rPoints,...
    rGridBase,radialEntropyLevels,entropyLevelWidths,entropyLevelStartingIndex,sigmaBarProfile,...
    thetaEigenFncIndicies,rEigenFncIndicies,besselZeroDataTable,besselDerivativeZeroDataTable);

%Now with the eigenvalues, we can calculate gamma_2, gamma_3, ..., gamma_N
%   (for a disc, we start with the assumption that gamma_1 = 0)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:thetaBasisSize

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    for m = 1:rBasisSize

        %This matrix consists of the values of the m-th Bessel function 
        % corresponding to O_n(theta). This means that for O_n(theta) =
        % a_n*cos(n*theta) + b_n*sin(n*theta), the m-th corresponding r
        % eigenfunction is J_(n,m)(r) = J_n(sigmaBarProfile*lambda_(n,m)*r) where
        % J_n denotes the n-th order Bessel function and lambda_(n,m) is such
        % that sigmaBarProfile*lambda_(n,m)*R is the m-th derivative zero. This is 
        % so we have the Neumann cond satisified,
        % J'_n(sigmaBarProfile*lambda*R)=0, where R is the radius of the disc.

        for l = 1:radialEntropyLevels

            for rLoop = entropyLevelStartingIndex(l):entropyLevelStartingIndex(l+1)-1

                JScaledValues(n,m,rLoop) = betaMatrix(n,m,l) * ...
                    besselj(thetaEigenFncIndicies(n),eigenFreqMatrix(n,m) * sigmaBarProfile(l) * ...
                    (rGridBase(rLoop)-gammaMatrix(n,m,l)));
                %weird plot came from having just gammaMatrix(l) or gammaMatrix(n,l)
            end

        end
            JScaledVectorLoop(:) = JScaledValues(n,m,:);
            %JScaledVectorLoop(:) = JScaledValues(n,m+1,:);
            %JScaledValuesNorm(n,m) = sqrt(trapz(rGridBase,(JScaledVectorLoop.*JScaledVectorLoop).*rGridBase));
            %rEigenFunctionMatrix(n,m,:) = JScaledValues(n,m,:)/JScaledValuesNorm(n,m);
            rEigenFunctionMatrix(n,m,:) = JScaledValues(n,m,:);

    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We now have an array of values for each R_(n,m) on rGrid. So we just need
%   to similarly get arrays of values of each O_n on thetaGrid to take the
%   product and find the full p-basis (and v-basis).

thetaEigenFunctionMatrix = zeros(thetaBasisSize,thetaPoints);
 
thetaEigenFunctionMatrix(1,:) = (1/sqrt(2*pi))*cos((0)*ones(1,thetaPoints));



%Check if k-mode is a cosine variant eigenfunction, then the whole
%   basis is cosine variant eigenfunctions
       
%Scaling factor to normalize eigenfunctions different for
%   theta-index = 0 eigenfunctions. So set these (would be n = 1) and
%   then just loop starting from n = 2

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
for n = 2:thetaBasisSize

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if pertTrigIndex == 0

        thetaEigenFunctionMatrix(n,:) = (1/sqrt(pi))*cos(thetaEigenFncIndicies(n)*thetaGrid);
        %thetaEigenFunctionMatrix(n,:) = cos((n-1)*(pertThetaIndex-1)*thetaGrid);

    else

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        if mod(n,2) == 0
            
            thetaEigenFunctionMatrix(n,:) = (1/sqrt(pi))*cos(thetaEigenFncIndicies(n)*thetaGrid);

        end

        if mod(n,2) == 1

            thetaEigenFunctionMatrix(n,:) = (1/sqrt(pi))*sin(thetaEigenFncIndicies(n)*thetaGrid);

        end
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%The pEigenBasis matrix and the others will have dimensions
% N x M x thetaPoints x rPoints
% the first 2 indicies specify the (n,m) eigenfunction we want the
% matrix of values for

%So pEigenBasis(2,3,:,:) gives us the matrix of values of the appropriate 
% version of the p-basis eigenfunction phi_(2,3)/sigmaBarProfile on the
% discretization of the disk rGrid x thetaGrid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if pertThetaIndex ~= 1

    pEigenBasis = zeros(thetaBasisSize,rBasisSize,thetaPoints,rPoints);

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    for n = 1:thetaBasisSize

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for m = 1:rBasisSize

            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            for i = 1:thetaPoints

                for l = 1:radialEntropyLevels

                    %pEigenBasis(n,m,i,entropyLevelStartingIndex(l):entropyLevelStartingIndex(l+1)-1) = 0 + ...
                    %    (rEigenFunctionMatrix(n,m,entropyLevelStartingIndex(l):entropyLevelStartingIndex(l+1)-1) * ...
                    %    thetaEigenFunctionMatrix(n,i) / sigmaBarProfile(l));
                    pEigenBasis(n,m,i,entropyLevelStartingIndex(l):entropyLevelStartingIndex(l+1)-1) = 0 + ...
                        rEigenFunctionMatrix(n,m,entropyLevelStartingIndex(l):entropyLevelStartingIndex(l+1)-1) * ...
                        thetaEigenFunctionMatrix(n,i);
                    %vEigenBasis(n,m,i,entropyLevelStartingIndex(l):entropyLevelStartingIndex(l+1)-1) = 0 + ...
                    %     pEigenBasis(n,m,i,entropyLevelStartingIndex(l):entropyLevelStartingIndex(l+1)-1) ...
                    %    * (sigmaBarProfile(l)^2);
                    vEigenBasis(n,m,i,entropyLevelStartingIndex(l):entropyLevelStartingIndex(l+1)-1) = 0 + ...
                         pEigenBasis(n,m,i,entropyLevelStartingIndex(l):entropyLevelStartingIndex(l+1)-1) ...
                        * (sigmaBarProfile(l)^2);
                end

            end
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    %J_0(0) = 1, J'_0(0) = 0 while J_n(0) = 0, J'_n(0) ~= 0 for all n >= 1, 
    %   so we will have one family of radially constant eigenfunctions
    %   one of these also has constant angular part, so we treat the single
    %   constant eigenfunction separately. Normalized, this eigenfunction 
    %   will just be 1/|D|
    for i0 = 1:thetaPoints
        for l0 = 1:radialEntropyLevels
            for q0 = entropyLevelStartingIndex(l0):entropyLevelStartingIndex(l0+1)-1
            pEigenFncZero(i0,q0) = 0 + ...
                (sqrt(2))*(1/(sqrt(pi)*domainRadius))* ...
                (1/sigmaBarProfile(l0));
            vEigenFncZero(i0,q0) = 0 + ...
                pEigenFncZero(i0,q0) * ...
                (sigmaBarProfile(l0)^2);
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    for m = 1:rBasisSize

        for l1 = 1:radialEntropyLevels

            pEigenBasis(1,m,1,entropyLevelStartingIndex(l):entropyLevelStartingIndex(l+1)-1) = 0 + ...
                rEigenFunctionMatrix(1,m,entropyLevelStartingIndex(l):entropyLevelStartingIndex(l+1)-1) .* ...
                thetaEigenFunctionMatrix(1,:) / sigmaBarProfile(l1);
            vEigenBasis(1,m,1,entropyLevelStartingIndex(l):entropyLevelStartingIndex(l+1)-1) = 0 + ...
                pEigenBasis(1,m,1,entropyLevelStartingIndex(l):entropyLevelStartingIndex(l+1)-1) * ...
                (sigmaBarProfile(l1)^2);

            pEigenFncZero(1,entropyLevelStartingIndex(l):entropyLevelStartingIndex(l+1)-1) = 0 + ...
                (1/(sqrt(pi)*domainRadius)) * (1/sigmaBarProfile(l1)) * ...
                ones(entropyLevelStartingIndex(l):entropyLevelStartingIndex(l+1)-1);
            vEigenFncZero(1,entropyLevelStartingIndex(l):entropyLevelStartingIndex(l+1)-1) = 0 + ...
                pEigenFncZero(1,entropyLevelStartingIndex(l):entropyLevelStartingIndex(l+1)-1) * ...
                (sigmaBarProfile(l1)^2);

        end

    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for l = 1:radialEntropyLevels
    pBasisWeight(entropyLevelStartingIndex(l):entropyLevelStartingIndex(l+1)-1) = sigmaBarProfile(l)^2;
    vBasisWeight(entropyLevelStartingIndex(l):entropyLevelStartingIndex(l+1)-1) = sigmaBarProfile(l)^(-2);
end



for n1 = 1:thetaBasisSize

    for m1 = 1:rBasisSize
        pEigenFncTemp = zeros(thetaPoints,rPoints);
        pEigenFncTemp(:,:) = pEigenBasis(n1,m1,:,:);

        for i1 = 1:thetaPoints
            pEigenFncTempRay(:) = pEigenFncTemp(i1,:);
            rayVector(i1) = simps2DCopy(rGridBase,((pEigenFncTempRay.*pEigenFncTempRay).*rGridBase).*pBasisWeight);
        end

        pEigenFncTempNorm = simps2DCopy(thetaGrid,rayVector);

        pEigenBasis(n1,m1,:,:) = pEigenBasis(n1,m1,:,:)/sqrt(pEigenFncTempNorm);
        
        for l1 = 1:radialEntropyLevels
            vEigenBasis(n1,m1,:,entropyLevelStartingIndex(l):entropyLevelStartingIndex(l+1)-1) = pEigenBasis(n1,m1,:,entropyLevelStartingIndex(l):entropyLevelStartingIndex(l+1)-1)*(sigmaBarProfile(l)^2);
        end

    end

end




%The u-basis eigenfunctions would be a 2D matrix with 2D vector components, 
% however the u-eigenfunctions are not actually required for the NL 
% evolution. So for now the u-basis will not be calculated




end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


