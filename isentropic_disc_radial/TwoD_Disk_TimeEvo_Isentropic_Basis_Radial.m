function [pEigenBasis,pEigenFncZero,pBasisWeight,vEigenBasis,vEigenFncZero,vBasisWeight,eigenFreqMatrix,JZeroDerivLocations] = ...
    TwoD_Disk_TimeEvo_Isentropic_Basis_Radial(rBasisSize,sigmaBar,domainRadius,rGridBase,rPoints)


%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%
%\%%%\%%.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   %%\%%%\%
%\%%%\%% \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ /%%\%%%\%
%\%%%\%%  `-/-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   -`-'%%\%%%\%
%\%%%\%%-.___ %\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\% ("""%%\%%%\%
%\%%%\%%   __)%\%\      2-D TIME EVOLUTION - DISC - ISENTROPIC       \%\%) "")%%\%%%\%
%\%%%\%%  (__,%\%\                  RADIAL VARIANT                   \%\%,-"" %%\%%%\%
%\%%%\%%  (__)%\%\        p,v EIGENFUNCTION BASES GENERATION         \%\%)("")%%\%%%\%
%\%%%\%%  (__,%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%,-"" %%\%%%\%
%\%%%\%%.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   %%\%%%\%
%\%%%\%% \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ /%%\%%%\%
%\%%%\%%  `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   -`-'%%\%%%\%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%


%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %


%This function is to create the eigenfunction bases for p and v, 
% in the case of an isentropic 2-D disc, D = B_r^2(0), with boundary
% conditions u * n = 0, dp/dn = 0, where the theta-Index of the chosen
% k-mode to be fixed is 0 -> complete angular symmetry -> radial
% eigenfunction bases.


%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %




%The eigenfunctions are described by the Helmholtz equation, and thus 
% the bases for p and v would actually be infinite-dimensional. Since the 
% problem is 2-D, the eigenfunctions are indexed by 2 variables, an r-index
% and a theta-index. Normally the theta-indicies are given by
%   (pertThetaIndex-1)*[0:1:thetaBasisSize-1]
% but since our index is 1 (which corresponds to value 0) we only have the
% single r-index, m=1,...,M

rEigenFncIndicies = 1:1:rBasisSize;

%The eigenfunctions are of the form phi_nm = R_nm(r), where the radial
%   component is given by R_nm(r) = J_m(sigmaBar*lambda_nm*r) where
%   lambda_nm is such that sigmaBar*lambda_nm*radius is the location at
%   which the derivative of the mth order Bessel function, J_m', reaches
%   zero for the nth time.

%Tolerance with which to find the zero derivative locations
besselFncZeroDerivTolerance = 1e-12;

%initialize matrix for the (n,m)th eigenvalues
eigenFreqMatrix = zeros(rBasisSize,1);
rEigenFunctionMatrix = zeros(rPoints,rBasisSize);




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% FIND RADIAL EIGENFUNCTIONS %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load precalculated data table containing the location of zeros for the
% derivatives of Bessel functions that were done to machine precision.
temp1 = load('besselFunctionPrecalculatedDerivativeZeros.mat');
temp2 = fieldnames(temp1);
besselDerivativeZeroDataTable = temp1.(temp2{1});
clear temp1 temp2

JZeroDerivLocations = besselDerivativeZeroDataTable(1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m = 1:rBasisSize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    eigenFreqMatrix(m) = JZeroDerivLocations(m) / (sigmaBar*domainRadius);

    %This matrix consists of the values of the m-th Bessel function 
    % corresponding to O_n(theta). This means that for O_n(theta) =
    % a_n*cos(n*theta) + b_n*sin(n*theta), the m-th corresponding r
    % eigenfunction is J_(n,m)(r) = J_n(sigmaBar*lambda_(n,m)*r) where
    % J_n denotes the n-th order Bessel function and lambda_(n,m) is such
    % that sigmaBar*lambda_(n,m)*R is the m-th derivative zero. This is 
    % so we have the Neumann cond satisified,
    % J'_n(sigmaBar*lambda*R)=0, where R is the radius of the disc.

    rEigenFunctionMatrix(:,m) = besselj(0,rGridBase * eigenFreqMatrix(m) * sigmaBar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% NORMALIZE  EIGENBASIS MATRIX %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The pEigenBasis matrix and the others will have dimensions
% N x M x thetaPoints x rPoints
% the first 2 indicies specify the (n,m) eigenfunction we want the
% matrix of values for

%So pEigenBasis(2,3,:,:) gives us the matrix odf values of the appropriate 
% version of the p-basis eigenfunction phi_(2,3)/sigmaBar on the
% discretization of the disk rGrid x thetaGrid

pEigenBasis = zeros(rPoints,rBasisSize);
pEigenBasis = rEigenFunctionMatrix;

%J_0(0) = 1, J'_0(0) = 0 while J_n(0) = 0, J'_n(0) ~= 0 for all n >= 1, 
%   so we will have one family of radially constant eigenfunctions
%   one of these also has constant angular part, so we treat the single
%   constant eigenfunction separately. Normalized, this eigenfunction 
%   will just be 1/|D|
pEigenFncZero = ones(rPoints,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m1 = 1:rBasisSize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pEigenFncTemp = zeros(rPoints,1);
    pEigenFncTemp(:) = pEigenBasis(:,m1);

    pEigenFncTempNorm = simps2DCopy(rGridBase,((pEigenFncTemp.*pEigenFncTemp).*rGridBase));

    pEigenBasis(:,m1) = pEigenBasis(:,m1)/sqrt(pEigenFncTempNorm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pEigenFncZeroNorm = simps2DCopy(rGridBase,((pEigenFncZero.*pEigenFncZero).*rGridBase));
pEigenFncZero = pEigenFncZero/sqrt(pEigenFncZeroNorm);

pBasisWeight = sigmaBar^2;
vBasisWeight = sigmaBar^((-1)*2);

pEigenBasis = pEigenBasis * (sigmaBar^(-1));
vEigenBasis = pEigenBasis * (sigmaBar^2);

pEigenFncZero = pEigenFncZero/sigmaBar;
vEigenFncZero = pEigenFncZero*sigmaBar^2;

%The u-basis eigenfunctions will be 2-D matricies, however a matrix of the
% matricies of values of these eigenfunctions is not actually required for
% the NL evolution. So for now the u-basis will not be calculated




end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


