function eigenFreqMatrix = TwoD_Disk_TimeEvo_NonIsentropicEigenvalueCalc2(pertThetaIndex, ...
    pertRIndex,rPoints,rGrid,radialEntropyLevels,entropyLevelStartingIndex)


%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%
%\%%%\%%.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   %%\%%%\%
%\%%%\%% \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ /%%\%%%\%
%\%%%\%%  `-/-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   -`-'%%\%%%\%
%\%%%\%%-.___ %\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\% ("""%%\%%%\%
%\%%%\%%   __)%\%\    2-D TIME EVOLUTION - DISC - NON-ISENTROPIC     \%\%) "")%%\%%%\%
%\%%%\%%  (__)%\%\              EIGENVALUE CALCULATION               \%\%)("")%%\%%%\%
%\%%%\%%  (__,%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%,-"" %%\%%%\%
%\%%%\%%.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   %%\%%%\%
%\%%%\%% \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ /%%\%%%\%
%\%%%\%%  `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   -`-'%%\%%%\%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%

%This function calculates a specific index eigenvalue lambda_(n,m) for the
%   case of a disc with radially piecewise constant entropy,
%   s(r) = s_i   for   r_(i-1) < r < r_i   where   r_0 = 0, r_(N-1) = R
%   for the eigenfunction-eigenvalue problem
%   nabla^2[phi_nm/sigmaBar] = - lambda^2_(n,m) sigmaBar^2 (phi_nm/sigmaBar)

%%  METHOD:
%   Start with an


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The eigenfunctions are essentially described by the Helmholtz equation,
%   and thus the p and v bases would actually be infinite-dimensional. Since
%   the problem is 2-D, the eigenfunctions are indexed by 2 variables, an
%   r-index and a theta-index.
%So the finite eigenfunction bases for p and v are determined by two
%parameters N for the theta-index n and M for the r-index M:
%n = 0,...,N-1 , m=1,...,M

thetaEigenFncIndicies = (pertThetaIndex-1)*(0:1:thetaBasisSize-1);

rEigenFncIndicies = 1:1:rBasisSize;

%The eigenfunctions are of the form phi_nm = O_n(theta)*R_nm(r). The radial
%   component is given by R_nm(r) = J_m(sigmaBar*lambda_nm*r) where
%   lambda_nm is such that sigmaBar*lambda_nm*radius is the location at
%   which the derivative of the mth order Bessel function, J_m', reaches
%   zero for the nth time.

%Tolerance with which to find the zero derivative locations
besselFncZeroDerivTolerance = 1e-12;

%initialize matrix for the (n,m)th eigenvalues
eigenFreqMatrix = zeros(thetaBasisSize,rBasisSize);
JScaledValues = zeros(thetaBasisSize,rBasisSize,rPoints);
JScaledValuesNorm = zeros(thetaBasisSize,rBasisSize);
rEigenFunctionMatrix = zeros(thetaBasisSize,rBasisSize,rPoints);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Unstored second variable is JZeroDerivValues
[JZeroDerivLocations,~] = BessDerivZerosBisect2(thetaEigenFncIndicies,rEigenFncIndicies,besselFncZeroDerivTolerance);

%Now, with an ordered lists of the locations of zero derivatives of J_n(x),
% R'_(n,m)(r_0) = J'_n(lambda_(n,m)*sigmaBar*r_0) -> r_m = lambda_(n,m)*sigmaBar*r_0


for n = 1:thetaBasisSize

    for m = 1:rBasisSize

        

        eigenFreqMatrix(n,m) = JZeroDerivLocations(n,m) / (sigmaBar*domainRadius);

        %This matrix consists of the values of the m-th Bessel function 
        % corresponding to O_n(theta). This means that for O_n(theta) =
        % a_n*cos(n*theta) + b_n*sin(n*theta), the m-th corresponding r
        % eigenfunction is J_(n,m)(r) = J_n(sigmaBar*lambda_(n,m)*r) where
        % J_n denotes the n-th order Bessel function and lambda_(n,m) is such
        % that sigmaBar*lambda_(n,m)*R is the m-th derivative zero. This is 
        % so we have the Neumann cond satisified,
        % J'_n(sigmaBar*lambda*R)=0, where R is the radius of the disc.

        JScaledValues(n,m,:) = besselj(thetaEigenFncIndicies(n),rGridBase * eigenFreqMatrix(n,m) * sigmaBar);

        JScaledVectorLoop(:) = JScaledValues(n,m,:);
        %JScaledVectorLoop(:) = JScaledValues(n,m+1,:);

        JScaledValuesNorm(n,m) = sqrt(trapz(rGridBase,(JScaledVectorLoop.*JScaledVectorLoop).*rGridBase));
        rEigenFunctionMatrix(n,m,:) = JScaledValues(n,m,:)/JScaledValuesNorm(n,m);

    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We now have an array of values for each R_(n,m) on rGrid. So we just need
%   to similarly get arrays of values of each O_n on thetaGrid to take the
%   product and find the full p-basis (and v-basis).

%Solving the theta equation when separating variables just gives
%   O_n(theta) = a*cos(n*theta) + b*sin(n*theta)

%NOTE 9/16/25:
%Currently disregarding sine variant of eigenfunctions because depending on
%   whether the starting k-mode choice of an eigenfunction is a cosine or
%   sine variant. The basis of eigenfunctions needed will be those that we
%   can get by interacting the $k$-mode with itself (any number of times)

%E.g. if we are starting with an eigenfunction with angular term
%   cos(2w), then we will need to construct the p-eigenbasis as all
%   eigenfunctions that have angular part 1,cos(2w),cos(4w),cos(6w),...
%   The list of possible theta-indicies start at n=0, however matlab indexing
%   starts at 1, so when we have theta index n, this really corresponds to
%   (n-1) angular symmetries, that is cos((n-1)w) or sin(n-1)w)

thetaEigenFunctionMatrix = zeros(thetaBasisSize,thetaPoints);
 
thetaEigenFunctionMatrix(1,:) = (1/sqrt(2*pi))*cos((0)*ones(1,thetaPoints));

    if pertThetaIndex ~= 1

        if pertTrigIndex == 1

            %Scaling factor to normalize eigenfunctions different for
            %   theta-index = 0 eigenfunctions. So set these (would be n = 1) and
            %   then just loop starting from n = 2

            for n = 2:thetaBasisSize

                if pertTrigIndex == 1

                    thetaEigenFunctionMatrix(n,:) = (1/sqrt(pi))*cos((n-1)*(pertThetaIndex-1)*thetaGrid);

                end

            end

        %If we want our fixed k-mode phi_k = R_nm(r)O_m(phi) to be a sine
        %   variant eigenfunction O_m(phi) = cos(m*phi) or sin(m*phi)
        %   then the iteraction terms will have terms of the form
        %   cos(j*m*phi) or sin(j*m*phi) depending on whether j is even
        %   or odd [even j -> sin^2 * ... * sin^2 -> cos(2phi), cos(0)
        %       (NOTE: combining any of these interaction terms into new
        %       iteraction terms does not create any new trig terms)
        else


            if mod(n,2) == 0
            
                thetaEigenFunctionMatrix(n,:) = (1/sqrt(pi))*cos((n-1)*(pertThetaIndex-1)*thetaGrid);

            end

            if mod(n,2) == 1

                thetaEigenFunctionMatrix(n,:) = (1/sqrt(pi))*sin((n-1)*(pertThetaIndex-1)*thetaGrid);

            end

      
        end

    end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The pEigenBasis matrix and the others will have dimensions
% N x M x thetaPoints x rPoints
% the first 2 indicies specify the (n,m) eigenfunction we want the
% matrix of values for

%So pEigenBasis(2,3,:,:) gives us the matrix of values of the appropriate 
% version of the p-basis eigenfunction phi_(2,3)/sigmaBar on the
% discretization of the disk rGrid x thetaGrid

if pertThetaIndex ~= 1

    pEigenBasis = zeros(thetaBasisSize,rBasisSize,thetaPoints,rPoints);

    for n = 1:thetaBasisSize

        for m = 1:rBasisSize

            for i = 1:thetaPoints

                pEigenBasis(n,m,i,:) = rEigenFunctionMatrix(n,m,:) * thetaEigenFunctionMatrix(n,i) / sigmaBar;

            end

        end

    end

    %J_0(0) = 1, J'_0(0) = 0 while J_n(0) = 0, J'_n(0) ~= 0 for all n >= 1, 
    %   so we will have one family of radially constant eigenfunctions
    %   one of these also has constant angular part, so we treat the single
    %   constant eigenfunction separately. Normalized, this eigenfunction 
    %   will just be 1/|D|
    pEigenFncZero = (1/(sqrt(pi)*domainRadius))*ones(thetaPoints,rPoints)*(1/sigmaBar);
    vEigenFncZero = (1/(sqrt(pi)*domainRadius))*ones(thetaPoints,rPoints)*sigmaBar;


else

    for m = 1:rBasisSize

        pEigenBasis(1,m,1,:) = rEigenFunctionMatrix(1,m,:) .* thetaEigenFunctionMatrix(1,:) / sigmaBar;

    end

    pEigenFncZero = (1/(sqrt(pi)*domainRadius))*ones(1,rPoints)*(1/sigmaBar);
    vEigenFncZero = (1/(sqrt(pi)*domainRadius))*ones(1,rPoints)*sigmaBar;

end


%the v-basis is just the p-basis scaled by sigmaBar^2. Since this is for
%the isentropic case, this wouldn't change anything after normalizing.
%However, for consistency, we will keep the bases as follows:
% pBasis = {phi_nm/sigmaBar}, orthogonal w.r.t. w_1(x) = sigmaBar^2(x)
% vBasis = {phi_nm*sigmaBar}, orthogonal w.r.t. w_2(x) = sigmaBar^(-2)(x)
vEigenBasis = sigmaBar^2 * pEigenBasis;
pBasisWeight = sigmaBar^2;
vBasisWeight = sigmaBar^(-2);

%The u-basis eigenfunctions will be 2-D matricies, however a matrix of the
% matricies of values of these eigenfunctions is not actually required for
% the NL evolution. So for now the u-basis will not be calculated




end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


