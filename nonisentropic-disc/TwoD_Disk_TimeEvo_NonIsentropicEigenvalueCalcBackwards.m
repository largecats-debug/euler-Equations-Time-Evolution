
function [lambda,gamma,beta] = TwoD_Disk_TimeEvo_NonIsentropicEigenvalueCalcBackwards(rPoints,rGrid, ...
    radialEntropyLevels,entropyLevelWidths,entropyLevelStartingIndex,sigmaBarProfile,thetaEigenFncIndicies, ...
    rEigenFncIndicies,besselZeroDataTable,besselDerivativeZeroDataTable)


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

% METHOD:
%   Looking at this eigenvalue problem on each subinterval of constant
%   entropy (i.e. on each [r_(i-1),r_i] where i = 1,...,N) we take
%   R_nm(r) = beta_i J_n( lambda_nm * sigmaBar * (r-gamma_i) ) on each
%   subinterval [r_(i-1),r_i], i = 1,...,N.
%   
%   We immediately have gamma_1 = 0 so that the eigenfunctions will be
%   continuous at the origin, leaving us with 2N unknowns, gamma_i, beta_i,
%   lambda_nm
%
%   We get N equations each from imposing continuity of p and u at the
%   entropy jumps, [p]|_(r_i) = [u]|_(r_i) = 0
%
%   [p] = 0 
% -> beta_i J_n( lambda_nm * sigmaBar_i * (r_i - gamma_i) )
%    = beta_(i+1) J_n( lambda_nm * sigmaBar_(i+1) * (r_i - gamma_(i+1)) )
%   
%   [u] = 0 
% -> beta_i J_n( lambda_nm * sigmaBar_i * r_i - gamma_i )
%    lambda_nm sigmaBar_i beta_i d_r [J_n(lambda_nm*sigmaBar_i*(x_i - gamma_i))]
%     = lambda_nm sigmaBar_(i+1) beta_(i+1) d_r [ J_n(lambda_nm*sigmaBar_(i+1)*(r_i - gamma_(i+1) )) ]
%   
%    Divide the [u] equation by the [p] equation:
% -> sigmaBar_i d_r[ log( J_n(lambda_nm sigmaBar_i (r_i - gamma_i)) ) ]
%     = sigmaBar_(i+1) d_r[ log( J_n(lambda_nm sigmaBar_(i+1) r_i - gamma_(i+1)) ) ]
%
%   For each i = 1,...,N, the above equation can be solved for
% ->gamma_(i+1) = gamma_(i+1)(gamma_i,lambda_nm) = r_(i+1) - 
%   - (1/(lambda_nm sigmaBar_(i+1)) * d_r[J_n(lambda_nm sigmaBar_N (r_N - gamma_N(lambda_nm)) )]^(-1) [(sigmaBar_(i+1)/sigmaBar_i) d_r[J_n(lambda_nm sigmaBar_N (r_N - gamma_N(lambda_nm)) )][lambda_nm sigmaBar_i (r_i-gamma_i)] ]
%   However, Bessel functions (in particular, we are using Bessel functions
%   of the first kind) are oscillatory, so we need to make use of the
%   function "unwrap" when defining each gamma_i(lambda_nm)
%
%   The condition that u(T,x) = 0 translates to the condition that
%   J_n'( lambda_nm sigmaBar_N (x_N - gamma_N) ) = 0, so each lambda_nm by
%   solving lambda_nm * sigmaBar_N * (x_N - gamma_N) = j_nm where j_nm is
%   the mth zero of the derivative of J_n(r)
%
%
%
%   For each fixed theta-index n, iterate m and find lambda_nm by:
%   1. Fix n
%   2. Fix m
%   3. Find {j_nm}, the mth zero of J_n'(r)
%   4. derive gamma_N = gamma_N(lambda_nm) by nesting anonymous functions
%   5. Solve lambda_nm * sigmaBar_N * (x_N - gamma_N(lambda_nm)) = j_nm
%       for each m = 1,...,rBasisSize
%   6. iterate over m
%   7. iterate over n


%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%    function intervalIndex = intervalOfAngle(angle,n)%%
%
%        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        if angle < besselZeroDataTable(n,1)%%
%
%            intervalIndex = 1;%
%
%            return
%
%       end
%        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        for j1 = 1:length(besselZeroDataTable(n,:))-1
%
%            intervalIndex = j1+1;
%
%            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%            if angle > besselZeroDataTable(n,j1) && angle < besselZeroDataTable(n,j1+1)
%
%                return
%
%            end
%            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
%        end
%        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%    function inverseValue = besselLogDerivativeInverse(nIndex,value,intervalIndex,besselLogDerivativeFncHandle)
%
%        fAbs = @(y) abs(besselLogDerivativeFncHandle(y)-value);
%
%        if intervalIndex == 1
%
%            inverseValue = fminbnd(fAbs, 0, besselZeroDataTable(nIndex,intervalIndex)-0.00001);
%
%        
%        else
%
%            inverseValue = fminbnd(fAbs, besselZeroDataTable(nIndex,intervalIndex-1)+0.00001, besselZeroDataTable(nIndex,intervalIndex)-0.00001);
%
%        end
%
%    end

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

thetaBasisSize = length(thetaEigenFncIndicies);
rBasisSize = length(rEigenFncIndicies);

lambda = zeros(thetaBasisSize,rBasisSize);
gamma = zeros(thetaBasisSize,rBasisSize,radialEntropyLevels);
beta = ones(thetaBasisSize,rBasisSize,radialEntropyLevels);
besselLogDerivativeFnc = @(q) 0;
options = optimoptions('fsolve','FunctionTolerance',1e-6,'StepTolerance',1e-10,'OptimalityTolerance',1e-9,'MaxFunctionEvaluations',1e7);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n1 = 1:thetaBasisSize

    n1Value = thetaEigenFncIndicies(n1);

    clear besselLogDerivativeFnc
    clear fLoop fLoopInverse fLoopIndex

    fLoop = cell(radialEntropyLevels,1);
    fLoopInverse = cell(radialEntropyLevels,1);
    fLoopIndex = cell(radialEntropyLevels,1);


    besselLogDerivativeFnc = @(q) ((1/2)*(besselj(n1Value-1,q)-besselj(n1Value+1,q)))/(besselj(n1Value,q));
    %fLoop{1} = @(x) besselLogDerivativeFnc(x*sigmaBarProfile(1)*(rGrid(entropyLevelStartingIndex(2))-0));


    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    for m1 = 1:rBasisSize

        %lambda(n1,m1) = besselLogDerivativeFnc(besselDerivativeZeroDataTable(n1,m1)/(sigmaBarProfile(radialEntropyLevels)*rGrid(rPoints)));
        %angleTemp = lambda(n1,m1)*sigmaBarProfile(radialEntropyLevels)*rGrid(entropyLevelStartingIndex(i1));
        %angleTempIndex = intervalOfAngle(angleTemp,n1,besselZeroDataTable);
        %angleTempInverse = besselLogDerivativeInverse(n1,(sigmaBarProfile(i1+1)/sigmaBarProfile(i1))*angleTemp,angleTempIndex,besselLogDerivativeFnc,besselZeroDataTable);
        lambda(n1,m1) = besselDerivativeZeroDataTable(n1Value+1,m1)/(rGrid(rPoints)*sigmaBarProfile(radialEntropyLevels))
        %gammaFunctions{n1,1} = @(x) 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i1 = radialEntropyLevels-1:-1:1

            angleTemp = lambda(n1,m1)*sigmaBarProfile(radialEntropyLevels)*(rGrid(entropyLevelStartingIndex(i1))-gamma(n1,m1,i1+1));
            angleTempIndex = intervalOfAngle(angleTemp,n1Value+1,besselZeroDataTable);
            angleTempInverse = besselLogDerivativeInverse(n1Value+1,(sigmaBarProfile(i1+1)/sigmaBarProfile(i1))*angleTemp,angleTempIndex,besselLogDerivativeFnc,besselZeroDataTable);
            gamma(n1,m1,i1) = rGrid(entropyLevelStartingIndex(i1)) - ((1/(lambda(n1,m1)*sigmaBarProfile(i1)))*angleTempInverse);
            beta(n1,m1,i1) = (beta(n1,m1,i1+1)*besselj(n1Value,lambda(n1,m1)*sigmaBarProfile(i1+1)*(rGrid(entropyLevelStartingIndex(i1+1))-gamma(n1,m1,i1+1))))/(besselj(n1Value,lambda(n1,m1)*sigmaBarProfile(i1)*(rGrid(entropyLevelStartingIndex(i1+1))-gamma(n1,m1,i1))));

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %fLoop = @(x) abs(  (1/2) * ( besselj( n-1,x*sigmaBarProfile(radialEntropyLevels)*(rGrid(rPoints)-fLoop(x)) ) - besselj( n+1,x*sigmaBarProfile(radialEntropyLevels)*(rGrid(rPoints)-fLoop(x)) ) )  );
        %finalAngleDiff = @(x) x * sigmaBarProfile(radialEntropyLevels) * (rGrid(rPoints)-gammaFunctions{n1,radialEntropyLevels}(x)) - besselDerivativeZeroDataTable(n1,m1);

%if n1 == 1 && m1 == 1
%   fLoopReturn = @(x) fLoop(x);
%end

        %lambda(n1,m1) = fminbnd( fLoop , (besselDerivativeZeroDataTable(n1,m1)/(sigmaBarProfile(radialEntropyLevels)*rGrid(rPoints)))*0.8 , (besselDerivativeZeroDataTable(n1,m1)/(sigmaBarProfile(radialEntropyLevels)*rGrid(rPoints)))*1.2 );
        %lambda(n1,m1) = fminbnd( fLoop , (besselDerivativeZeroDataTable(n1,m1)/(sigmaBarProfile(1)*entropyLevelWidths(1)+sigmaBarProfile(2)*entropyLevelWidths(2)))*0.8 , (besselDerivativeZeroDataTable(n1,m1)/(sigmaBarProfile(1)*entropyLevelWidths(1)+sigmaBarProfile(2)*entropyLevelWidths(2)))*1.2 );
        %if m1 == 1
            %lambda(n1,m1) = fsolve(finalAngleDiff,(besselDerivativeZeroDataTable(n1,m1)/(sigmaBarProfile(2)*rGrid(rPoints))),options);
        %else
        %    lambda(n1,m1) = fsolve(fLoop,lambda(n1,1)*m1,options);
        %end

        %gammaLoop = besselLogDerivativeFnc(lambda(n1,m1)*sigmaBarProfile(1)*(rGrid(entropyLevelStartingIndex(2))-0));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %for i2 = 2:radialEntropyLevels
        %    gammaLoop = besselLogDerivativeFnc(lambda(n1,m1)*sigmaBarProfile(i2-1)*(rGrid(entropyLevelStartingIndex(i2))-gamma(n1,m1,i2-1)));
        %    gammaLoopIndex = intervalOfAngle(gammaLoop,n1,besselZeroDataTable);
        %    gammaLoop = besselLogDerivativeInverse(n1,gammaLoop,gammaLoopIndex,besselLogDerivativeFnc,besselZeroDataTable);
        %    gamma(n1,m1,i2) = rGrid(entropyLevelStartingIndex(i2)) - (1/(lambda(n1,m1)*sigmaBarProfile(i2))) * gammaLoop;
        %    beta(n1,m1,i2) = (beta(n1,m1,i2-1)*besselj(n1Value,lambda(n1,m1)*sigmaBarProfile(i2-1)*(rGrid(entropyLevelStartingIndex(i2))-gamma(n1,m1,i2-1))))/(besselj(n1Value,lambda(n1,m1)*sigmaBarProfile(i2)*(rGrid(entropyLevelStartingIndex(i2))-gamma(n1,m1,i2))));
        %end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%