    %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %


%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%
%%%\%%%"-.___ ("""-.___ ("""-.___ ("""-.___ ("""-.___ ("""-.___ ("""-.___ ("""-\%%%\%%
%%%\%%%)   __) "")   __) "")   __) "")   __) "")   __) "")   __) "")   __) "") \%%%\%%
%%%\%%%   (__,-""   (__,-""   (__,-""   (__,-""   (__,-""   (__,-""   (__,-""  \%%%\%%
%%%\%%%"-.___ %\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\% ("""-\%%%\%%
%%%\%%%)   __)%\%\ 2-D TIME EVOLUTION - DISC - NON ISENTROPIC - MAIN \%\%) "") \%%%\%%
%%%\%%%   (__,%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%,-""  \%%%\%%
%%%\%%%"-.___ ("""-.___ ("""-.___ ("""-.___ ("""-.___ ("""-.___ ("""-.___ ("""-\%%%\%%
%%%\%%%)   __) "")   __) "")   __) "")   __) "")   __) "")   __) "")   __) "") \%%%\%%
%%%\%%%   (__,-""   (__,-""   (__,-""   (__,-""   (__,-""   (__,-""   (__,-""  \%%%\%%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%
%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%Set constants for the isentropic disc in question (currently hardcoded to
%   best approximate air at sea level (at some specified temperature)
domainRadius = 1;
domainSize = pi*domainRadius*domainRadius;


%Set the number of eigenfunctions that will be used to approximate the
%   infinite dimensional bases of eigenfunctions described by our EVP.
%   Since we are working with a 2-D spatial domain, the eigenfunctions are 
%   indexed by 2 variables, phi_nm = O_n(theta)*R_nm(r)
%   So, specificy size of basis from number of r-indicies and theta-indicies
rBasisSize = 8;
thetaBasisSize=8;
%NOTE: the eigenfunctions are actually indexed by 3 variables. This is
%   because for each n=1,...,N, O_n(theta) can be both cos((n-1)*theta) and
%   sin((n-1)*theta). If including both variants we would index the
%   eigenfunction basis by phi_nml where again n=1,...,N and m=1,...,M, plus
%   the index l=0,1 where l=0 is for cos and l=1 is for sin.

%If the eigenfunction to be perturbed is phi_{k,k',l}, then we can only
%   make use of an eigenfunction basis that is a subset of the form
%   { phi_{n*k,m,0} | n=0,...,N-1 , m=1,...,M } if l=0
%   { phi_{n*k,m,mod(n,2) | n=0,...,N-1 , m=1,...,M} if l=1
%   so our basis of eigenfunctions will not need a third index.

%Set the eigenfunction that we will perturb to find a solution of the
%   compressible Euler equations.
%The p-basis consists of functions phi_{n,m,l} where:
%   n = theta-index (number of angular symmetries + 1, cos((n-1)Theta)/sin((n-1)Theta)
%   m = r-index (number of radial symmetries)
%   l = 0 or 1 (picks whether the theta eigenfunction is cos (1) or sin (2)
pertThetaIndex = 3;
pertRIndex = 2;
pertTrigIndex = 1;

%because matlab indexing starts at 1 instead of 0, having for example
%   thetaIndex=3, rIndex=3, trigIndex = 1 corresponds to
%   cos(2*theta)*J_2(lambda_(2,3)*sigmaBar*r)

%Set number of points for 2D grid on disk by number of radial points and
%   number of angular points, "r-points" and "theta-points"
rPoints = 500;
thetaPoints = 300;
rGrid = linspace(0,domainRadius,rPoints);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

pBar = 101325;
pConstIntMatrix = pBar*ones(thetaPoints,rPoints);

%If our fixed k-mode is theta-index zero eigenfunction 
%   (i.e. phi_k = phi_k*(r)) then we only need to go through the whole
%   computational process for a constant angle slice of the domain.
%   (i.e. either a disk or annulus) So, we check this, and if this is the
%   case then we set thetaBasisSize = 1 and prompt the user on if they
%   would like to set a new rBasisSize
%       NOTE 9/28/25: Option to perform calculation on annulus still being
%       worked on
if pertThetaIndex == 1

    thetaBasisSize = 1;

    clear prompt
    prompt = ['The theta-index of the chosen k-mode was 0, so thetaBasisSize has been set to 1', newline, 'Set a new (or the same, if still desired) value for rBasisSize: '];
    rBasisSize = input([prompt,newline]);

    prompt = ['This means that after the computations, the resolution of the [0,2pi] discretization', newline, 'that surfaces are plotted on after the computations will be thetaPoints.', newline, 'Set a new (or the same, if still desired) value for thetaPoints: '];
    thetaPoints = input([prompt,newline]);

    %We still define a thetaGrid so that after the computations are run, we
    %   have a mesh covering the whole domain on which to recompose and plot
    %   the matricies we have the values of on a single ray
    thetaGridFull = linspace(0,2*pi,thetaPoints);

else

    thetaGrid = linspace(0,2*pi/(pertThetaIndex-1),thetaPoints);
    thetaGridFull = linspace(0,2*pi,(pertThetaIndex-1)*thetaPoints);

end


%Set size of time grid on [0,T], i.e. the number of finite difference steps
%   that will be taken to calculate the nonlinear evolution
timePoints = 200;


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\%%%\
%entropy levels are determined by specifying temperatures for each level,
% where the temperature in the first region is taken as the reference
% temperature. Prompt user for these n temperature values (degrees F)
disp("The different entropy levels are calculated by specifying a temperature for each region and assuming the pressure is a constant 1 atm (101325 Pa)" ...
    + newline + "The entropy levels are calculated with respect to the base temperature in the first region")
%We use a data table for the density of air at 1 atm of pressure at various
% temperatures to estimate the density in each region of a different
% temperature. Invert these for the specific volume in each region.
% The change in entropy from the reference state in the first region is
% then given by:
% s_2 - s_1 = c_P log(T_2/T_1)

radialTemperatureProfile = zeros(1,radialEntropyLevels);
radialAirDensityProfile = zeros(1,radialEntropyLevels);
radialSpecificVolumeProfile = zeros(1,radialEntropyLevels);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:radialEntropyLevels
    clear prompt
    prompt = ['Specify the temperature (in degrees Fahrenheit) for region ', int2str(i), '(between 50 and 80 F): '];
    radialTemperatureProfile(i) = input([prompt,newline])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

backgroundTemperature = radialTemperatureProfile(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seaLevelAirCp = [1003.84,1003.87,1003.89,1003.92,1003.94,1003.97,1004.00,1004.02,1004.05,1004.08,1004.11,1004.13,1004.16,1004.19,1004.22,1004.25,1004.28,1004.30,1004.33,1004.36,1004.39,1004.42,1004.45,1004.48,1004.51,1004.55,1004.58,1004.61,1004.64,1004.67,1004.70];
seaLevelAirCv = [716.79,716.82,716.84,716.87,716.89,716.92,716.95,716.97,717.00,717.03,717.06,717.08,717.11,717.14,717.17,717.20,717.23,717.25,717.28,717.31,717.34,717.37,717.40,717.43,717.46,717.50,717.53,717.56,717.59,717.62,717.65];
seaLevelAirGamma = seaLevelAirCp./seaLevelAirCv;
seaLevelAirDensity = [1.246,1.244,1.241,1.239,1.237,1.234,1.232,1.229,1.227,1.225,1.222,1.22,1.218,1.215,1.213,1.211,1.208,1.206,1.204,1.202,1.199,1.197,1.195,1.193,1.19,1.188,1.186,1.184,1.182,1.179,1.177];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

radialTemperatureLowerIndex = floor(radialTemperatureProfile)-49;
radialTemperatureDecimal = mod(radialTemperatureProfile,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:radialEntropyLevels

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    if radialTemperatureDecimal(i) ~= 0
        radialAirDensityProfile(i) = (1-radialTemperatureDecimal(i))*seaLevelAirDensity(radialTemperatureLowerIndex(i)) + ...
            radialTemperatureDecimal(i)*seaLevelAirDensity(radialTemperatureLowerIndex(i)+1);
        radialCpProfile(i) = (1-radialTemperatureDecimal(i))*seaLevelAirCp(radialTemperatureLowerIndex(i)) + ...
            radialTemperatureDecimal(i)*seaLevelAirCp(radialTemperatureLowerIndex(i)+1);
        radialCvProfile(i) = (1-radialTemperatureDecimal(i))*seaLevelAirCv(radialTemperatureLowerIndex(i)) + ...
            radialTemperatureDecimal(i)*seaLevelAirCv(radialTemperatureLowerIndex(i)+1);
        radialRProfile(i) = radialCpProfile(i) - radialCvProfile(i);
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    if radialTemperatureDecimal(i) == 0
        radialAirDensityProfile(i) = seaLevelAirDensity(radialTemperatureLowerIndex(i));
        radialCpProfile(i) = seaLevelAirCp(radialTemperatureLowerIndex(i));
        radialCvProfile(i) = seaLevelAirCv(radialTemperatureLowerIndex(i));
        radialRProfile(i) = radialCpProfile(i) - radialCvProfile(i);
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    radialGammaProfile(i) = (radialCpProfile(i)/radialCvProfile(i));
    radialSpecificVolumeProfile(i) = 1/radialAirDensityProfile(i);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vBar = radialSpecificVolumeProfile(1);
radialEntropyProfile(1) = 0;
entropyWidthsSum = 0;
entropyLevelWidths = zeros(radialEntropyLevels,1);
sigmaBarProfile = zeros(radialEntropyLevels,1);

clear prompt
prompt = 'Specify the length for entropy level 1 (in meters): ';
entropyLevelWidths(1) = input([prompt,newline]);
entropyWidthsSum = entropyWidthsSum + entropyLevelWidths(1);
sigmaBarProfile(1) = (sqrt(vBar)/(sqrt(radialGammaProfile(1))*sqrt(pBar)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 2:radialEntropyLevels

    radialEntropyProfile(i) = ((radialCpProfile(1)+radialCpProfile(i))/2) * ...
        log(radialTemperatureProfile(i)/radialTemperatureProfile(1)) + ...
        ((radialRProfile(i) + radialRProfile(1))/2) * ...
        log(radialSpecificVolumeProfile(i)/radialSpecificVolumeProfile(1));

    sigmaBarProfile(i) = (sqrt(vBar)/(sqrt(radialGammaProfile(1))*sqrt(pBar))) * exp((-1)*(radialEntropyProfile(i)/radialCpProfile(1)));

    clear prompt
    prompt = ['Specify the length for entropy level ', int2str(i), ' (in meters): '];
    entropyLevelWidths(i) = input([prompt,newline]);
    entropyWidthsSum = entropyWidthsSum + entropyLevelWidths(i);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

entropyLevelStartingIndex = zeros(radialEntropyLevels,1);
entropyLevelStartingIndex(1) = 1;
entropyIntervalPoints = zeros(radialEntropyLevels,1);

entropyWidthsSumLoop = entropyWidthsSum;
rPointsLoop = rPoints;
floorLoop = floor((entropyLevelWidths(1)/entropyWidthsSum)*rPointsLoop);
entropyIntervalPoints(1) = floorLoop;
rPointsLoop = rPointsLoop - entropyIntervalPoints(1);
entropyWidthsSumLoop = entropyWidthsSumLoop - entropyLevelWidths(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i1 = 2:radialEntropyLevels

    entropyLevelStartingIndex(i1) = entropyLevelStartingIndex(i1-1) + entropyIntervalPoints(i1-1);

    floorLoop = floor((entropyLevelWidths(i1)/entropyWidthsSumLoop)*rPointsLoop);
    entropyIntervalPoints(i1) = floorLoop;
    rPointsLoop = rPointsLoop - entropyIntervalPoints(i1);
    entropyWidthsSumLoop = entropyWidthsSumLoop - entropyLevelWidths(i1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%for j1 = 2:radialEntropyLevels
%    floorRatioLoop = floor(entropyLevelWidths(j1)/entropyWidthsSum);
%    entropyIntervalPoints(j1) = entropyIntervalPoints(j1) + floorRatioLoop;
%    if mod(floorRatioLoop*rPoints,1) ~= 0
%        entropyIntervalPoints(j1+1) = 1;
%    end
%    for j2 = 1:j1-1
%        entropyLevelStartingIndex(j1) = entropyLevelStartingIndex(j1) + entropyIntervalPoints(j2);
%    end
%    entropyLevelStartingIndex(j1) = entropyLevelStartingIndex(j1) + 1;
%end




%Now generate the eigenfunction basis for p and v
[pEigenBasis,pEigenFncZero,pBasisWeight,vEigenBasis,vEigenFncZero,vBasisWeight,JScaledValues,...
    eigenFreqMatrix,JZeroDerivLocations,JScaledValuesNorm] = ...
    TwoD_Disk_TimeEvo_NonIsentropicBasis(rBasisSize,thetaBasisSize,sigmaBarConst,domainRadius,rGrid, ...
    rPoints,thetaGrid,thetaPoints,pertThetaIndex,pertTrigIndex,radialEntropyLevels,entropyLevelStartingIndex);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Again since we only use theta-indicies that are integer multiples of the
%   pertThetaIndex, i.e. n*[0:N], to separately save the eigenfunction to be
%   perturbed, we want the eigenfunction on the 2nd row and pertRIndex-th 
%   column of the basis matrix
pEigenFncPertTemp(:,:) = pEigenBasis(2,pertRIndex,:,:);
































