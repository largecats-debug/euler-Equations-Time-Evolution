%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%##\ .--. \##\ .-'. \##\ .--. \##\ .--. \##\ .--. \##\ .--. \##\ .`-. \##\.--%%%%
%%%%%:::::.\::::::::.\::::::::.\::::::::.\::::::::.\::::::::.\::::::::.\:::::::%%%%%
%%%%-' \##\ `--' \##\ `.-' \##\ `--' \##\ `--' \##\ `--' \##\ `-.' \##\ `--' \##%%%%
%%%%% \#\ /""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7 \#\ %%%%%
%%%%\  .--/""7__/""                                             /""7__/""7\ .--.%%%%
%%%%%:::::/""7__/" 2-D TIME EVOLUTION - DISC - ISENTROPIC - MAIN \"7__/""7:::::%%%%%
%%%%-'  \%/""7__/""                                             /""7__/""7-' \#\%%%%
%%%%% \#\ /""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7 \#\.%%%%%
%%%%##\ .--. \##\ .-'. \##\ .--. \##\ .--. \##\ .--. \##\ .--. \##\ .`-. \##\.--%%%%
%%%%%:::::.\::::::::.\::::::::.\::::::::.\::::::::.\::::::::.\::::::::.\:::::::%%%%%
%%%%-' \##\ `--' \##\ `.-' \##\ `--' \##\ `--' \##\ `--' \##\ `-.' \##\ `--' \##%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


%Set constants for the isentropic disc in question (currently hardcoded to
%   best approximate air at sea level (at some specified temperature)
domainRadius = 1;
domainSize = pi*domainRadius*domainRadius;


%Set number of points for 2D grid on disk by number of radial points and
% number of angular points, "r-points" and "theta-points"
rPoints = 500;
thetaPoints = 300;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
useUniformRGrid = 2;

while useUniformRGrid ~= 0 && useUniformRGrid ~= 1

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    clear prompt
    prompt = [newline,'Would you like to use a uniform discretization for the discretization the radial interval [0,R]',newline,'(Input 1 for yes, 0 for no. If no, you will be prompted to make a choice regarding a nonuniform discretization)'];
    useUniformRGrid = input([prompt,newline]);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

end

if useUniformRGrid

    rGrid = linspace(0,domainRadius,rPoints);

else

    gridTransformExp = 2;

    while gridTransformExp < 0 || gridTransformExp > 1

   
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        clear prompt
        prompt = [newline,'Input an a value for the exponent a, where the uniform discretization linspace(0,domainRadius,rPoints) will be transformed according to f(x) = x^a', newline, '(So, choose 0 < a < 1 to cluster points closer to the edge of the disc)'];
        gridTransformExp = input([prompt,newline]);
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    end

    rGrid = linspace(0,domainRadius,rPoints);
    rGrid = rGrid.^gridTransformExp;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Set size of time grid on [0,T], i.e. the number of finite difference steps
%   that will be taken to calculate the nonlinear evolution
timePoints = 200;


%Set the number of eigenfunctions that will be used to approximate the
% infinite dimensional bases of eigenfunctions described by the our EVP.
% Since we are working with a 2-D spatial domain, the eigenfunctions are 
% indexed by 2 variables, phi_nm = O_n(theta)*R_nm(r)
% So, specificy size of basis from number of r-indicies and theta-indicies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if promptUserEigenBasisSizes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    clear prompt
    prompt = [newline,'Input the number of theta-indicies to use in constructing our finite eigenbasis', newline, '(reccomended thetaBasisSize is ~16 at the higher end)'];
    thetaBasisSize = input([prompt,newline]);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    clear prompt
    prompt = [newline,'Input the number of r-indicies to use in constructing our finite eigenbasis', newline, '(reccomended rBasisSize is ~12 at the higher end)'];
    rBasisSize = input([prompt,newline]);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rBasisSize = 10;
    thetaBasisSize=20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%NOTE: the eigenfunctions are actually indexed by 3 variables. This is
% because for each n=1,...,N, O_n(theta) can be both cos((n-1)*theta) and
% sin((n-1)*theta). If including both variants we would index the
% eigenfunction basis by phi_nml where again n=1,...,N and m=1,...,M, plus
% the index l=0,1 where l=0 is for cos and l=1 is for sin.

%If the eigenfunction to be perturbed is phi_{k,k',l}, then we can only
% make use of an eigenfunction basis that is a subset of the form
% { phi_{n*k,m,0} | n=0,...,N-1 , m=1,...,M } if l=0
% { phi_{n*k,m,mod(n,2) | n=0,...,N-1 , m=1,...,M} if l=1
% so our basis of eigenfunctions will not need a third index.


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


%Set the eigenfunction that we will perturb to find a solution of the
% compressible Euler equations.
%The p-basis consists of functions phi_{n,m,l} where:
% n = theta-index (number of angular symmetries + 1, cos((n-1)Theta)/sin((n-1)Theta)
% m = r-index (number of radial symmetries)
% l = 0 or 1 (picks whether the theta eigenfunction is cos (1) or sin (2)

%Because matlab indexing starts at 1 instead of 0, having for example
% thetaIndex=3, rIndex=3, trigIndex = 1 corresponds to
% cos(2*theta)*J_2(lambda_(2,3)*sigmaBar*r)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if promptUserKModeIndicies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    clear prompt
    prompt = [newline,'Input the desired theta-index for the fixed k-mode', newline, '(Valid options are 0, 1, ..., N-1 where N is the thetaBasisSize previously set'];
    pertThetaIndex = input([prompt,newline]);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    clear prompt
    prompt = [newline,'Input the desired r-index for the fixed k-mode', newline, '(Should be on the lower end of the range of values 1, 2, ..., M where M is the rBasisSize previously set'];
    pertRIndex = input([prompt,newline]);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pertThetaIndex = 3;
    pertRIndex = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pertTrigIndex = 0;


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


%If our fixed k-mode is theta-index zero eigenfunction 
% (i.e. phi_k = phi_k*(r)) then we only need to go through the whole
% computational process for a constant angle slice of the domain. This is
% easier to just handle with a separate main script for the radial case

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if pertThetaIndex == 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    run("TwoD_Disk_TimeEvo_Isentropic_Main_Radial")

    return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    thetaGrid = linspace(0,2*pi/(pertThetaIndex),thetaPoints);
    thetaGridFull = linspace(0,2*pi,(pertThetaIndex)*thetaPoints);
    thetaPointsFull = length(thetaGridFull);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NOTE MAY 4:
%   Will now use the fact that phi_k = cos(k*w)f(r) ->
%   -> p = sum[ cos(k*n*w)f_n(r) ] -> p (2pi/k)-periodic (angularly)
%   So we need only use the wedge of the circle from 0 to 2pi/k throughout
%   the problem.


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


pBar = 101325;
pConstIntMatrix = pBar*ones(thetaPoints,rPoints);

%Prompt user for a temperature with which to determine the baseline air
%   density to go along with the baseline air pressure pBar = 101325 Pa

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:radialEntropyLevels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    clear prompt
    prompt = [newline,'Specify the temperature (in degrees Fahrenheit, between 50 and 80 F)', newline, 'to determine baseline air density (and thus specific volume): '];
    backgroundTemperature = input([prompt,newline]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seaLevelAirCp = [1003.84,1003.87,1003.89,1003.92,1003.94,1003.97,1004.00,1004.02,1004.05,1004.08,1004.11,1004.13,1004.16,1004.19,1004.22,1004.25,1004.28,1004.30,1004.33,1004.36,1004.39,1004.42,1004.45,1004.48,1004.51,1004.55,1004.58,1004.61,1004.64,1004.67,1004.70];
seaLevelAirCv = [716.79,716.82,716.84,716.87,716.89,716.92,716.95,716.97,717.00,717.03,717.06,717.08,717.11,717.14,717.17,717.20,717.23,717.25,717.28,717.31,717.34,717.37,717.40,717.43,717.46,717.50,717.53,717.56,717.59,717.62,717.65];
seaLevelAirGamma = seaLevelAirCp./seaLevelAirCv;
seaLevelAirDensity = [1.246,1.244,1.241,1.239,1.237,1.234,1.232,1.229,1.227,1.225,1.222,1.22,1.218,1.215,1.213,1.211,1.208,1.206,1.204,1.202,1.199,1.197,1.195,1.193,1.19,1.188,1.186,1.184,1.182,1.179,1.177];

rhoBar = (1-mod(backgroundTemperature,1))*seaLevelAirDensity(floor(backgroundTemperature)-49) + ...
    mod(backgroundTemperature,1) * seaLevelAirDensity(floor(backgroundTemperature)-49 + 1);
Cp = (1-mod(backgroundTemperature,1))*seaLevelAirCp(floor(backgroundTemperature)-49) + ...
    mod(backgroundTemperature,1) * seaLevelAirCp(floor(backgroundTemperature)-49 + 1);
Cv = (1-mod(backgroundTemperature,1))*seaLevelAirCv(floor(backgroundTemperature)-49) + ...
    mod(backgroundTemperature,1) * seaLevelAirCv(floor(backgroundTemperature)-49 + 1);
gamma = Cp/Cv;
vBar = 1/rhoBar;

%Set sigmaBar = sigmaBar(pBar,s)  ( the inverse of the wave speed c )
%   by sigmaBar = sqrt(-dv/dp)|_{p=pBar} with equation of state
%   v(p,s) = vBar * [ (p/pBar)^(-1/gamma) ] * [e^(-s/Cp)]
%   where the last term can be ignored, as this is isentropic, s(x) = 0
sigmaBarConst = (1/sqrt(gamma))*sqrt(vBar)*(1/sqrt(pBar));


%Now generate the eigenfunction basis for p and v (which are the same
%   because the v-basis is sigmaBar^2*pBasis, but since this is the
%   isentropic case, sigmaBar^2=const does nothing after normalizing)
[pEigenBasis,pEigenFncZero,pBasisWeight,vEigenBasis,vEigenFncZero,vBasisWeight,eigenFreqMatrix,~] = ...
    TwoD_Disk_TimeEvo_Isentropic_Basis(rBasisSize,thetaBasisSize,sigmaBarConst,domainRadius,rGrid, ...
    rPoints,thetaGrid,thetaPoints,pertThetaIndex,pertTrigIndex);


%Again since we only use theta-indicies that are integer multiples of the
%   pertThetaIndex, i.e. n*[0:N], to separately save the eigenfunction to be
%   perturbed, we want the eigenfunction on the 2nd row and pertRIndex-th 
%   column of the basis matrix
pEigenFncPertTemp(:,:) = pEigenBasis(2,pertRIndex,:,:);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotPerturbationEigenFnc = false;

if plotPerturbationEigenFnc == true

[thetaMeshGridPlot,rMeshGridPlot] = meshgrid(thetaGrid,rGrid);
mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot), ...
    transpose(pEigenFncPertTemp));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


%Set percentage of pBar chosen k-mode eigenfunction should perturb the
%   background state by:

%-*-%-*-%-*-%-*-%-*
   alpha = 0.04;  %
%-*-%-*-%-*-%-*-%-*

perturbationCoeff = (alpha/abs((min(min(pEigenFncPertTemp)))))*pBar;

%Set the time period so that the eigenfunction we are perturbing the
%   background pressure by is a solution of the linearization around pBar.
%   We consider our solution as being 2T-periodic so set T = pi/lambda_k.
timePeriod = pi/eigenFreqMatrix(2,pertRIndex);
timeGrid = linspace(0,timePeriod,timePoints);
timeStepSize = timeGrid(2) - timeGrid(1);

%Set initial data for u, which will just be u(0,x) = 0 -> u_k(0) = 0
uIntDecompMatrixEvo = zeros(thetaBasisSize,rBasisSize);

%The unstored variables in the NL evolution are the matrix and decomp 
%   variables for v in case thsoe are needed for something later:
%   [pNLMatrix,pNLDecompMatrix,pNLZeroDecomp,vNLMatrix,vNLDecompMatrix,vNLZeroDecomp,uNLDecompMatrix]


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if performBasicIteration == true
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    if performInitialGuessCorrection == true
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

        pNLMatrix = zeros(thetaPoints,rPoints,timePoints,numberIterations+2);
        pNLDecompMatrix = zeros(thetaBasisSize,rBasisSize,timePoints,numberIterations+2);
        pNLZeroDecomp = zeros(timePoints,numberIterations+2);
        vNLMatrix = zeros(thetaPoints,rPoints,timePoints,numberIterations+2);
        vNLDecompMatrix = zeros(thetaBasisSize,rBasisSize,timePoints,numberIterations+2);
        vNLZeroDecomp = zeros(timePoints,numberIterations+2);
        uNLDecompMatrix = zeros(thetaBasisSize,rBasisSize,timePoints,numberIterations+2);
        pIntMatrixIteration = zeros(thetaPoints,rPoints,numberIterations+2);

        pIntMatrixIteration(:,:,1) = pConstIntMatrix + perturbationCoeff*pEigenFncPertTemp;

        [pIntMatrixIteration(:,:,2),pNLDecompMatrix(:,:,1,2),pNLZeroDecomp(1,2)] = ...
            TwoD_Disk_TimeEvo_Isentropic_InitialCorrection(pIntMatrixIteration(:,:,1),pEigenBasis, ...
            pEigenFncZero,pBasisWeight,rBasisSize,thetaBasisSize,rGrid,rPoints,thetaGrid,thetaPoints, ...
            pEigenFncPertTemp,eigenFreqMatrix,perturbationCoeff,pertThetaIndex,pertRIndex, ...
            pertTrigIndex,pBar,vBar,sigmaBarConst,gamma,timePeriod);

        tStart_NLEvolution = tic;
        [pNLMatrix(:,:,:,1),pNLDecompMatrix(:,:,:,1),pNLZeroDecomp(:,1),vNLMatrix(:,:,:,1),vNLDecompMatrix(:,:,:,1),vNLZeroDecomp(:,1),uNLDecompMatrix(:,:,:,1)] = ...
            TwoD_Disk_TimeEvo_Isentropic_NLEvolution(pIntMatrixIteration(:,:,1),uIntDecompMatrixEvo,rBasisSize,rGrid,rPoints, ...
            thetaBasisSize,thetaGrid,thetaPoints,timeGrid,timePoints,timeStepSize,pEigenBasis,pEigenFncZero,...
            pBasisWeight,vEigenBasis,vEigenFncZero,vBasisWeight,eigenFreqMatrix,sigmaBarConst,gamma,pertRIndex,pertThetaIndex,pBar,vBar);
        elapsed_NLEvolution = toc(tStart_NLEvolution)

        tStart_NLEvolution = tic;
        [pNLMatrix(:,:,:,2),pNLDecompMatrix(:,:,:,2),pNLZeroDecomp(:,2),vNLMatrix(:,:,:,2),vNLDecompMatrix(:,:,:,2),vNLZeroDecomp(:,2),uNLDecompMatrix(:,:,:,2)] = ...
            TwoD_Disk_TimeEvo_Isentropic_NLEvolution(pIntMatrixIteration(:,:,2),uIntDecompMatrixEvo,rBasisSize,rGrid,rPoints, ...
            thetaBasisSize,thetaGrid,thetaPoints,timeGrid,timePoints,timeStepSize,pEigenBasis,pEigenFncZero,...
            pBasisWeight,vEigenBasis,vEigenFncZero,vBasisWeight,eigenFreqMatrix,sigmaBarConst,gamma,pertRIndex,pertThetaIndex,pBar,vBar);
        elapsed_NLEvolution = toc(tStart_NLEvolution)

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    if performInitialGuessCorrection == false
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

        pNLMatrix = zeros(thetaPoints,rPoints,timePoints,numberIterations+1);
        pNLDecompMatrix = zeros(thetaBasisSize,rBasisSize,timePoints,numberIterations+1);
        pNLZeroDecomp = zeros(timePoints,numberIterations+1);
        vNLMatrix = zeros(thetaPoints,rPoints,timePoints,numberIterations+1);
        vNLDecompMatrix = zeros(thetaBasisSize,rBasisSize,timePoints,numberIterations+1);
        vNLZeroDecomp = zeros(timePoints,numberIterations+1);
        uNLDecompMatrix = zeros(thetaBasisSize,rBasisSize,timePoints,numberIterations+1);
        pIntMatrixIteration = zeros(thetaPoints,rPoints,numberIterations+1);

        pIntMatrixIteration(:,:,1) = pConstIntMatrix + perturbationCoeff*pEigenFncPertTemp;

        tStart_NLEvolution = tic;
        [pNLMatrix(:,:,:,1),pNLDecompMatrix(:,:,:,1),pNLZeroDecomp(:,1),vNLMatrix(:,:,:,1),vNLDecompMatrix(:,:,:,1),vNLZeroDecomp(:,1),uNLDecompMatrix(:,:,:,1)] = ...
            TwoD_Disk_TimeEvo_Isentropic_NLEvolution(pIntMatrixIteration(:,:,1),uIntDecompMatrixEvo,rBasisSize,rGrid,rPoints, ...
            thetaBasisSize,thetaGrid,thetaPoints,timeGrid,timePoints,timeStepSize,pEigenBasis,pEigenFncZero,...
            pBasisWeight,vEigenBasis,vEigenFncZero,vBasisWeight,eigenFreqMatrix,sigmaBarConst,gamma,pertRIndex,pertThetaIndex,pBar,vBar);
        elapsed_NLEvolution = toc(tStart_NLEvolution)

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if performBasicIteration == false
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pNLMatrix = zeros(thetaPoints,rPoints,timePoints);
    pNLDecompMatrix = zeros(thetaBasisSize,rBasisSize,timePoints);
    pNLZeroDecomp = zeros(timePoints,1);
    uNLDecompMatrix = zeros(thetaBasisSize,rBasisSize,timePoints);
    pIntMatrix = zeros(thetaPoints,rPoints);

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    if performInitialGuessCorrection == true
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

        pNLMatrixCorrected = zeros(thetaPoints,rPoints,timePoints);
        pNLDecompMatrixCorrected = zeros(thetaBasisSize,rBasisSize,timePoints);
        pNLZeroDecompCorrected = zeros(timePoints,1);
        uNLDecompMatrixCorrected = zeros(thetaBasisSize,rBasisSize,timePoints);
        pIntMatrix = pConstIntMatrix + perturbationCoeff*pEigenFncPertTemp;
        pIntMatrixCorrected = zeros(thetaPoints,rPoints);

        [pIntMatrixCorrected(:,:),pNLDecompMatrixCorrected(:,:,1),pNLZeroDecompCorrected(1,1)] = ...
            TwoD_Disk_TimeEvo_Isentropic_InitialCorrection(pIntMatrix,pEigenBasis, ...
            pEigenFncZero,rBasisSize,thetaBasisSize,rGrid,rPoints,thetaGrid, ...
            thetaPoints,pEigenFncPertTemp,eigenFreqMatrix,perturbationCoeff,pertThetaIndex, ...
            pertRIndex,pertTrigIndex,pBar,sigmaBarConst,gamma);

        tStart_NLEvolution = tic;
        [pNLMatrix,pNLDecompMatrix,pNLZeroDecomp,~,~,~,uNLDecompMatrix] = ...
            TwoD_Disk_TimeEvo_Isentropic_NLEvolution(pIntMatrix,uIntDecompMatrixEvo, ...
            rBasisSize,rGrid,rPoints,thetaBasisSize,thetaGrid,thetaPoints,timeGrid, ...
            timePoints,timeStepSize,pEigenBasis,pEigenFncZero,pBasisWeight,vEigenBasis, ...
            vEigenFncZero,vBasisWeight,eigenFreqMatrix,sigmaBarConst,gamma,pertRIndex,pertThetaIndex,pBar,vBar);
        elapsed_NLEvolution = toc(tStart_NLEvolution)

        tStart_NLEvolution = tic;
        [pNLMatrixCorrected,pNLDecompMatrixCorrected,pNLZeroDecompCorrected,~,~,~,uNLDecompMatrixCorrected] = ...
            TwoD_Disk_TimeEvo_Isentropic_NLEvolution(pIntMatrixCorrected, ...
            uIntDecompMatrixEvo,rBasisSize,rGrid,rPoints,thetaBasisSize,thetaGrid, ...
            thetaPoints,timeGrid,timePoints,timeStepSize,pEigenBasis,pEigenFncZero,...
            pBasisWeight,vEigenBasis,vEigenFncZero,vBasisWeight,eigenFreqMatrix, ...
            sigmaBarConst,gamma,pertRIndex,pertThetaIndex,pBar,vBar);
        elapsed_NLEvolution = toc(tStart_NLEvolution)

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    if performInitialGuessCorrection == false
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

        pIntMatrix = pConstIntMatrix + perturbationCoeff*pEigenFncPertTemp;

        tStart_NLEvolution = tic;
        [pNLMatrix,pNLDecompMatrix,pNLZeroDecomp,~,~,~,uNLDecompMatrix] = ...
            TwoD_Disk_TimeEvo_Isentropic_NLEvolution(pIntMatrix,uIntDecompMatrixEvo, ...
            rBasisSize,rGrid,rPoints,thetaBasisSize,thetaGrid,thetaPoints,timeGrid, ...
            timePoints,timeStepSize,pEigenBasis,pEigenFncZero,pBasisWeight,vEigenBasis, ...
            vEigenFncZero,vBasisWeight,eigenFreqMatrix,sigmaBarConst,gamma,pertRIndex,pertThetaIndex,pBar,vBar);
        elapsed_NLEvolution = toc(tStart_NLEvolution)

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 






%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
if performBasicIteration == true
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    %Our iteration by which we refine a solution to the Euler equations is:
    %   P^-1*f(y^(k)_0) = P^-1*Df(y^(k)_0)[r^(k)]
    %   where f(y^(k)_0) is a column vector representation of u^(k)(T,x),
    %   Df(y^(k)_0)[r^k] is a matrix mapping P^(k)(0,x) to u^(k)(T,x),
    %   and P is a suitable preconditioning matrix

    %NOTE 9/18/25:
    %   In the case of the 2D disc (at least for isentropic?) the only resonant
    %   eigenfunction to any given eigenfunction is that with the same theta-index
    %   and r-index but opposite trig function sin/cos. Would this be considered
    %   resonance (or some kind of it)? What is their iteraction term?
    %What we can say for now is this:
    %   this means we never need to consider how to separately handle resonant
    %   eigenfunctions as every eigenfunction we are using will be nonresonant
    %   to our chosen starting k-mode (at least in isentropic, research needed)

    smallDivisorsBasicNoZero = zeros(thetaBasisSize,rBasisSize);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for n = 1:thetaBasisSize 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        for m = 1:rBasisSize
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

            smallDivisorsBasicNoZero(n,m) = sin(eigenFreqMatrix(n,m)*timePeriod);

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        end
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %For every j-mode for j ~= k (our chosen fixed mode), we are able to solve
    %   the j-mode equation with the j-mode.
    %For the k-mode equation however, we do not have a k-mode variable to solve
    %   it with. Instead this becomes a bifurcation problem, and we use the
    %   0-mode to solve the k-mode equation.
    %   From Duhamel, we find that the k-mode effect due to adding
    %   1*phi_(00) to our initial guess pBar + c*phi_k is
    %   f(pBar + c*phi_k * phi_0) = f(pBar) + Df(pBar)[c*phi_k + phi_0]
    %   + c^2/2*D^2f(pBar)[phi_k,phi_k] + 1/2 D^2f(pBar)[phi_0,phi_0]
    %   + c D^2f(pBar)[phi_0,phi_k]

    %Only have alpha*D^2f(pBar)[phi_0,phi_k] which was found to be
    %   alpha*D^2f(pBar)[phi_0,phi_k] = -alpha*(pi/2)*|phi_0|*v''(pBar)
    %   = -(pi/2)*[(1+gamma)/(gamma*pBar)]*sigmaBar^2*alpha
    %   Here we are assuming isentropic so phi_0 is a scalar, so |phi_0| is just
    %   meant to mean the value of the zero eigenfunction at any point
    
    %NOTE APRIL 25: pBar vs pBar in calculation of these second derivative
    %   constants/small divisors. Why are we linearizing around pBar instead of
    %   the pBar we are actually using? Seems like we need to use pBar though
    %   not important until we consider non-isentropic since in isentropic pBar =
    %   pBar though
    %We store delta_(0,0) as just <D^2f(pBar)[phi_(00),phi_k],psi_k>
    %   as this is the same value as <D^2f(pBar)[phi_(00),phi_j],psi_j> for phi_j
    %   that are resonant with phi_k (lambda_k = lambda_j), so the alpha will be
    %   added on where needed

    %smallDivisorZeroMode = (-1) * ...
        %(((1+gamma)*vBar*eigenFreqMatrix(2,pertRIndex)*timePeriod)/(2*(gamma^2)*(pBar^2))) * ...
        %pEigenFncZero(1,1);

    smallDivisorZeroMode = (-1) * ...
        ( ( (1+gamma)*vBar*eigenFreqMatrix(2,pertRIndex)*timePeriod ) / ( 2*(gamma^2)*(pBar) ) ) * ...
        ( 1/sigmaBarConst );


    %NOTE APRIL 25: Change this to loop over arbitrary nubmer of iterations
    %   after testing single basic iteration. Otherwise need to set up function
    %   for linearization evolution around nonlienar state

    zeroCorrectionBasic = zeros(1,numberIterations);
    nonResCorrectionBasic = zeros(thetaBasisSize,rBasisSize,numberIterations);
    smallDivisorsPBarPlusZ = zeros(thetaBasisSize,rBasisSize,numberIterations);
    nonResCorrectionPBarPlusZ = zeros(thetaBasisSize,rBasisSize,numberIterations);
    pEigenFncIterationLoop(:,:) = zeros(thetaPoints,rPoints);



    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~performInitialGuessCorrection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    
    
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        for i = 1:numberIterations
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

            tStart_SolutionIteration = tic;

            pIntMatrixIteration(:,:,i+1) = pIntMatrixIteration(:,:,i);
    
            %We start by finding the corrections for the psi_k term, 
            %   which will be solved with the zero eigenfunction (where phi_k is our
            %   perturbed eigenfunction)
            %The phi_k equation for which we pick the z in +z*phi_(00) to solve is
            %   z*alpha * D^2f(pBar)[phi_(00),phi_k] = -<f(p(0,x)),psi_k>
    
            zeroCorrectionBasic(i) = (-1)*uNLDecompMatrix(2,pertRIndex,timePoints,i) /...
                (perturbationCoeff * smallDivisorZeroMode);
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for n = 1:thetaBasisSize
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
                for m = 1:rBasisSize
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
       
                    pEigenFncIterationLoop(:,:) = pEigenBasis(n,m,:,:);

                    %For the nonresonant modes, the basic iteration will just
                    %   attempt to correct the psi_j terms by choosing beta_j such that
                    %   beta_j*<Df(pBar)[phi_j],psi_j> = -<f(p(0,x)),psi_j>
                    %   where p(0,x) is the previous guess for the initial data that is
                    %   being improved in this iteration

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if n ~= 2 || m ~= pertRIndex
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        nonResCorrectionBasic(n,m,i)=(-1)*uNLDecompMatrix(n,m,timePoints,i)/ ...
                            smallDivisorsBasicNoZero(n,m);

                        pIntMatrixIteration(:,:,i+1) = pIntMatrixIteration(:,:,i+1) + ...
                            nonResCorrectionBasic(n,m,i) * pEigenFncIterationLoop(:,:);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    else
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    nonResCorrectionBasic(n,m,i)=0;

                    pIntMatrixIteration(:,:,i+1) = pIntMatrixIteration(:,:,i+1) + (zeroCorrectionBasic(i)*...
                        pEigenFncZero(:,:));

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
                end
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            [pNLMatrix(:,:,:,i+1),pNLDecompMatrix(:,:,:,i+1),pNLZeroDecomp(:,i+1),vNLMatrix(:,:,:,i+1),vNLDecompMatrix(:,:,:,i+1),vNLZeroDecomp(:,i+1),uNLDecompMatrix(:,:,:,i+1)] = ...
                TwoD_Disk_TimeEvo_Isentropic_NLEvolution(pIntMatrixIteration(:,:,i+1),uIntDecompMatrixEvo,rBasisSize,rGrid,rPoints, ...
                thetaBasisSize,thetaGrid,thetaPoints,timeGrid,timePoints,timeStepSize,pEigenBasis,pEigenFncZero,...
                pBasisWeight,vEigenBasis,vEigenFncZero,vBasisWeight,eigenFreqMatrix,sigmaBarConst,gamma,pertRIndex,pertThetaIndex,pBar,vBar);


            elapsed_SolutionIteration = toc(tStart_SolutionIteration)

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        end
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 




    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 





    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if performInitialGuessCorrection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        for i = 2:numberIterations+1
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

            tStart_SolutionIteration = tic;

            pIntMatrixIteration(:,:,i+1) = pIntMatrixIteration(:,:,i);

            %We start by finding the corrections for the psi_k term, 
            %   which will be solved with the zero eigenfunction (where phi_k is our
            %   perturbed eigenfunction)
            %The phi_k equation for which we pick the z in +z*phi_(00) to solve is
            %   z*alpha * D^2f(pBar)[phi_(00),phi_k] = -<f(p(0,x)),psi_k>

            zeroCorrectionBasic(i) = (-1)*uNLDecompMatrix(2,pertRIndex,timePoints,i) /...
                (perturbationCoeff * smallDivisorZeroMode);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for n = 1:thetaBasisSize
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
                for m = 1:rBasisSize
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
       
                     pEigenFncIterationLoop(:,:) = pEigenBasis(n,m,:,:);

                    %For the nonresonant modes, the basic iteration will just
                    %   attempt to correct the psi_j terms by choosing beta_j such that
                    %   beta_j*<Df(pBar)[phi_j],psi_j> = -<f(p(0,x)),psi_j>
                    %   where p(0,x) is the previous guess for the initial data that is
                    %   being improved in this iteration

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if n ~= 2 || m ~= pertRIndex
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        nonResCorrectionBasic(n,m,i)=(-1)*uNLDecompMatrix(n,m,timePoints,i)/ ...
                            smallDivisorsBasicNoZero(n,m);

                        pIntMatrixIteration(:,:,i+1) = pIntMatrixIteration(:,:,i+1) + ...
                            nonResCorrectionBasic(n,m,i) * pEigenFncIterationLoop(:,:);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    else
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        nonResCorrectionBasic(n,m,i)=0;

                        pIntMatrixIteration(:,:,i+1) = pIntMatrixIteration(:,:,i+1) + (zeroCorrectionBasic(i)*...
                            pEigenFncZero(:,:));

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
                end
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            [pNLMatrix(:,:,:,i+1),pNLDecompMatrix(:,:,:,i+1),pNLZeroDecomp(:,i+1),~,~,~,uNLDecompMatrix(:,:,:,i+1)] = ...
                TwoD_Disk_TimeEvo_Isentropic_NLEvolution(pIntMatrixIteration(:,:,i+1),uIntDecompMatrixEvo,rBasisSize,rGrid,rPoints, ...
                thetaBasisSize,thetaGrid,thetaPoints,timeGrid,timePoints,timeStepSize,pEigenBasis,pEigenFncZero,...
                pBasisWeight,vEigenBasis,vEigenFncZero,vBasisWeight,eigenFreqMatrix,sigmaBarConst,gamma,pertRIndex,pertThetaIndex,pBar,vBar);

            elapsed_SolutionIteration = toc(tStart_SolutionIteration)

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        end
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%for n = 1:thetaPoints

%    p231Ray(n,:) = pEigenBasis(3,3,1,n,:);
%    p241Ray(n,:) = pEigenBasis(3,4,1,n,:);
%    p131Ray(n,:) = pEigenBasis(2,3,1,n,:);
%    p232Ray(n,:) = pEigenBasis(3,3,2,n,:);

%    pFncDiffRRayInt(n) = trapz(rGrid,(p231Ray(n,:).*p241Ray(n,:)).*rGrid);
%    pFncDiffThetaRayInt(n) = trapz(rGrid,(p231Ray(n,:).*p131Ray(n,:)).*rGrid);
%    pFncDiffTrigRayInt(n) = trapz(rGrid,(p231Ray(n,:).*p232Ray(n,:)).*rGrid);
%    pFncSameRayInt(n) = trapz(rGrid,(p231Ray(n,:).*p231Ray(n,:)).*rGrid);

%end

%pFncDiffRDoubleInt = trapz(thetaGrid,pFncDiffRRayInt);
%pFncDiffRInnerProd = trapz(thetaGrid,pFncDiffRRayInt);
%pFncDiffTrigInnerProd = trapz(thetaGrid,pFncDiffTrigRayInt);
%pFncSameInnerProd = trapz(thetaGrid,pFncSameRayInt);


%for n = 1:thetaPoints
%J11CosCosVectorRayInt(n) = (1/pi)*cos(thetaGrid(n))*cos(thetaGrid(n))*trapz(rGrid,(J11VectorTemp.*J11VectorTemp).*rGrid);
%J11SinSinVectorRayInt(n) = (1/pi)*sin(thetaGrid(n))*sin(thetaGrid(n))*trapz(rGrid,(J11VectorTemp.*J11VectorTemp).*rGrid);
%J11CossSinVectorRayInt(n) = (1/pi)*sin(thetaGrid(n))*cos(thetaGrid(n))*trapz(rGrid,(J11VectorTemp.*J11VectorTemp).*rGrid);
%end


%J11CosCosInnerProd = trapz(thetaGrid,J11CosCosVectorRayInt)
%J11SinSinInnerProd = trapz(thetaGrid,J11SinSinVectorRayInt)
%J11CosSinInnerProd = trapz(thetaGrid,J11CosSinVectorRayInt)


%clear pNLMatrix pNLDecompMatrix vNLMatrix vNLDecompMatrix



%for n = 1:thetaBasisSize
%    for m = 1:rBasisSize
%        JScaledVectorLoop(:) = JScaledValues(n,m,:);

%        JScaledValuesNorm(n,m) = sqrt(trapz(rGrid,(JScaledVectorLoop.*JScaledVectorLoop).*rGrid));
%        JScaledValuesNormalized(n,m,:) = JScaledValues(n,m,:)/JScaledValuesNorm(n,m);
%    end
%end

%JInnerProdChoice(:) = JScaledValuesNormalized(2,3,:);
%JInnerProdCheck = trapz(rGrid,(JInnerProdChoice.*JInnerProdChoice).*rGrid)