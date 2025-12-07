%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%##\ .--. \##\ .-'. \##\ .--. \##\ .--. \##\ .--. \##\ .--. \##\ .`-. \##\.--%%%%
%%%%%:::::.\::::::::.\::::::::.\::::::::.\::::::::.\::::::::.\::::::::.\:::::::%%%%%
%%%%-' \##\ `--' \##\ `.-' \##\ `--' \##\ `--' \##\ `--' \##\ `-.' \##\ `--' \##%%%%
%%%%% \#\ /""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7 \#\ %%%%%
%%%%\  .--/""7__/""                                             /""7__/""7\ .--.%%%%
%%%%%:::::/""7__/"   2-D TIME EVOLUTION - TIME SERIES ANALYSIS   \"7__/""7:::::%%%%%
%%%%\  .--/""7__/""  RECONSTRUCT P/CREATE FIXED POINT PROFILES  /""7__/""7\ .--.%%%%
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


clear prompt
prompt = [newline,'Input the r-index for the position at which to save the pressure values at over time (must be between 1 and rPoints = ',int2str(rPoints),'):'];
listenerRIndex = input([prompt,newline]);

clear prompt
prompt = [newline,'Input the theta-index for the position at which to save the pressure values at over time (must be between 1 and thetaPoints = ',int2str(thetaPoints),'):'];
listenerThetaIndex = input([prompt,newline]);

clear prompt
prompt = [newline,'Input the number of iteration steps to save the pressure profiles up to (that is, the pressure profiles will be saved for step 0 to step j):'];
jStep = input([prompt,newline]);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 



%pressureProfileFullStepJ = zeros(thetaPoints,rPoints,timePoints);
%pressureProfileFullStep0 = zeros(thetaPoints,rPoints,timePoints);
eigenFncLoop = zeros(thetaPoints,rPoints);
pressureProfilesFixedFullPeriod = zeros(timePoints+timePoints-1,jStep+1);
pressureProfilesFixedFullPeriodNoBackground = zeros(timePoints+timePoints-1,jStep+1);

for j = 1:timePoints

    for i = 0:jStep

        pressureProfilesFixedFullPeriod(j,i+1) = pressureProfilesFixedFullPeriod(j,i+1) + pNLZeroDecomp(1,i+1) + pEigenFncZero(listenerThetaIndex,listenerRIndex);

    end

end


for n = 1:thetaBasisSize

    for m = 1:rBasisSize

        for j = 1:timePoints

            for i = 0:jStep

            pressureProfilesFixedFullPeriod(j,i+1) = pressureProfilesFixedFullPeriod(j,i+1) + pNLDecompMatrix(n,m,j,i+1) * pEigenBasis(n,m,listenerThetaIndex,listenerRIndex);
            pressureProfilesFixedFullPeriodNoBackground(j,i+1) = pressureProfilesFixedFullPeriodNoBackground(j,i+1) + pNLDecompMatrix(n,m,j,i+1) * pEigenBasis(n,m,listenerThetaIndex,listenerRIndex);
 
            end

        end
    
    end

end


timePointsFullPeriod = timePoints + (timePoints-1);
timeGridFullPeriod = zeros(1,timePoints + (timePoints-1));
timeGridFullPeriod(1:timePoints) = timeGrid(1:timePoints);
timeGridFullPeriod(timePoints+1:timePointsFullPeriod) = timeGrid(2:timePoints) + timePeriod;

for i = 0:jStep

    for j = 1:timePoints-1

        pressureProfilesFixedFullPeriodNoBackground(timePoints+j,i+1) = pressureProfilesFixedFullPeriodNoBackground(timePoints-j,i+1);
        pressureProfilesFixedFullPeriod(timePoints+j,i+1) = pressureProfilesFixedFullPeriod(timePoints-j,i+1);

    end

end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear prompt
prompt = [newline,'Create a longer version of these pressure profiles to generate audio files with? (input 0 or 1):'];
createProfileForAudioFiles = input([prompt,newline]);

if createProfileForAudioFiles

    clear prompt
    prompt = [newline,'Input the number of periods to extend the pressure profiles to (1024 periods creates a file ~1 second long for average choice of time period):'];
    numberPeriods = input([prompt,newline]);

    timePointsExtended = timePointsFullPeriod + ((timePointsFullPeriod-1) * (numberPeriods-1) );

    pressureProfilesFixedExtended = zeros(timePointsExtended,jStep+1);
    pressureProfilesFixedExtended(1:timePointsFullPeriod,1:jStep+1) = pressureProfilesFixedFullPeriodNoBackground(1:timePointsFullPeriod,1:jStep+1);

    timeGridExtended = linspace(0,numberPeriods*2*timePeriod,timePointsExtended);

    
    for i = 1:numberPeriods-1

        pressureProfilesFixedExtended(timePointsFullPeriod + 1 + ((i-1)*(timePointsFullPeriod-1)):timePointsFullPeriod + (i*(timePointsFullPeriod-1)),1:jStep+1) = ...
            pressureProfilesFixedFullPeriodNoBackground(2:timePointsFullPeriod,1:jStep+1);
        
    end


    numberResamplePoints = round(numberPeriods*(2*timePeriod)/(1/96000));
    resampleGrid = linspace(0,timePeriod*numberPeriods,numberResamplePoints);
    pressureProfileFixedExtendedResampled = zeros(numberResamplePoints,jStep+1);

    for l = 0:jStep

        pressureProfileFixedExtendedResampled(1:numberResamplePoints,l+1) = interp1(timeGridExtended,pressureProfilesFixedExtended(1:timePointsExtended,l+1),resampleGrid);

    end


end


%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %


















