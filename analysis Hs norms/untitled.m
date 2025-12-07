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
prompt = [newline,'Input the iteration step to create a fully reconstructed pressure profile for (that is, the pressure profiles will be saved for step 0 to step j):'];
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

        pressureProfilesFixedFullPeriod(j,i+1) = pressureProfilesFixedFullPeriod(j,i+1) + pNLZeroDecomp(1,i+1);

    end

end


for n = 1:thetaBasisSize

    for m = 1:rBasisSize

        eigenFncLoop(:,:) = pEigenBasis(n,m,:,:);

        for j = 1:timePoints

            for i = 0:jStep

            pressureProfilesFixedFullPeriod(j,i+1) = pressureProfilesFixedFullPeriod(j,i+1) + pNLDecompMatrix(n,m,j,i+1) * pEigenBasis(n,m,listenerThetaIndex,listenerRIndex);
            pressureProfilesFixedFullPeriodNoBackground = pressureProfilesFixedFullPeriodNoBackground(j,i+1) + pNLDecompMatrix(n,m,j,i+1) * pEigenBasis(n,m,listenerThetaIndex,listenerRIndex);
 
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



%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %


















