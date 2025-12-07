


timePointsFullPeriod = timePoints + (timePoints-1);
timeGridFullPeriod = zeros(1,timePoints + (timePoints-1));
timeGridFullPeriod(1:timePoints) = timeGrid(1:timePoints);
timeGridFullPeriod(timePoints+1:timePointsFullPeriod) = timeGrid(2:timePoints) + timePeriod;

pressureProfileFullStepJ = zeros(timePointsFullPeriod,1);
pressureProfileFullStepJ(1:timePoints) = pressureProfileStepJ(1:timePoints);
pressureProfileFullStep0 = zeros(timePointsFullPeriod,1);
pressureProfileFullStep0(1:timePoints) = pressureProfileStep0(1:timePoints);

for j = 1:timePoints-1

    pressureProfileFullStepJ(timePoints+j) = pressureProfileStepJ(timePoints-j);
    pressureProfileFullStep0(timePoints+j) = pressureProfileStep0(timePoints-j);

end



numberPeriods = 1024;

pressureProfileExtendedStep0 = zeros(timePointsFullPeriod + ((timePointsFullPeriod-1) * (numberPeriods-1) ),1);
pressureProfileExtendedStep0(1:timePointsFullPeriod) = pressureProfileFullStep0(1:timePointsFullPeriod);
pressureProfileExtendedStepJ = zeros(timePointsFullPeriod + ((timePointsFullPeriod-1) * (numberPeriods-1) ),1);
pressureProfileExtendedStepJ(1:timePointsFullPeriod) = pressureProfileFullStepJ(1:timePointsFullPeriod);

%timeGridExtended = zeros(1,timePointsFullPeriod + ((timePointsFullPeriod-1) * (numberPeriods-1) ));
%timeGridExtended(1:timePointsFullPeriod) = timeGridFullPeriod(1:timePointsFullPeriod);
timeGridExtended = linspace(0,numberPeriods*2*timePeriod,length(pressureProfileExtendedStep0));
for i = 1:numberPeriods-1

    pressureProfileExtendedStep0(timePointsFullPeriod + 1 + ((i-1)*(timePointsFullPeriod-1)):timePointsFullPeriod + (i*(timePointsFullPeriod-1))) = ...
        pressureProfileFullStep0(2:timePointsFullPeriod);

    pressureProfileExtendedStepJ(timePointsFullPeriod + 1 + ((i-1)*(timePointsFullPeriod-1)):timePointsFullPeriod + (i*(timePointsFullPeriod-1))) = ...
        pressureProfileFullStepJ(2:timePointsFullPeriod);

    %timeGridExtended(1 + (i*timePointsFullPeriod):timePointsFullPeriod + (i*timePointsFullPeriod)) = ...
        %timeGridFullPeriod(1:timePointsFullPeriod) + (i*timePeriod);

end


numberResamplePoints = round(numberPeriods*(2*timePeriod)/(1/96000));
resampleGrid = linspace(0,timePeriod*numberPeriods,numberResamplePoints);

pressureResampleStep0 = interp1(timeGridExtended,pressureProfileExtendedStep0,resampleGrid);
pressureResampleStepJ = interp1(timeGridExtended,pressureProfileExtendedStepJ,resampleGrid);




%audiowrite('Nov_07_2025_K32Alpha004_Step0Audio512.wav',pressureResampleStep0/max(abs(pressureResampleStep0)),96000)
%audiowrite('Nov_07_2025_K32Alpha004_Step8Audio512.wav',pressureResampleStepJ/max(abs(pressureResampleStepJ)),96000)
%audiowrite('Nov_07_2025_K32Alpha004_Step8Audio1024.wav',pressureResampleStepJ/max(abs(pressureResampleStepJ)),96000)
%audiowrite('Nov_07_2025_K32Alpha004_Step0Audio1024.wav',pressureResampleStep0/max(abs(pressureResampleStep0)),96000)









