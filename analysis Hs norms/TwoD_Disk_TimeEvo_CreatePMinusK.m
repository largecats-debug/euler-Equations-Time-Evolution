%specify which iteration step to look at the pressure profile of
pStep = 8;

pStepMinusK = zeros(thetaPoints,rPoints,timePoints);
eigenFncLoop = zeros(thetaPoints,rPoints);

%pStepMinusK = pStepMinusK + pNLZeroDecomp(1,pStep)*pEigenFncZero;

for n = 1:thetaBasisSize

    for m = 1:rBasisSize

        if n == 2 && m == pertRIndex
            continue
        end

        eigenFncLoop(:,:) = pEigenBasis(n,m,:,:);

        for j = 1:timePoints

            pStepMinusK(:,:,j) = pStepMinusK(:,:,j) + pNLDecompMatrix(n,m,j,pStep) * eigenFncLoop(:,:);

        end
    
    end

end


pStepMinusKFull = zeros(thetaPointsFull,rPoints,timePoints);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for l2 = 1:timePoints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    for j1 = 0:pertThetaIndex-1
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j2 = 1:length(thetaGrid)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            pStepMinusKFull((j1)*thetaPoints+j2,[1:rPoints],l2) = pStepMinusK(j2,[1:rPoints],l2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%pStepMinus0KFullStep8 = pStepMinusKFull - pBar;









