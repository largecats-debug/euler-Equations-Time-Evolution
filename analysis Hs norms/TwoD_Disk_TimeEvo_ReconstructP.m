%specify which iteration step to look at the pressure profile of
jStep = 8;

%specify the location on the domain that we want to look at
pThetaIndex = 220; %from 1 to thetaPoints
pRIndex = 379; %from 1 to rPoints

pressureProfileStep0 = zeros(timePoints,1);
pressureProfileStepJ = zeros(timePoints,1);
%pressureProfileFullStepJ = zeros(thetaPoints,rPoints,timePoints);
%pressureProfileFullStep0 = zeros(thetaPoints,rPoints,timePoints);
eigenFncLoop = zeros(thetaPoints,rPoints);

for n = 1:thetaBasisSize

    for m = 1:rBasisSize

        eigenFncLoop(:,:) = pEigenBasis(n,m,:,:);

        for j = 1:timePoints

            pressureProfileStepJ(j) = pressureProfileStepJ(j) + pNLDecompMatrix(n,m,j,jStep) * pEigenBasis(n,m,pThetaIndex,pRIndex);
            pressureProfileStep0(j) = pressureProfileStep0(j) + pNLDecompMatrix(n,m,j,1) * pEigenBasis(n,m,pThetaIndex,pRIndex);

            %pressureProfileOmegaStepJ(:,:,j) = pressureProfileOmegaStepJ(:,:,j) + pNLDecompMatrix(n,m,j,jStep) * eigenFncLoop(:,:);
            %pressureProfileOmegaStep0(:,:,j) = pressureProfileOmegaStep0(:,:,j) + pNLDecompMatrix(n,m,j,1) * eigenFncLoop(:,:);

        end
    
    end

end