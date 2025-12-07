%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


set(0,'DefaultFigureWindowStyle','docked')

%[thetaMeshGridPlot,rMeshGridPlot] = meshgrid(thetaGrid,rGrid);
[thetaMeshGridPlot,rMeshGridPlot] = meshgrid(thetaGridFull,rGrid);

figTemp3 = figure

iterationStep = 7;

v3 = VideoWriter('Oct_30_2025_Linear_Int10_DiffFromLinear_SideBySide_K23Alpha004_4Plots.mp4', 'MPEG-4');
open(v3);


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


pNLMatrixFullInt = zeros(length(thetaGridFull),rPoints,timePoints);
pNLMatrixFullStep = zeros(length(thetaGridFull),rPoints,timePoints);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for l2 = 1:timePoints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    for j1 = 1:pertThetaIndex-1
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j2 = 1:length(thetaGrid)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            pNLMatrixFullInt((j1-1)*thetaPoints+j2,[1:rPoints],l2) = pNLMatrix(j2,[1:rPoints],l2,1);
            pNLMatrixFullStep((j1-1)*thetaPoints+j2,[1:rPoints],l2) = pNLMatrix(j2,[1:rPoints],l2,iterationStep);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    pEigenFncPertFull = zeros(thetaPointsFull,rPoints);

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    for j1 = 1:pertThetaIndex-1
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j2 = 1:length(thetaGrid)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            pEigenFncPertFull((j1-1)*thetaPoints+j2,1:rPoints) = pEigenFncPertTemp(j2,1:rPoints);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 




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



linearSolutionSlice(:,:)=zeros(thetaPointsFull,rPoints);
nonlinearSolutionSlice(:,:)=zeros(thetaPointsFull,rPoints);
nonlinearDiffFromLinear(:,:)=zeros(thetaPointsFull,rPoints);

linearSolutionSlice(:,:) = pNLMatrixFullInt(:,:,1);
nonlinearSolutionSlice(:,:) = pNLMatrixFullStep(:,:,1);
nonlinearDiffFromLinear(:,:) = nonlinearSolutionSlice(:,:) - linearSolutionSlice(:,:);
nonlinearDiffFromWaveEqnSol(:,:) = nonlinearSolutionSlice(:,:) - pBar - perturbationCoeff*pEigenFncPertFull(:,:);

fractionArray = linspace(0,1,timePoints);



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


s1 = subplot(1,4,1);
mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(linearSolutionSlice));
axis manual
axis([-1 1 -1 1 97000 106000])
view(-37.5,20);
title(['Nonlinear evolution of $p^{(0)}(x,t)$'],'Interpreter','latex','FontSize',9.5 );


s2 = subplot(1,4,2);
mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(nonlinearSolutionSlice));
axis manual
axis([-1 1 -1 1 97000 106000])
view(-37.5,20);
title(['Nonlinear evolution of $p^{(8)}(x,t)$'],'Interpreter','latex','FontSize',9.5 );
%xlabel(['$\bar{p}=',num2str(pBar),', \alpha \bar{p} = ',num2str(perturbationCoeff),', T=\frac{\pi}{\lambda_{2,3}}=$',num2str(timePeriod)],'Interpreter','latex','FontSize',10);
xlabel(['$\bar{p}=',num2str(pBar),',\: \alpha \bar{p} = ',num2str(perturbationCoeff),',\: T=\frac{\pi}{\lambda_{2,3}}=$',num2str(timePeriod)],'Interpreter','latex','FontSize',10);


s3 = subplot(1,4,3);
mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(nonlinearDiffFromLinear));
axis manual
axis([-1 1 -1 1 -280 280])
view(-37.5,20);
title(['Difference of $p^{(8)}(x,t)$ and solution starting from $p^{(0)}(0,x) = \bar{p} + (\alpha \bar{p}) \phi_{2,3}$'],'Interpreter','latex','FontSize',9.5 );
text(0.5,-0.5,700,['$t=0 \cdot T$'],'Interpreter','latex','BackgroundColor','white','FontSize',12,'FontWeight','bold','Color','black')


s3 = subplot(1,4,4);
mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(nonlinearDiffFromWaveEqnSol));
axis manual
axis([-1 1 -1 1 -8200 8200])
view(-37.5,20);
title(['Difference of $p^{(8)}(x,t)$ and the linear wave equation solution $P^{(0)}(x,t) = \bar{p} + (\alpha \bar{p}) \cos(\lambda_{2,3}t) \phi_{2,3}(x)$'],'Interpreter','latex','FontSize',9.5 );
text(1.2,-1.2,0,['$t=0 \cdot T$'],'Interpreter','latex','BackgroundColor','white','FontSize',12,'FontWeight','bold','Color','black')


frameTemp = getframe(figTemp3);
writeVideo(v3,frameTemp); 
writeVideo(v3,frameTemp);
writeVideo(v3,frameTemp);


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for l = 2:timePoints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %pIntStepNoKModeTemp = eigenFncRecompNoKMode2DDisc(pNLDecompMatrix(:,:,1,l),pNLZeroDecomp(1,l),...
     %   pEigenBasis,pEigenFncZero,rBasisSize,thetaBasisSize,rPoints,thetaPoints,pertThetaIndex,pertRIndex);


    linearSolutionSlice(:,:) = pNLMatrixFullInt(:,:,l);
    nonlinearSolutionSlice(:,:) = pNLMatrixFullStep(:,:,l);
    nonlinearDiffFromLinear(:,:) = nonlinearSolutionSlice(:,:) - linearSolutionSlice(:,:);
    %nonlinearDiffFromWaveEqnSol(:,:) = nonlinearSolutionSlice(:,:) - (pBar - perturbationCoeff*pEigenFncPertFull(:,:)*cos(eigenFreqMatrix(2,pertRIndex)*(pertThetaIndex-1)*timeGrid(l)));
    nonlinearDiffFromWaveEqnSol(:,:) = nonlinearSolutionSlice(:,:) - (pBar - perturbationCoeff*pEigenFncPertFull(:,:)*cos(eigenFreqMatrix(2,pertRIndex)*timeGrid(l)));
    
    
    s1 = subplot(1,4,1);
mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(linearSolutionSlice));
axis manual
axis([-1 1 -1 1 97000 106000])
view(-37.5,20);
title(['Nonlinear evolution of $p^{(0)}(x,t)$'],'Interpreter','latex','FontSize',9.5 );


s2 = subplot(1,4,2);
mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(nonlinearSolutionSlice));
axis manual
axis([-1 1 -1 1 97000 106000])
view(-37.5,20);
title(['Nonlinear evolution of $p^{(8)}(x,t)$'],'Interpreter','latex','FontSize',9.5 );
%xlabel(['$\bar{p}=',num2str(pBar),', \alpha \bar{p} = ',num2str(perturbationCoeff),', T=\frac{\pi}{\lambda_{2,3}}=$',num2str(timePeriod)],'Interpreter','latex','FontSize',10);
xlabel(['$\bar{p}=',num2str(pBar),',\: \alpha \bar{p} = ',num2str(perturbationCoeff),',\: T=\frac{\pi}{\lambda_{2,3}}=$',num2str(timePeriod)],'Interpreter','latex','FontSize',10);


s3 = subplot(1,4,3);
mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(nonlinearDiffFromLinear));
axis manual
axis([-1 1 -1 1 -280 280])
view(-37.5,20);
title(['Difference of $p^{(8)}(x,t)$ and solution starting from $p^{(0)}(0,x) = \bar{p} + (\alpha \bar{p}) \phi_{2,3}$'],'Interpreter','latex','FontSize',9.5 );
text(0.5,-0.5,700,['$t=0 \cdot T$'],'Interpreter','latex','BackgroundColor','white','FontSize',12,'FontWeight','bold','Color','black')


s3 = subplot(1,4,4);
mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(nonlinearDiffFromWaveEqnSol));
axis manual
axis([-1 1 -1 1 -8200 8200])
view(-37.5,20);
title(['Difference of $p^{(8)}(x,t)$ and the linear wave equation solution $P^{(0)}(x,t) = \bar{p} + (\alpha \bar{p}) \cos(\lambda_{2,3}t) \phi_{2,3}(x)$'],'Interpreter','latex','FontSize',9.5 );
text(1.2,-1.2,0,['$t=$',num2str(fractionArray(l)),'$\cdot T$'],'Interpreter','latex','BackgroundColor','white','FontSize',12,'FontWeight','bold','Color','black')


    frameTemp = getframe(figTemp3);
    writeVideo(v3,frameTemp);
    writeVideo(v3,frameTemp);
    writeVideo(v3,frameTemp);
    %writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);
    %writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);
    %writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for l = 1:timePoints-2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %pIntStepNoKModeTemp = eigenFncRecompNoKMode2DDisc(pNLDecompMatrix(:,:,1,l),pNLZeroDecomp(1,l),...
     %   pEigenBasis,pEigenFncZero,rBasisSize,thetaBasisSize,rPoints,thetaPoints,pertThetaIndex,pertRIndex);

    linearSolutionSlice(:,:) = pNLMatrixFullInt(:,:,timePoints-l);
    nonlinearSolutionSlice(:,:) = pNLMatrixFullStep(:,:,timePoints-l);
    nonlinearDiffFromLinear(:,:) = nonlinearSolutionSlice(:,:) - linearSolutionSlice(:,:);
    %nonlinearDiffFromWaveEqnSol(:,:) = nonlinearSolutionSlice(:,:) - (pBar - perturbationCoeff*pEigenFncPertFull(:,:)*cos(eigenFreqMatrix(2,pertRIndex)*(pertThetaIndex-1)*timeGrid(timePoints-l)));\
    nonlinearDiffFromWaveEqnSol(:,:) = nonlinearSolutionSlice(:,:) - (pBar - perturbationCoeff*pEigenFncPertFull(:,:)*cos(eigenFreqMatrix(2,pertRIndex)*timeGrid(timePoints-l)));
    
    s1 = subplot(1,4,1);
mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(linearSolutionSlice));
axis manual
axis([-1 1 -1 1 97000 106000])
view(-37.5,20);
title(['Nonlinear evolution of $p^{(0)}(x,t)$'],'Interpreter','latex','FontSize',9.5 );


s2 = subplot(1,4,2);
mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(nonlinearSolutionSlice));
axis manual
axis([-1 1 -1 1 97000 106000])
view(-37.5,20);
title(['Nonlinear evolution of $p^{(8)}(x,t)$'],'Interpreter','latex','FontSize',9.5 );
xlabel(['$\bar{p}=',num2str(pBar),',\: \alpha \bar{p} = ',num2str(perturbationCoeff),',\: T=\frac{\pi}{\lambda_{2,3}}=$',num2str(timePeriod)],'Interpreter','latex','FontSize',10);


s3 = subplot(1,4,3);
mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(nonlinearDiffFromLinear));
axis manual
axis([-1 1 -1 1 -280 280])
view(-37.5,20);
title(['Difference of $p^{(8)}(x,t)$ and solution starting from $p^{(0)}(0,x) = \bar{p} + (\alpha \bar{p}) \phi_{2,3}$'],'Interpreter','latex','FontSize',9.5 );


s3 = subplot(1,4,4);
mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(nonlinearDiffFromWaveEqnSol));
axis manual
axis([-1 1 -1 1 -8200 8200])
view(-37.5,20);
title(['Difference of $p^{(8)}(x,t)$ and the linear wave equation solution $P^{(0)}(x,t) = \bar{p} + (\alpha \bar{p}) \cos(\lambda_{2,3}t) \phi_{2,3}(x)$'],'Interpreter','latex','FontSize',9.5 );
text(1.2,-1.2,0,['$t=$',num2str(1.0+fractionArray(l)),'$\cdot T$'],'Interpreter','latex','BackgroundColor','white','FontSize',12,'FontWeight','bold','Color','black')




    frameTemp = getframe(figTemp3);
    writeVideo(v3,frameTemp);
    writeVideo(v3,frameTemp);
    writeVideo(v3,frameTemp);

    
    %writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);
    %writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);
    %writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close(v3)





%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %









