set(0,'DefaultFigureWindowStyle','docked')

[thetaMeshGridPlot,rMeshGridPlot] = meshgrid(thetaGrid,rGrid);

figTemp3 = figure

iterationStep = 9;

v3 = VideoWriter('Sept_28_2025_Linear_Int8_DiffFromLinear_SideBySide_K14Alpha004.mp4', 'MPEG-4');
open(v3);

%pIntStepNoKModeTemp = eigenFncRecompNoKMode2DDisc(pOrigNLDecompMatrix(:,:,1),pOrigNLZeroDecomp(1,1),...
 %   pEigenBasis,pEigenFncZero,rBasisSize,thetaBasisSize,rPoints,thetaPoints,pertThetaIndex,pertRIndex);

linearSolutionSlice(:,:)=zeros(thetaPoints,rPoints);
nonlinearSolutionSlice(:,:)=zeros(thetaPoints,rPoints);
nonlinearDiffFromLinear(:,:)=zeros(thetaPoints,rPoints);

linearSolutionSlice(:,:) = pNLMatrixFull(1,:,:,1);
nonlinearSolutionSlice(:,:) = pNLMatrixFull(1,:,:,iterationStep);
nonlinearDiffFromLinear(:,:) = nonlinearSolutionSlice(:,:) - linearSolutionSlice(:,:);

fractionArray = linspace(0,1,timePoints);



for l1 = 1:11

    for l2 = 1:timePoints

    for j1 = 1:pertThetaIndex-1

        for j2 = 1:length(thetaGridFull)

            pNLMatrixFull(l2,(j1-1)*thetaPoints+j2,[1:rPoints],l1) = pNLMatrix(l2,j2,[1:rPoints],l1);

        end

    end

    end

end




s1 = subplot(1,3,1);
mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(linearSolutionSlice));

axis manual
axis([-1 1 -1 1 97000 106000])
%view(-37.5,27);
%view(-39.5324,11.0115);
view(-37.5,18);


title(['Nonlinear evolution of $p^{(0)}(t,x)$'],'Interpreter','latex' );



s2 = subplot(1,3,2);
mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(nonlinearSolutionSlice));

axis manual
axis([-1 1 -1 1 97000 106000])
%view(-37.5,27);
%view(-39.5324,11.0115);
view(-37.5,18);


title(['Nonlinear evolution of $p^{(8)}(t,x)$'],'Interpreter','latex' );
xlabel(['$\bar{p}=',num2str(pBar),', \alpha \bar{p} = ',num2str(perturbationCoeff),', T=\frac{\pi}{\lambda_{(1,4)}}=$',num2str(timePeriod)],'Interpreter','latex');



s3 = subplot(1,3,3);
mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(nonlinearDiffFromLinear));

axis manual
axis([-1 1 -1 1 -600 600])
%view(-37.5,27);
%view(-39.5324,11.0115);
view(-37.5,18);


title(['Difference of nonlinear solution $p^{(8)}(t,x)$ from the linear solution $p^{(0)}(t,x)$, where $p^{(0)}(0,x) = \bar{p} + (\alpha \bar{p}) \phi_{(1,4)}$'],'Interpreter','latex' );
text(0.5,-0.5,700,['$t=0 \cdot T$'],'Interpreter','latex','BackgroundColor','white','FontSize',12,'FontWeight','bold','Color','black')




frameTemp = getframe(figTemp3);
writeVideo(v3,frameTemp); 
writeVideo(v3,frameTemp);
writeVideo(v3,frameTemp);


for l = 2:timePoints

    %pIntStepNoKModeTemp = eigenFncRecompNoKMode2DDisc(pNLDecompMatrix(:,:,1,l),pNLZeroDecomp(1,l),...
     %   pEigenBasis,pEigenFncZero,rBasisSize,thetaBasisSize,rPoints,thetaPoints,pertThetaIndex,pertRIndex);


    linearSolutionSlice(:,:) = pNLMatrixFull(l,:,:,1);
    nonlinearSolutionSlice(:,:) = pNLMatrixFull(l,:,:,iterationStep);
    nonlinearDiffFromLinear(:,:) = nonlinearSolutionSlice(:,:) - linearSolutionSlice(:,:);
    
    
    subplot(1,3,1);
    mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(linearSolutionSlice));

    axis manual
    axis([-1 1 -1 1 97000 106000])
    %view(-37.5,27);
    %view(-39.5324,11.0115);
    view(-37.5,18);


    title(['Nonlinear evolution of $p^{(0)}(t,x)$'],'Interpreter','latex' );



    subplot(1,3,2);
    mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(nonlinearSolutionSlice));

    axis manual
    axis([-1 1 -1 1 97000 106000])
    %view(-37.5,27);
    %view(-39.5324,11.0115);
    view(-37.5,18);

    
    title(['Nonlinear evolution of $p^{(8)}(t,x)$'],'Interpreter','latex' );
    xlabel(['$\bar{p}=',num2str(pBar),', \alpha \bar{p} = ',num2str(perturbationCoeff),', T=\frac{\pi}{\lambda_{(1,4)}}=$',num2str(timePeriod)],'Interpreter','latex');


    subplot(1,3,3);
    mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(nonlinearDiffFromLinear));

    axis manual
    axis([-1 1 -1 1 -600 600])
    %view(-37.5,27);
    %view(-39.5324,11.0115);
    view(-37.5,18);


    title(['Difference of nonlinear solution $p^{(8)}(t,x)$ from the linear solution $p^{(0)}(t,x)$, where $p^{(0)}(0,x) = \bar{p} + (\alpha \bar{p}) \phi_{(1,4)}$'],'Interpreter','latex' );
    text(0.5,-0.5,700,['$t=$',num2str(fractionArray(l)),'$\cdot T$'],'Interpreter','latex','BackgroundColor','white','FontSize',12,'FontWeight','bold','Color','black')


    frameTemp = getframe(figTemp3);
    writeVideo(v3,frameTemp);
    writeVideo(v3,frameTemp);
    writeVideo(v3,frameTemp);
    %writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);
    %writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);
    %writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);

end

for l = 1:timePoints-2

    %pIntStepNoKModeTemp = eigenFncRecompNoKMode2DDisc(pNLDecompMatrix(:,:,1,l),pNLZeroDecomp(1,l),...
     %   pEigenBasis,pEigenFncZero,rBasisSize,thetaBasisSize,rPoints,thetaPoints,pertThetaIndex,pertRIndex);

    linearSolutionSlice(:,:) = pNLMatrixFull(timePoints-l,:,:,1);
    nonlinearSolutionSlice(:,:) = pNLMatrixFull(timePoints-l,:,:,iterationStep);
    nonlinearDiffFromLinear(:,:) = nonlinearSolutionSlice(:,:) - linearSolutionSlice(:,:);
    
    subplot(1,3,1);
    mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(linearSolutionSlice));

    axis manual
    axis([-1 1 -1 1 97000 106000])
    %view(-37.5,27);
    %view(-39.5324,11.0115);
    view(-37.5,18);


    title(['Nonlinear evolution of $p^{(0)}(t,x)$'],'Interpreter','latex' );



    subplot(1,3,2);
    mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(nonlinearSolutionSlice));

    axis manual
    axis([-1 1 -1 1 97000 106000])
    %view(-37.5,27);
    %view(-39.5324,11.0115);
    view(-37.5,18);

    
    title(['Nonlinear evolution of $p^{(8)}(t,x)$'],'Interpreter','latex' );
    xlabel(['$\bar{p}=',num2str(pBar),', \alpha \bar{p} = ',num2str(perturbationCoeff),', T=\frac{\pi}{\lambda_{(1,4)}}=$',num2str(timePeriod)],'Interpreter','latex');


    subplot(1,3,3);
    mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(nonlinearDiffFromLinear));

    axis manual
    axis([-1 1 -1 1 -600 600])
    %view(-37.5,27);
    %view(-39.5324,11.0115);
    view(-37.5,18);


    title(['Difference of nonlinear solution $p^{(8)}(t,x)$ from the linear solution $p^{(0)}(t,x)$, where $p^{(0)}(0,x) = \bar{p} + (\alpha \bar{p}) \phi_{(1,4)}$'],'Interpreter','latex' );
    text(0.5,-0.5,700,['$t=$',num2str(1.0+fractionArray(l)),'$\cdot T$'],'Interpreter','latex','BackgroundColor','white','FontSize',12,'FontWeight','bold','Color','black')




    frameTemp = getframe(figTemp3);
    writeVideo(v3,frameTemp);
    writeVideo(v3,frameTemp);
    writeVideo(v3,frameTemp);

    
    %writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);
    %writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);
    %writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);writeVideo(v4,frameTemp);

end

close(v3)

