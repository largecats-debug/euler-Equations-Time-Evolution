[thetaMeshGridPlot,rMeshGridPlot] = meshgrid(thetaGridFull,rGrid);


pIntMatrixTemp2(:,:) = zeros((pertThetaIndex-1)*thetaPoints,rPoints);

for i = 1:11

    for j1 = 1:pertThetaIndex-1

        for j2 = 1:length(thetaGridFull)

            pNLMatrixFull(1,(j1-1)*thetaPoints+j2,[1:rPoints],i) = pNLMatrix(1,j2,[1:rPoints],i);

        end

    end

    pIntMatrixTemp2(:,:) = pNLMatrixFull(1,:,:,i);

    figure

    mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(pIntMatrixTemp2))

    %NOTE: have to manually change the subscripts for the \phi term in the
    %   title of the figure to whatever pertThetaIndex and pertRIndex are.
    %   (cannot currently use something like int2str inside a string)

    if i == 1
        title(['Initial data $p^{(0)}(t,x) = \bar{p} + \alpha \bar{p} \phi_{(2,3)}$'],'Interpreter','latex' );
        subtitle(['$\bar{p}=',num2str(pBar),', \alpha=',num2str(perturbationCoeff/pBar),', T=\frac{\pi}{\lambda_{(2,3)}}=$',num2str(timePeriod)],'Interpreter','latex' );
    end

    if i ~= 1
        title([int2str(i),'-th iterated initial data,starting from $p^{(0)}(0,x) = \bar{p} + \alpha \bar{p} \phi_{(2,3)}$'],'Interpreter','latex' );
        subtitle(['$\bar{p}=',num2str(pBar),', \alpha=',num2str(perturbationCoeff/pBar),', T=\frac{\pi}{\lambda_{(2,3)}}=$',num2str(timePeriod)],'Interpreter','latex' );
    end

end