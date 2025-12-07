
[thetaMeshGridPlot,rMeshGridPlot] = meshgrid(thetaGrid,rGrid);
figTemp = figure;
pIntMatrixTemp(:,:) = zeros(thetaPoints,rPoints);

v4 = VideoWriter('Sept_28_2025_IntDataIteration_K14Alpha004_11Steps.mp4', 'MPEG-4');
open(v4);

pIntMatrixTemp(:,:) = pNLMatrix(1,:,:,1);
s = mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(pIntMatrixTemp));

title(['Initial $p^{(i)}(0,x)=\bar{p} + \alpha \bar{p} \phi_{(1,3)}$'],'Interpreter','latex');
    subtitle(['$\bar{p}=',num2str(pBar),', ',' \alpha=',num2str(perturbationCoeff/pBar),',  T=\frac{\pi}{\lambda_{(4,1)}}=$',num2str(timePeriod)],'Interpreter','latex' );
text(0.5,-0.5,106000,"$p^{(0)}(0,x)$",'Interpreter','latex','BackgroundColor','white','FontSize',16,'FontWeight','bold','Color','black')
axis manual
%axis([-1 1 -1 1 97000,106000]);
axis([-1 1 -1 1 96000,107000]);
view(-37.5,17);

frame = getframe(figTemp);
    writeVideo(v4,frame);
    writeVideo(v4,frame);
    writeVideo(v4,frame);
    writeVideo(v4,frame);
    writeVideo(v4,frame);
    writeVideo(v4,frame);

for i = 2:11

    %NOTE: have to manually change the subscripts for the \phi term in the
    %   title of the figure to whatever pertThetaIndex and pertRIndex are.
    %   (cannot currently use something like int2str inside a string)

    pIntMatrixTemp(:,:) = pNLMatrix(1,:,:,i);
    clf
    s = mesh(rMeshGridPlot.*cos(thetaMeshGridPlot),rMeshGridPlot.*sin(thetaMeshGridPlot),transpose(pIntMatrixTemp));

    %title([int2str(i),'-th initial data,starting from $p^{(0)}(0,x) = \bar{p} + \alpha \bar{p} \phi_{(4,1)}$'],'Interpreter','latex' );
    title(['Initial data $p^{(i)}(0,x)$ iteration, starting from $p^{(0)}(0,x)=\bar{p} + \alpha \bar{p} \phi_{(1,3)}$'],'Interpreter','latex');
    subtitle(['$\bar{p}=',num2str(pBar),', ',' \alpha=',num2str(perturbationCoeff/pBar),',  T=\frac{\pi}{\lambda_{(4,1)}}=$',num2str(timePeriod)],'Interpreter','latex' );
    %axis([-1 1 -1 1 97000,106000]);
    axis([-1 1 -1 1 96000,107000]);
    view(-37.5,17);

    %if i == 1
    %    text(0.5,-0.5,108000,"$p^{(1)}(0,x)$",'Interpreter','latex','BackgroundColor','white','FontSize',16,'FontWeight','bold','Color','black')
    %end

    if i == 2
        text(0.5,-0.5,106000,"$p^{(1)}(0,x)$",'Interpreter','latex','BackgroundColor','white','FontSize',16,'FontWeight','bold','Color','black')
    end

    if i == 3
        text(0.5,-0.5,106000,"$p^{(2)}(0,x)$",'Interpreter','latex','BackgroundColor','white','FontSize',16,'FontWeight','bold','Color','black')
    end

    if i == 4
        text(0.5,-0.5,106000,"$p^{(3)}(0,x)$",'Interpreter','latex','BackgroundColor','white','FontSize',16,'FontWeight','bold','Color','black')
    end

    if i == 5
        text(0.5,-0.5,106000,"$p^{(4)}(0,x)$",'Interpreter','latex','BackgroundColor','white','FontSize',16,'FontWeight','bold','Color','black')
    end

    if i == 6
        text(0.5,-0.5,106000,"$p^{(5)}(0,x)$",'Interpreter','latex','BackgroundColor','white','FontSize',16,'FontWeight','bold','Color','black')
    end

    if i == 7
        text(0.5,-0.5,106000,"$p^{(6)}(0,x)$",'Interpreter','latex','BackgroundColor','white','FontSize',16,'FontWeight','bold','Color','black')
    end

    if i == 8
        text(0.5,-0.5,106000,"$p^{(7)}(0,x)$",'Interpreter','latex','BackgroundColor','white','FontSize',16,'FontWeight','bold','Color','black')
    end

    if i == 9
        text(0.5,-0.5,106000,"$p^{(8)}(0,x)$",'Interpreter','latex','BackgroundColor','white','FontSize',16,'FontWeight','bold','Color','black')
    end

    if i == 10
        text(0.5,-0.5,106000,"$p^{(9)}(0,x)$",'Interpreter','latex','BackgroundColor','white','FontSize',16,'FontWeight','bold','Color','black')
    end

    if i == 11
        text(0.5,-0.5,106000,"$p^{(10)}(0,x)$",'Interpreter','latex','BackgroundColor','white','FontSize',16,'FontWeight','bold','Color','black')
    end

    frame = getframe(figTemp);
    writeVideo(v4,frame);
    writeVideo(v4,frame);
    writeVideo(v4,frame);
    writeVideo(v4,frame);
    writeVideo(v4,frame);
    writeVideo(v4,frame);
    writeVideo(v4,frame);
    writeVideo(v4,frame);
    writeVideo(v4,frame);
    writeVideo(v4,frame);

    

end

close(v4);
