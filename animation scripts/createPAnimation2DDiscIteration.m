set(0,'DefaultFigureWindowStyle','docked')
figTemp2 = figure;

iterationStep = 7;

[thetaGridAnimate,rGridAnimate] = meshgrid(thetaGrid,rGrid);
%For drawing a mesh on a figure, we want the valueMatrix that is to be
%   plotted in mesh(xGrid,yGrid,valueMatrix)
pIntMatrixAnimation(:,:) = pNLMatrix(1,:,:,iterationStep);
fractionArray = linspace(0,1,timePoints);
%Create x-y mesh from the theta-r mesh that we discretized the disc into
s = mesh(rGridAnimate.*cos(thetaGridAnimate),rGridAnimate.*sin(thetaGridAnimate),transpose(pIntMatrixAnimation));

%NOTE: have to manually change the subscripts for the \phi term in the
%   title of the figure to whatever pertThetaIndex and pertRIndex are.
%   (cannot currently use something like int2str inside a string)
title(['Nonlinear evolution of $p^{(7)}(t,x)$, where $p^{(1)}(0,x) = \bar{p} + \alpha \bar{p} \phi_{(3,4)}$'],'Interpreter','latex' );
subtitle(['$\bar{p}=',num2str(pBar),', \alpha=',num2str(perturbationCoeff/pBar),', T=\frac{\pi}{\lambda_{(3,4)}}=$',num2str(timePeriod)],'Interpreter','latex' );
text(0.5,-0.5,108000,"$t=0 \cdot T$",'Interpreter','latex','BackgroundColor','white','FontSize',16,'FontWeight','bold','Color','black')
%NOTE: to set the axes, which is especially needed for side by side
%   animations of figures, need to manually check max(matrixToPlot) and
%   min(matrixToPlot) and then set axes just past those values.

axis manual
axis([-1 1 -1 1 96000 108000])
view(-37.5,30);
%axis tight;

v4 = VideoWriter('Sept_19_2025_7th_Iterate_NLEvo.mp4', 'MPEG-4');
open(v4);

frame = getframe(figTemp2);
writeVideo(v4,frame);
writeVideo(v4,frame);
writeVideo(v4,frame);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 2:timePoints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %NOTE: If animation is too fast, must manually come into this script
    %   and uncomment the line below. This will just mean for each time step
    %   snapshot of the matrix we are plotting, we save that image twice.
    %   writeVideo(v4,frame);

    pNLMatrixSliceTemp(:,:) = pNLMatrix(k,:,:,iterationStep);
    %pNLMatrixSliceTemp(:,:) = pOrigNLMatrix(k,:,:);
    
    clf
    mesh(rGridAnimate.*cos(thetaGridAnimate),rGridAnimate.*sin(thetaGridAnimate),transpose(pNLMatrixSliceTemp));
    axis manual
    axis([-1 1 -1 1 96000 108000])
    view(-37.5,30);
    title(['Nonlinear evolution of $p^{(7)}(t,x)$, where $p^{(1)}(0,x) = \bar{p} + \alpha \bar{p} \phi_{(3,4)}$'],'Interpreter','latex' );
    subtitle(['$\bar{p}=',num2str(pBar),', \alpha=',num2str(perturbationCoeff/pBar),', T=\frac{\pi}{\lambda_{(3,4)}}=$',num2str(timePeriod)],'Interpreter','latex' );
    text(0.5,-0.5,108000,['$t=$',num2str(fractionArray(k)),'$\cdot T$'],'Interpreter','latex','BackgroundColor','white','FontSize',16,'FontWeight','bold','Color','black')

    frame = getframe(figTemp2);
    writeVideo(v4,frame);
    writeVideo(v4,frame);
    writeVideo(v4,frame);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    close(v4);



%view(3);
%grid on
%colormap jet;
%xlabel('\bf{x}','fontsize',18);
%ylabel('\bf{y}','fontsize',18);
%zlabel('\bf{u (x,y,t)}','fontsize',18);
%s = surf(x,y,zeros(size(x)),'Facealpha','1'); % temporary surface to get the handle


%for i=1:numel(T)
%    t=T(i);
%    b=cos(t);
%    A=log(b);
%    s.ZData = c1./((x./(A1+3*t).^1./3+A-(A2+y)./(A1+3*t).^1./3-A4)*(A1+3*t).^1./3);
%    drawnow;
%    
%    frame = getframe(fig);
%    writeVideo(v, frame)
%end
%close(v);

%view(axes1,[-32.1630230572161 18.5776076555021]);













