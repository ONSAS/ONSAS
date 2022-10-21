%md# Static Von-Mises Truss example

close all; clear all; addpath( genpath( [ pwd '/../../src'] ) );
% scalar parameters
E = 210e9 ; L = sqrt(2) ;
A = 7.854e-05 ; ang1 = 65 ;
Kplas = 21e8;
sigma_Y_0 = 250e6 ;

global Rstress Rstrain Rstrainacum ;

x = [0 1 2];
y = [0 1 0];
z = [0 0 0];
s = [1 3];
t = [2 2];
weights = [2 2];
G = graph(s,t,weights);
figure(1);
p=plot(G,'XData',x,'YData',y,'linewidth',2);
p.NodeColor = 'r';
p.MarkerSize = 5;
title('Plasticidad / Reticulado de von Mises');
hold on;

%{
x = [0 1 2] ;
y = [0 1 0] ;
z = [0 0 0] ;

figure(1);
plot(x,y,'o-','markersize',4,'markerfacecolor','b','linewidth',2);
title('Plasticidad / Reticulado de von Mises');
hold on;
%}

% x and z coordinates of node 2
x2 = 1 ;
z2 = 1 ;

materials.hyperElasModel  = 'isotropicHardening' ;
materials.hyperElasParams = [ E Kplas sigma_Y_0 ] ;

elements(1).elemType = 'node' ;
elements(2).elemType = 'truss';
elements(3).elemType = 'truss';
elements(2).elemCrossSecParams = { 'circle' , 0.010} ;
elements(3).elemCrossSecParams = { 'circle' , 0.012} ;

% sqrt(A*4/pi)

boundaryConds(1).imposDispDofs = [ 1 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 ] ;

boundaryConds(2).imposDispDofs =  3 ;
boundaryConds(2).imposDispVals =  0 ;
boundaryConds(2).loadsCoordSys = 'global'         ;
boundaryConds(2).loadsTimeFact = @(t) 2.8e4*t     ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 1 0 ]  ;

initialConds                = struct() ;

mesh.nodesCoords = [   0  0   0 ; ...
                      x2  0  z2 ; ...
                    2*x2  0   0 ] ;

mesh.conecCell = cell(5,1) ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 1 0  3   ] ;
mesh.conecCell{ 3, 1 } = [ 0 1 2 0  2   ] ;
mesh.conecCell{ 4, 1 } = [ 1 2 0 0  1 2 ] ;
mesh.conecCell{ 5, 1 } = [ 1 3 0 0  2 3 ] ;

analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.deltaT        =   2e-2 ;
analysisSettings.finalTime     =   1 ;

analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   15   ;

analysisSettings.posVariableLoadBC = 2 ;

otherParams.problemName = 'static_plastic_von_mises_truss';
otherParams.plots_format = 'vtk' ;
otherParams.plots_deltaTs_separation = 2 ;

[matUs, loadFactorsMat ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

deltas = [-matUs(6+1,:)' -matUs(6+5,:)'] ;

    for t = 1:(analysisSettings.finalTime/analysisSettings.deltaT)
        x1(t)=Rstrain(t,1);
        y1(t)=Rstress(t,1);
    end
    for t = 1:(analysisSettings.finalTime/analysisSettings.deltaT)
        x2(t)=Rstrain(t,2);
        y2(t)=Rstress(t,2);
    end
        figure(2);
        plot(x1,y1,'b','LineWidth',1.5);
        title('Plasticidad / Barra 1 --algoritmo de retorno-- \sigma(\epsilon)');
        ylabel('Tensión \sigma');
        xlabel('Deformación \epsilon');
        hold on;
        figure(3);
        plot(x2,y2,'r','LineWidth',1.5);
        title('Plasticidad / Barra 2 --algoritmo de retorno-- \sigma(\epsilon)');
        ylabel('Tensión \sigma');
        xlabel('Deformación \epsilon');
        hold on;
        
    for t = 1:(analysisSettings.finalTime/analysisSettings.deltaT)
        x1(t)=t;
        y1(t)=Rstrainacum(t,1);
    end
    for t = 1:(analysisSettings.finalTime/analysisSettings.deltaT)
        x2(t)=t;
        y2(t)=Rstrainacum(t,2);
    end
        figure(4);
        plot(x1,y1,'b','LineWidth',1.5);
        title('Deformación Plástica Acumulada / Barra 1 \epsilon_p acumulada');
        ylabel('Deformación acumulada\epsilon_p');
        xlabel('Tiempo');
        hold on;
        figure(5);
        plot(x2,y2,'r','LineWidth',1.5);
        title('Deformación Plástica Acumulada / Barra 2 \epsilon_p acumulada');
        ylabel('Deformación acumulada\epsilon_p');
        xlabel('Tiempo');
        hold on;

fprintf('Desplazamiento vertical:\n');
fprintf('\n');
disp(deltas(length(deltas),2));
fprintf('\n');

x = [0 (1-deltas(length(deltas),1))/1.05 2]; % se amplifica el desplazamiento
y = [0 (1-deltas(length(deltas),2))/1.05 0]; % para su visualización
z = [0 0 0];        
s = [1 3];
t = [2 2];
weights = [2 2];
G = graph(s,t,weights);
figure(1);
p=plot(G,'XData',x,'YData',y,'linewidth',2,'EdgeColor','r');
p.NodeColor = 'r';
p.MarkerSize = 5;
p.LineStyle = '--';
title('Plasticidad / Reticulado de von Mises');
hold on;

%{
x = [0 (1+deltas(length(deltas),1))/1.05 2]; % se amplifica el desplazamiento
y = [0 (1+deltas(length(deltas),2))/1.05 0]; % para su visualización
z = [0 0 0];       
figure(1);
plot(x,y,'ro--','markersize',4,'markerfacecolor','b','linewidth',2);
title('Plasticidad / Reticulado de von Mises');
hold on;
%}