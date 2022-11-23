%md# Static von Mises Truss example

close all; clear; addpath( genpath( [ pwd '/../../src'] ) );
format long;
% scalar parameters
E = 210e9 ;
Kplas = 21e8;
sigma_Y_0 = 250e6 ;

global Rstress Rstrain Rstrainacum ;

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

boundaryConds(1).imposDispDofs = [ 1 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 ] ;

boundaryConds(2).imposDispDofs =  3 ;
boundaryConds(2).imposDispVals =  0 ;
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsTimeFact = @(t) (t<1)*28000*t+(t>=1 & t<=2)*(28000-56300*(t-1))+ (t>2 & t<=3)*(-28300 + 57300*(t-2)) + (t>3 & t<=4)*(29000 - 58000*(t-3));
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 1 0 ]  ;

% Rango=[0,28000,28000,-28300,-28300,29000,29000,-29000];

initialConds = struct() ;

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
analysisSettings.deltaT        =   1/140 ;
analysisSettings.finalTime     =   4 ;

analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   15   ;

analysisSettings.posVariableLoadBC = 2 ;

otherParams.problemName = 'static_plastic_von_mises_truss';
otherParams.plots_format = 'vtk' ;
otherParams.plots_deltaTs_separation = 2 ;

[matUs, loadFactorsMat ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

deltas = [-matUs(6+1,:)' -matUs(6+5,:)'] ;

    for t = 1:(floor(analysisSettings.finalTime/analysisSettings.deltaT))
        x1(t)=Rstrain(t,1);
        y1(t)=Rstress(t,1);
    end

    for t = 1:(floor(analysisSettings.finalTime/analysisSettings.deltaT))
        x2(t)=Rstrain(t,2);
        y2(t)=Rstress(t,2);
    end
        figure(2);
        plot(x1,y1,'LineWidth',1.5);
        title('Plasticidad / Barra 1 --ONSAS-- \sigma(\epsilon)');
        ylabel('Tensión \sigma');
        xlabel('Deformación \epsilon');
        hold on;
        figure(3);
        plot(x2,y2,'LineWidth',1.5);
        title('Plasticidad / Barra 2 --ONSAS-- \sigma(\epsilon)');
        ylabel('Tensión \sigma');
        xlabel('Deformación \epsilon');
        hold on;
        fprintf('\n');

% solución analítica
%{
syms dw1;
syms dw2;
z2=1;
Ae=pi*0.010^2/4;
eqn1= (2*z2*dw1 + (dw1)^2)/(2*((dw1+z2)^2+(z2)^2)) == (sigma_Y_0)/E + (E+Kplas)/(E*Kplas)*(Fu/2/(z2+dw1)/Ae*sqrt((dw1+z2)^2+(z2)^2)-sigma_Y_0);
S1 = vpasolve(eqn1,dw1);
eqn2= (2*(z2+S1(1))*dw2 + (dw2)^2)/(2*((dw2+z2+S1(1))^2+(z2+S1(1))^2)) == Fu/2/(z2+S1(1)+dw2)/Ae*sqrt((S1(1)+z2+dw2)^2+(z2+S1(1))^2)/E;
S2 = vpasolve(eqn2,dw2);
%}
% fprintf('| Desplazamiento vertical\n');
% fprintf(' |Solución analítica %.15f\n',S1(1));
fprintf('| Desplazamiento vertical               %.15f\n',abs(deltas(length(deltas),2)));
% fprintf('\n');
%}
% fprintf('| Deformación plástica acumulada alfa\n');
% alfaa=(E+Kplas)/(E*Kplas)*(Fu/2/(z2+S1(1))/Ae*sqrt((S1(1)+z2)^2+(z2)^2)-sigma_Y_0);
fprintf('| Deformación plástica acumulada (alfa) %.15f\n',Rstrainacum(1));