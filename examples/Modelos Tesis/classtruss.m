%md# Classical Truss

close all; clear; addpath( genpath( [ pwd '/../../src'] ) );
% scalar parameters (N/cm2)
E = 210e3 ;
Kplas = 529.5 ;
sigma_Y_0 = 123.6 ;
Fu = -1 ;

global Rstress Rstrain Rstrainacum Largo ;

% x and z coordinates of node 2
x2 = 220 ;
z2 = 127 ;

materials.hyperElasModel  = 'isotropicHardening' ;
materials.hyperElasParams = [ E Kplas sigma_Y_0 ] ;

elements(1).elemType = 'node' ;
elements(2).elemType = 'truss';
elements(3).elemType = 'truss';
elements(2).elemCrossSecParams = { 'circle' , sqrt(7*4/pi)} ;
elements(3).elemCrossSecParams = { 'circle' , sqrt(7*4/pi)} ;

boundaryConds(1).imposDispDofs = [ 1 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 ] ;
boundaryConds(2).imposDispDofs = 3 ;
boundaryConds(2).imposDispVals = 0 ;
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsTimeFact = @(t) Fu*t ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 1 0 ]  ;

initialConds = struct() ;

mesh.nodesCoords = [   0  0   0 ; ...
                      x2  0  z2 ; ...
                    2*x2  0   0 ] ;

% MEBI [Material Element Boundary_Conditions Initial_Conditions]

mesh.conecCell = cell(5,1) ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 1 3   ] ;
mesh.conecCell{ 3, 1 } = [ 0 1 2 2   ] ;
mesh.conecCell{ 4, 1 } = [ 1 2 0 1 2 ] ;
mesh.conecCell{ 5, 1 } = [ 1 3 0 2 3 ] ;

analysisSettings.deltaT        =       1 ;
analysisSettings.finalTime     =    1000 ;
analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   15   ;

otherParams.problemName       = 'ClassicalTruss_log_strain_arcLength' ;
analysisSettings.methodName   = 'arcLength' ;

analysisSettings.incremArcLen = 0.3 ;
analysisSettings.iniDeltaLamb = boundaryConds(2).loadsTimeFact(1)/1000 ;
analysisSettings.posVariableLoadBC = 2 ;

otherParams.problemName = 'classtruss';
otherParams.plots_format = 'vtk' ;
otherParams.plots_deltaTs_separation = 2 ;

[matUs, loadFactorsMat ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

deltas = [-matUs(6+1,:)' -matUs(6+5,:)'] ;

    for t = 1:(floor(analysisSettings.finalTime/analysisSettings.deltaT))
        x1(t)=Rstrain(t,1);
        y1(t)=Rstress(t,1);
        % Fuerza aplicada vertical
        F(t)=(z2-deltas(t,2))./Largo(t,2).*7.*y1(t);
    end

    for t = 1:(floor(analysisSettings.finalTime/analysisSettings.deltaT))
        x2(t)=Rstrain(t,2);
        y2(t)=Rstress(t,2);
    end

        figure(1);
        plot(-deltas(2:end,2),F,'LineWidth',1.5);
        title('Plasticidad / Barras --ONSAS-- F_{apl}(\delta)');
        ylabel('Fuerza Vertical F_{apl}');
        xlabel('Desplazamiento \delta');
        hold on;
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        figure(2);
        plot(x1,-deltas(2:end,2),'LineWidth',1.5);
        title('Plasticidad / Barras --ONSAS-- \sigma(\epsilon)');
        ylabel('Desplazamiento \delta');
        xlabel('Deformación \epsilon');
        hold on;
        figure(3);
        plot(x1,y1,'LineWidth',1.5);
        title('Plasticidad / Barras --ONSAS-- \sigma(\epsilon)');
        ylabel('Tensión \sigma');
        xlabel('Deformación \epsilon');
        hold on;
        fprintf('\n');
%{
% solución analítica
syms dw1;
syms dw2;
z2=1;
Ae=pi*0.010^2/4;
eqn1= (2*z2*dw1 + (dw1)^2)/(2*((dw1+z2)^2+(z2)^2)) == (sigma_Y_0)/E + (E+Kplas)/(E*Kplas)*(Fu/2/(z2+dw1)/Ae*sqrt((dw1+z2)^2+(z2)^2)-sigma_Y_0);
S1 = vpasolve(eqn1,dw1);
eqn2= (2*(z2+S1(1))*dw2 + (dw2)^2)/(2*((dw2+z2+S1(1))^2+(z2+S1(1))^2)) == Fu/2/(z2+S1(1)+dw2)/Ae*sqrt((S1(1)+z2+dw2)^2+(z2+S1(1))^2)/E;
S2 = vpasolve(eqn2,dw2);
%}
fprintf('|Desplazamiento vertical\n');
% fprintf('|Solución analítica %.15f\n',S1);
fprintf('|Solución numérica  %.15f\n',abs(deltas(length(deltas),2)));
% fprintf('\n');
%}
% fprintf('| Deformación plástica acumulada alfa\n');
% alfaa=(E+Kplas)/(E*Kplas)*(Fu/2/(z2+S1(1))/Ae*sqrt((S1(1)+z2)^2+(z2)^2)-sigma_Y_0);
fprintf('|Deformación plástica acumulada (alfa) %.15f\n',Rstrainacum(1));
% scalar parameters (N/cm2)
E = 210e3 ;
Kplas = 1093 ;
sigma_Y_0 = 330.3 ;
Fu = -1 ;

global Rstress Rstrain Rstrainacum Largo ;

% x and z coordinates of node 2
x2 = 220 ;
z2 = 127 ;

materials.hyperElasModel  = 'isotropicHardening' ;
materials.hyperElasParams = [ E Kplas sigma_Y_0 ] ;

elements(1).elemType = 'node' ;
elements(2).elemType = 'truss';
elements(3).elemType = 'truss';
elements(2).elemCrossSecParams = { 'circle' , sqrt(7*4/pi)} ;
elements(3).elemCrossSecParams = { 'circle' , sqrt(7*4/pi)} ;

boundaryConds(1).imposDispDofs = [ 1 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 ] ;
boundaryConds(2).imposDispDofs = 3 ;
boundaryConds(2).imposDispVals = 0 ;
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsTimeFact = @(t) (t<=1)*Fu*t ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 1 0 ]  ;

initialConds = struct() ;

mesh.nodesCoords = [   0  0   0 ; ...
                      x2  0  z2 ; ...
                    2*x2  0   0 ] ;

mesh.conecCell = cell(5,1) ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 1 3   ] ;
mesh.conecCell{ 3, 1 } = [ 0 1 2 2   ] ;
mesh.conecCell{ 4, 1 } = [ 1 2 0 1 2 ] ;
mesh.conecCell{ 5, 1 } = [ 1 3 0 2 3 ] ;

analysisSettings.deltaT        =       1 ;
analysisSettings.finalTime     =    1000 ;
analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   15   ;

otherParams.problemName       = 'ClassicalTruss_log_strain_arcLength' ;
analysisSettings.methodName   = 'arcLength' ;

analysisSettings.incremArcLen = 0.3 ;
analysisSettings.iniDeltaLamb = boundaryConds(2).loadsTimeFact(1)/1000 ;
analysisSettings.posVariableLoadBC = 2 ;

otherParams.problemName = 'classtruss';
otherParams.plots_format = 'vtk' ;
otherParams.plots_deltaTs_separation = 2 ;

[matUs, loadFactorsMat ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

deltas = [-matUs(6+1,:)' -matUs(6+5,:)'] ;

    for t = 1:(floor(analysisSettings.finalTime/analysisSettings.deltaT))
        x1(t)=Rstrain(t,1);
        y1(t)=Rstress(t,1);
        % Fuerza aplicada vertical
        F(t)=(z2-deltas(t,2))./Largo(t,2).*7.*y1(t);
    end

    for t = 1:(floor(analysisSettings.finalTime/analysisSettings.deltaT))
        x2(t)=Rstrain(t,2);
        y2(t)=Rstress(t,2);
    end

        figure(1);
        plot(-deltas(2:end,2),F,'LineWidth',1.5);
        title('Plasticidad / Barras --ONSAS-- F_{apl}(\delta)');
        ylabel('Fuerza Vertical F_{apl}');
        xlabel('Desplazamiento \delta');
        hold on;
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        figure(2);
        plot(x1,-deltas(2:end,2),'LineWidth',1.5);
        title('Plasticidad / Barras --ONSAS-- \sigma(\epsilon)');
        ylabel('Desplazamiento \delta');
        xlabel('Deformación \epsilon');
        hold on;
        figure(3);
        plot(x1,y1,'LineWidth',1.5);
        title('Plasticidad / Barras --ONSAS-- \sigma(\epsilon)');
        ylabel('Tensión \sigma');
        xlabel('Deformación \epsilon');
        hold on;
        fprintf('\n');
%{
% solución analítica
syms dw1;
syms dw2;
z2=1;
Ae=pi*0.010^2/4;
eqn1= (2*z2*dw1 + (dw1)^2)/(2*((dw1+z2)^2+(z2)^2)) == (sigma_Y_0)/E + (E+Kplas)/(E*Kplas)*(Fu/2/(z2+dw1)/Ae*sqrt((dw1+z2)^2+(z2)^2)-sigma_Y_0);
S1 = vpasolve(eqn1,dw1);
eqn2= (2*(z2+S1(1))*dw2 + (dw2)^2)/(2*((dw2+z2+S1(1))^2+(z2+S1(1))^2)) == Fu/2/(z2+S1(1)+dw2)/Ae*sqrt((S1(1)+z2+dw2)^2+(z2+S1(1))^2)/E;
S2 = vpasolve(eqn2,dw2);
%}
fprintf('|Desplazamiento vertical\n');
% fprintf('|Solución analítica %.15f\n',S1);
fprintf('|Solución numérica  %.15f\n',abs(deltas(length(deltas),2)));
% fprintf('\n');
%}
% fprintf('| Deformación plástica acumulada alfa\n');
% alfaa=(E+Kplas)/(E*Kplas)*(Fu/2/(z2+S1(1))/Ae*sqrt((S1(1)+z2)^2+(z2)^2)-sigma_Y_0);
fprintf('|Deformación plástica acumulada (alfa) %.15f\n',Rstrainacum(1));
% scalar parameters (N/cm2)
E = 206e3 ;
Kplas = 1545.1 ;
sigma_Y_0 = 258.3 ;
Fu = -1 ;

global Rstress Rstrain Rstrainacum Largo ;

% x and z coordinates of node 2
x2 = 220 ;
z2 = 127 ;

materials.hyperElasModel  = 'isotropicHardening' ;
materials.hyperElasParams = [ E Kplas sigma_Y_0 ] ;

elements(1).elemType = 'node' ;
elements(2).elemType = 'truss';
elements(3).elemType = 'truss';
elements(2).elemCrossSecParams = { 'circle' , sqrt(7*4/pi)} ;
elements(3).elemCrossSecParams = { 'circle' , sqrt(7*4/pi)} ;

boundaryConds(1).imposDispDofs = [ 1 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 ] ;
boundaryConds(2).imposDispDofs = 3 ;
boundaryConds(2).imposDispVals = 0 ;
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsTimeFact = @(t) (t<=1)*Fu*t ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 1 0 ]  ;

initialConds = struct() ;

mesh.nodesCoords = [   0  0   0 ; ...
                      x2  0  z2 ; ...
                    2*x2  0   0 ] ;

% MEBI [Material Element Boundary_Conditions Initial_Conditions]

mesh.conecCell = cell(5,1) ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 1 3   ] ;
mesh.conecCell{ 3, 1 } = [ 0 1 2 2   ] ;
mesh.conecCell{ 4, 1 } = [ 1 2 0 1 2 ] ;
mesh.conecCell{ 5, 1 } = [ 1 3 0 2 3 ] ;

analysisSettings.deltaT        =       1 ;
analysisSettings.finalTime     =    1000 ;
analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   15   ;

otherParams.problemName       = 'ClassicalTruss_log_strain_arcLength' ;
analysisSettings.methodName   = 'arcLength' ;

analysisSettings.incremArcLen = 0.3 ;
analysisSettings.iniDeltaLamb = boundaryConds(2).loadsTimeFact(1)/1000 ;
analysisSettings.posVariableLoadBC = 2 ;

otherParams.problemName = 'classtruss';
otherParams.plots_format = 'vtk' ;
otherParams.plots_deltaTs_separation = 2 ;

[matUs, loadFactorsMat ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

deltas = [-matUs(6+1,:)' -matUs(6+5,:)'] ;

    for t = 1:(floor(analysisSettings.finalTime/analysisSettings.deltaT))
        x1(t)=Rstrain(t,1);
        y1(t)=Rstress(t,1);
        % Fuerza aplicada vertical
        F(t)=(z2-deltas(t,2))./Largo(t,2).*7.*y1(t);
    end

    for t = 1:(floor(analysisSettings.finalTime/analysisSettings.deltaT))
        x2(t)=Rstrain(t,2);
        y2(t)=Rstress(t,2);
    end

        figure(1);
        plot(-deltas(2:end,2),F,'LineWidth',1.5);
        title('Plasticidad / Barras --ONSAS-- F_{apl}(\delta)');
        ylabel('Fuerza Vertical F_{apl}');
        xlabel('Desplazamiento \delta');
        hold on;
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        figure(2);
        plot(x1,-deltas(2:end,2),'LineWidth',1.5);
        title('Plasticidad / Barras --ONSAS-- \sigma(\epsilon)');
        ylabel('Desplazamiento \delta');
        xlabel('Deformación \epsilon');
        hold on;
        figure(3);
        plot(x1,y1,'LineWidth',1.5);
        title('Plasticidad / Barras --ONSAS-- \sigma(\epsilon)');
        ylabel('Tensión \sigma');
        xlabel('Deformación \epsilon');
        hold on;
        fprintf('\n');
%{
% solución analítica
syms dw1;
syms dw2;
z2=1;
Ae=pi*0.010^2/4;
eqn1= (2*z2*dw1 + (dw1)^2)/(2*((dw1+z2)^2+(z2)^2)) == (sigma_Y_0)/E + (E+Kplas)/(E*Kplas)*(Fu/2/(z2+dw1)/Ae*sqrt((dw1+z2)^2+(z2)^2)-sigma_Y_0);
S1 = vpasolve(eqn1,dw1);
eqn2= (2*(z2+S1(1))*dw2 + (dw2)^2)/(2*((dw2+z2+S1(1))^2+(z2+S1(1))^2)) == Fu/2/(z2+S1(1)+dw2)/Ae*sqrt((S1(1)+z2+dw2)^2+(z2+S1(1))^2)/E;
S2 = vpasolve(eqn2,dw2);
%}
fprintf('|Desplazamiento vertical\n');
% fprintf('|Solución analítica %.15f\n',S1);
fprintf('|Solución numérica  %.15f\n',abs(deltas(length(deltas),2)));
% fprintf('\n');
%}
% fprintf('| Deformación plástica acumulada alfa\n');
% alfaa=(E+Kplas)/(E*Kplas)*(Fu/2/(z2+S1(1))/Ae*sqrt((S1(1)+z2)^2+(z2)^2)-sigma_Y_0);
fprintf('|Deformación plástica acumulada (alfa) %.15f\n',Rstrainacum(1));