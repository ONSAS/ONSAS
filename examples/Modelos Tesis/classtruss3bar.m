%md# Classical Truss Three Rods

close all; clear; addpath( genpath( [ pwd '/../../src'] ) );
format long;
% scalar parameters (N/cm2)
E = 210e3 ;
Kplas = 1093 ;
sigma_Y_0 = 330.3 ;
Fu = 1 ;

global Rstress Rstrain Rstrainacum Largo ;

% x and z coordinates of node 4
x2 = 300 ;
z2 = -400 ;

materials.hyperElasModel  = 'isotropicHardening' ;
materials.hyperElasParams = [ E Kplas sigma_Y_0 ] ;

elements(1).elemType = 'node' ;
elements(2).elemType = 'truss';
elements(3).elemType = 'truss';
elements(4).elemType = 'truss';
elements(2).elemCrossSecParams = { 'circle' , sqrt(1*4/pi)} ;
elements(3).elemCrossSecParams = { 'circle' , sqrt(1*4/pi)} ;
elements(4).elemCrossSecParams = { 'circle' , sqrt(1*4/pi)} ;

boundaryConds(1).imposDispDofs = [ 1 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 ] ;
boundaryConds(2).imposDispDofs = 3 ;
boundaryConds(2).imposDispVals = 0 ;
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsTimeFact = @(t) (t<=1)*(Fu)*t ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 1 0 ]  ;

initialConds = struct() ;

mesh.nodesCoords = [   0  0   0 ; ...
                      x2  0   0 ; ...
                    2*x2  0   0 ; ...
                      x2  0   z2 ] ;

% MEBI [Material Element Boundary_Conditions Initial_Conditions]

mesh.conecCell = cell(7,1) ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 1 0  2   ] ;
mesh.conecCell{ 3, 1 } = [ 0 1 1 0  3   ] ;
mesh.conecCell{ 4, 1 } = [ 0 1 2 0  4   ] ;
mesh.conecCell{ 5, 1 } = [ 1 2 0 0  1 4 ] ;
mesh.conecCell{ 6, 1 } = [ 1 3 0 0  2 4 ] ;
mesh.conecCell{ 7, 1 } = [ 1 4 0 0  3 4 ] ;

analysisSettings.deltaT        =      1 ;
analysisSettings.finalTime     =    100 ;
analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   15   ;

otherParams.problemName       = 'ClassicalTruss_log_strain_arcLength' ;
analysisSettings.methodName   = 'arcLength' ;

analysisSettings.incremArcLen = 0.064;
analysisSettings.iniDeltaLamb = boundaryConds(2).loadsTimeFact(1)/100 ;
analysisSettings.posVariableLoadBC = 2 ;

otherParams.problemName = 'classtruss' ;
otherParams.plots_format = 'vtk' ;
otherParams.plots_deltaTs_separation = 2 ;

[ matUs, loadFactorsMat ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

deltas = [-matUs(6*3+1,:)' -matUs(6*3+5,:)'] ;

  x1=zeros(floor(analysisSettings.finalTime/analysisSettings.deltaT),1);
  y1=zeros(floor(analysisSettings.finalTime/analysisSettings.deltaT),1);
  F=zeros(floor(analysisSettings.finalTime/analysisSettings.deltaT),1);
    for t = 1:(floor(analysisSettings.finalTime/analysisSettings.deltaT))
        x1(t)=Rstrain(t,1);
        y1(t)=Rstress(t,1);
        % Fuerza aplicada vertical
        F(t)=(z2-deltas(t,2))./Largo(t,2).*1.*y1(t);
    end
    x2=zeros(floor(analysisSettings.finalTime/analysisSettings.deltaT),1);
    y2=zeros(floor(analysisSettings.finalTime/analysisSettings.deltaT),1);
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

fprintf('|Desplazamiento vertical\n');
fprintf('|Solución numérica  %.15f\n',abs(deltas(length(deltas),2)));
fprintf('|Deformación plástica acumulada (alfa) %.15f\n',Rstrainacum(1));

% an Eternal Golden Braid