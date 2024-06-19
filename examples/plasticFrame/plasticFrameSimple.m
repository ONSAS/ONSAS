
% =========================================================================

% Elastoplastic analysis of plane frame
% Elements with embedded discontinuity
% Softening hinges

% =========================================================================

close all ; clear ;
addpath( genpath( [ pwd '/../../src'] ) ) ;

% assumed XY plane
% geometry
l = 3           ;   % m
Inertia = 1e-3  ;   % m^4
E = 28.6e3*1e3  ;   % KN/m^2 [KPa]
EI = E*Inertia  ;   % KN.m^2
A  = 0.10       ;   % m^2

% material
kh1 = 12450 ;       % KN.m^2
kh2 = 195   ;
Ks  = -2410 ;       % KN.m
nu = 0.3 ;          % Poisson's ratio

Mc = 100 ;
My = 245 ;
Mu = 265 ;          % KN.m

% /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\

otherParams              = struct() ;
otherParams.problemName  = 'plastic_2dframe' ;
% otherParams.plots_format = 'vtk' ;

materials             = struct() ;
materials.modelName   = 'plastic-2Dframe' ;
materials.modelParams = [ E Mc My Mu kh1 kh2 Ks nu ] ;

elements             = struct() ;
elements(1).elemType = 'node' ;

elements(2).elemType = 'frame' ;

elements(2).elemCrossSecParams = {'generic' ; [A 1 Inertia Inertia] } ;

boundaryConds                  = {} ;
% supports
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;

boundaryConds(2).imposDispDofs = [ 2 4 5 ] ;
boundaryConds(2).imposDispVals = [ 0 0 0 ] ;
% loads
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsBaseVals = [ 1 0 0 0 0 0 ] ;
boundaryConds(2).loadsTimeFact = @(t) t ;

boundaryConds(3).imposDispDofs = [ 2 4 5 ] ;
boundaryConds(3).imposDispVals = [ 0 0 0 ] ;


% The coordinates of the nodes of the mesh are given by the matrix:
mesh = {} ;
mesh.nodesCoords = [0  0       0    ; ...
					0  l       0    ; ...
					l  l       0]   ;

mesh.conecCell = {} ;

% nodes
mesh.conecCell{ 1, 1 }      = [ 0 1 1 1 ] ; % node 1 fixed end support
mesh.conecCell{ end+1, 1 }  = [ 0 1 3 2 ] ; % node 2
mesh.conecCell{ end+1, 1 }  = [ 0 1 2 3 ] ; % node 3 with load applied

% frame elements
mesh.conecCell{ end+1, 1 } = [ 1 2 0   1 2] ;
mesh.conecCell{ end+1, 1 } = [ 1 2 0   2 3] ;

initialConds = {} ;

analysisSettings                    = {} ;
analysisSettings.methodName         = 'arcLength' ;
analysisSettings.deltaT             = 1 ;
analysisSettings.incremArcLen       = 1e-3*ones(1,1000) ;
analysisSettings.finalTime          = length(analysisSettings.incremArcLen) ;
analysisSettings.iniDeltaLamb       = 1 ;
analysisSettings.posVariableLoadBC  = 2 ;
analysisSettings.stopTolDeltau      = 1e-14 ;
analysisSettings.stopTolForces      = 1e-8 ;
analysisSettings.stopTolIts         = 30 ;
analysisSettings.ALdominantDOF      = [2*6+1 -1] ;

[ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

[matUs, loadFactorsMat, modelSolutions ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;

rotations = matUs((2)*6+6,:) ;
displacements = matUs((2)*6+1,:) ; % node with vertical load applied
loadfactors = loadFactorsMat(:,2) ;

moments_hist = zeros(4,length(modelSolutions)) ;
for i =1:length(modelSolutions)
    aux = modelSolutions{i}.localInternalForces(1) ;
    moments_hist(:,i) = [ aux.Mz; aux.Mz2; aux.Mz3; aux.tM ] ;
end
Mn1_numericONSAS = moments_hist(1,:) ;
Mn2_numericONSAS = moments_hist(2,:) ;
Mn3_numericONSAS = moments_hist(3,:) ;
tMn_numericONSAS = moments_hist(4,:) ;

% plots

lw = 2 ; ms = 1 ; plotfontsize = 14 ;

figure('Name','Frame / Plasticity (load factors)','NumberTitle','off') ;
hold on, grid on

plot(abs(rotations), abs(loadfactors), '-x', 'linewidth', lw, 'markersize', ms, "Color", "#EDB120") ;
plot(abs(displacements), abs(loadfactors), '-x', 'linewidth', lw, 'markersize', ms, "Color", "#0072BD") ;

labx = xlabel('Generalized displacements in free node (m, rad)') ;
laby = ylabel('Forces') ;

legend('ONSAS [\theta]', 'ONSAS [y]', 'location', 'Southeast') ;

set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Frame / Plasticity (load factors)') ;

figure('Name','Frame / Plasticity (Moments)','NumberTitle','off') ;
hold on, grid on

plot(abs(displacements), abs(tMn_numericONSAS), '-x', 'linewidth', lw, 'markersize', ms, "Color", "#0072BD") ;

labx = xlabel('Generalized displacements in free node (m, rad)') ;
laby = ylabel('Moments tMn') ;

legend('ONSAS tMn [y]', 'location', 'Southeast') ;

set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Frame / Plasticity (load factors)') ;

print('-f1','../../../Tesis/tex/imagenes/DarvallMendisFrameLoadFactors.pdf','-dpdf') ;