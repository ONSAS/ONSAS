% =========================================================================

% Elastoplastic analysis of Darvall-Mendis Frame
% Elements with embedded discontinuity
% Softening hinges

% Elastic-Plastic-Softening Analysis of Plane Frames
% Peter LePoer Darvall and Priyantha Anumddha Mendis
% Dept. of Civ. Engrg., Monash Univ., Victoria, Australia
% Journal of Structural Engineering, Vol. III, No. 4, April, 1985.

% =========================================================================

close all ; clear ;
addpath( genpath( [ pwd '/../../src'] ) ) ;

% assumed XY plane
% geometry
l = 3           ;   % m
Inertia = 1*1e3 ;   % m^4
E = 20.7*1e6    ;   % KN/m^2 [KPa]
EI = E*Inertia  ;   % KN.m^2
A  = 0.10       ;   % m^2

% material

kh1 = 0 ;           % KN.m^2
kh2 = 0 ;
Ks  = -1089 ;       % KN.m
nu = 0.3 ;          % Poisson's ratio

Mc = 100 ;
My = 100 ;
Mu = 100 ;          % KN.m

% /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\

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

boundaryConds(3).imposDispDofs = [ 2 4 5 ] ;
boundaryConds(3).imposDispVals = [ 0 0 0 ] ;

% loads
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsBaseVals = [ 0 0 -1 0 0 0 ] ;
boundaryConds(2).loadsTimeFact = @(t) t ;

% The coordinates of the nodes of the mesh are given by the matrix:
mesh = {} ;
mesh.nodesCoords = [0  0       0    ; ...
					0  0.5*l   0    ; ...
					0  l       0    ; ...
              0.275*l  l       0    ; ...
               0.55*l  l       0    ; ...
              0.775*l  l       0    ; ...
                    l  l       0    ; ...
                    l  0.5*l   0    ; ...
					l  0       0]   ;

mesh.conecCell = {} ;

% nodes
mesh.conecCell{ 1, 1 } = [ 0 1 1 1 ] ; % node 1 fixed end support
mesh.conecCell{ 2, 1 } = [ 0 1 3 2 ] ; % node 2
mesh.conecCell{ 3, 1 } = [ 0 1 3 3 ] ; % node 3
mesh.conecCell{ 4, 1 } = [ 0 1 3 4 ] ; % node 4
mesh.conecCell{ 5, 1 } = [ 0 1 2 5 ] ; % node 5 with vertical load applied
mesh.conecCell{ 6, 1 } = [ 0 1 3 6 ] ; % node 6
mesh.conecCell{ 7, 1 } = [ 0 1 3 7 ] ; % node 7
mesh.conecCell{ 8, 1 } = [ 0 1 3 8 ] ; % node 8
mesh.conecCell{ 9, 1 } = [ 0 1 1 9 ] ; % node 9 fixed end support

% frame elements
mesh.conecCell{ 10, 1 } = [ 1 2 0   1 2] ;
mesh.conecCell{ 11, 1 } = [ 1 2 0   2 3] ;
mesh.conecCell{ 12, 1 } = [ 1 2 0   3 4] ;
mesh.conecCell{ 13, 1 } = [ 1 2 0   4 5] ;
mesh.conecCell{ 14, 1 } = [ 1 2 0   5 6] ;
mesh.conecCell{ 15, 1 } = [ 1 2 0   6 7] ;
mesh.conecCell{ 16, 1 } = [ 1 2 0   7 8] ;
mesh.conecCell{ 17, 1 } = [ 1 2 0   8 9] ;

initialConds = {} ;

analysisSettings                    = {} ;
analysisSettings.methodName         = 'arcLength' ;
analysisSettings.deltaT             = 1 ;
analysisSettings.incremArcLen       = 1e-4*ones(1,100) ;
analysisSettings.finalTime          = length(analysisSettings.incremArcLen) ;
analysisSettings.iniDeltaLamb       = 1 ;
analysisSettings.posVariableLoadBC  = 2 ;
analysisSettings.stopTolDeltau      = 1e-14 ;
analysisSettings.stopTolForces      = 1e-8 ;
analysisSettings.stopTolIts         = 30 ;
analysisSettings.ALdominantDOF      = [2*6-3 -1] ;

otherParams              = struct() ;
otherParams.problemName  = 'plastic_2dframe' ;
% otherParams.plots_format = 'vtk' ;

[ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

[matUs, loadFactorsMat, ~ ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;

rotations = matUs((4+1)*6,:) ;
displacements = matUs((4+1)*6-3,:) ; % node with vertical load applied
loadfactors = loadFactorsMat(:,2) ;

% GRAPHICS

lw = 2 ; ms = 1 ; plotfontsize = 14 ;

figure('Name','Darvall-Mendis Frame / Plasticity (load factors)','NumberTitle','off') ;
hold on, grid on

plot(abs(rotations), loadfactors, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#EDB120") ;
plot(abs(displacements), loadfactors, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#0072BD") ;

labx = xlabel('Generalized displacements in free node (m, rad)') ;
laby = ylabel('Forces') ;

legend('ONSAS (8 elem) [\theta]', 'ONSAS (8 elem) [y]', 'location', 'Southeast') ;

set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Darvall-Mendis Frame / Plasticity (load factors)') ;

print('-f1','../../../Tesis/tex/imagenes/DarvallMendisFrameLoadFactors.pdf','-dpdf') ;