% =========================================================================

% Elastoplastic analysis of Darvall-Mendis Frame
% Elements with embedded discontinuity
% Softening hinges

% Elastic-Plastic-Softening Analysis of Plane Frames
% Peter LePoer Darvall and Priyantha Anumddha Mendis
% Dept. of Civ. Engrg., Monash Univ., Victoria, Australia
% Journal of Structural Engineering, Vol. III, No. 4, April, 1985

% Numerical modeling of softening hinges in thin Euler–Bernoulli beams
% Francisco Armero, David Ehrlich / University of California, Berkeley

% =========================================================================

close all ; clear ;
addpath( genpath( [ pwd '/../../src'] ) ) ;

% assumed XY plane
% geometry
l = 3.048       ;   % m
Inertia = 1e-3  ;   % m^4
E = 2.068e7     ;   % KN/m^2 [KPa]
EI = E*Inertia  ;   % KN.m^2
A  = 0.103      ;   % m^2

% material
kh1 = 1 ;       % KN.m^2
kh2 = 1 ;

% Ks  = -1e-6           ;   % almost zero

a  = -0.06              ;
Ks = a*EI/10/l          ;   % KN.m

nu  = 0.3               ;   % Poisson's ratio

Mc_columns = 158.06 ;
My_columns = 158.12 ;

Mc_beams = 169.16 ;
My_beams = 169.32 ;

Mu_columns = 158.18    ;   % KN.m
Mu_beams   = 169.48    ;   % KN.m

% /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\

materials             = struct() ;
materials(1).modelName   = 'plastic-2Dframe' ;
materials(1).modelParams = [ E Mc_beams My_beams Mu_beams kh1 kh2 Ks nu ] ;

materials(2).modelName   = 'plastic-2Dframe' ;
materials(2).modelParams = [ E Mc_columns My_columns Mu_columns kh1 kh2 Ks nu ] ;

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
boundaryConds(2).loadsBaseVals = [ 0 0 -1 0 0 0 ] ;
boundaryConds(2).loadsTimeFact = @(t) t ;

boundaryConds(3).imposDispDofs = [ 2 4 5 ] ;
boundaryConds(3).imposDispVals = [ 0 0 0 ] ;

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
mesh.conecCell{ 1, 1 }     = [ 0 1 1 1 ] ; % node 1 fixed end support
mesh.conecCell{ 9, 1 }     = [ 0 1 1 9 ] ; % node 9 fixed end support

mesh.conecCell{ 2, 1 } = [ 0 1 3 2 ] ; % node 2
mesh.conecCell{ 3, 1 } = [ 0 1 3 3 ] ; % node 3
mesh.conecCell{ 4, 1 } = [ 0 1 3 4 ] ; % node 4
mesh.conecCell{ 5, 1 } = [ 0 1 2 5 ] ; % node 5 with vertical load applied
mesh.conecCell{ 6, 1 } = [ 0 1 3 6 ] ; % node 6
mesh.conecCell{ 7, 1 } = [ 0 1 3 7 ] ; % node 7
mesh.conecCell{ 8, 1 } = [ 0 1 3 8 ] ; % node 8

% frame elements
mesh.conecCell{ 10, 1 } = [ 2 2 0   1 2] ;
mesh.conecCell{ 11, 1 } = [ 2 2 0   2 3] ;
mesh.conecCell{ 12, 1 } = [ 1 2 0   3 4] ;
mesh.conecCell{ 13, 1 } = [ 1 2 0   4 5] ;
mesh.conecCell{ 14, 1 } = [ 1 2 0   5 6] ;
mesh.conecCell{ 15, 1 } = [ 1 2 0   6 7] ;
mesh.conecCell{ 16, 1 } = [ 2 2 0   7 8] ;
mesh.conecCell{ 17, 1 } = [ 2 2 0   8 9] ;

initialConds = {} ;

analysisSettings                    = {} ;
analysisSettings.methodName         = 'arcLength' ;
analysisSettings.deltaT             = 1 ;
analysisSettings.incremArcLen       = 1e-4*ones(1,500);
analysisSettings.finalTime          = length(analysisSettings.incremArcLen) ;
analysisSettings.iniDeltaLamb       = 1 ;
analysisSettings.posVariableLoadBC  = 2 ;
analysisSettings.stopTolDeltau      = 1e-14 ;
analysisSettings.stopTolForces      = 1e-8 ;
analysisSettings.stopTolIts         = 30 ;
analysisSettings.ALdominantDOF      = [3*6+3 -1] ;

otherParams              = struct() ;
otherParams.problemName  = 'plastic_2dframe' ;
% otherParams.plots_format = 'vtk' ;

[ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

[matUs, loadFactorsMat, modelSolutions ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;

displacements = matUs(4*6+3,:) ; % node with vertical load applied
loadfactors = loadFactorsMat(:,2) ;

moments_hist = zeros(4,length(modelSolutions)) ;
for i =1:length(modelSolutions)
    aux = modelSolutions{i}.localInternalForces(5) ;
    moments_hist(:,i) = [ aux.Mz; aux.Mz2; aux.Mz3; aux.tM ] ;
end
Mn1_numericONSAS = moments_hist(1,:) ;
Mn2_numericONSAS = moments_hist(2,:) ;
Mn3_numericONSAS = moments_hist(3,:) ;
tMn_numericONSAS = moments_hist(4,:) ;

% Plots

lw = 2 ; ms = 1 ; plotfontsize = 14 ;

figure('Name','Darvall-Mendis Frame / Plasticity (load factors)','NumberTitle','off') ;
hold on, grid on

plot(abs(displacements), loadfactors, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#0072BD") ;

labx = xlabel('Vertical displacements in free node (m, rad)') ;
laby = ylabel('\lambdaF') ;

legend('ONSAS \lambdaF [y]', 'location', 'Southeast') ;

set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Darvall-Mendis Frame / Plasticity (load factors)') ;

figure('Name','Darvall-Mendis Frame / Plasticity (Hinge Moment)','NumberTitle','off') ;
hold on, grid on

plot(abs(displacements), tMn_numericONSAS, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#D95319") ;

labx = xlabel('Generalized displacements in free node (m, rad)') ;
laby = ylabel('Hinge Moment') ;

legend('ONSAS tMn [y]', 'location', 'Southeast') ;

set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Darvall-Mendis Frame / Plasticity (Hinge Moment)') ;

print('-f1','../../../Tesis/tex/imagenes/DarvallMendisFrameLoadFactors.pdf','-dpdf') ;
print('-f2','../../../Tesis/tex/imagenes/DarvallMendisFrameHingeMoment.pdf','-dpdf') ;