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
kh1 = 1         ;   % KN.m^2
kh2 = 1         ;

% Ks    = -18000        ;   % KN.m
% a     = -0.04         ;
% Ks    = a*EI/10/l     ;   % KN.m

Ks  = -1e-6             ;   % almost zero

nu  = 0.3               ;   % Poisson's ratio

Mu_columns  = 158.18    ;   % KN.m
Mu_beams    = 169.48    ;   % KN.m

My = 1e18               ; % almost infinite
Mc = 1e18               ; % almost infinite

% /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\

materials             = struct() ;
materials(1).modelName   = 'plastic-2Dframe' ;
materials(1).modelParams = [ E Mc My Mu_beams kh1 kh2 Ks nu ] ;

materials(2).modelName   = 'plastic-2Dframe' ;
materials(2).modelParams = [ E Mc My Mu_columns kh1 kh2 Ks nu ] ;

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
mesh.conecCell{ end+1, 1 } = [ 0 1 1 9 ] ; % node 9 fixed end support

mesh.conecCell{ end+1, 1 } = [ 0 1 2 5 ] ; % node 5 with vertical load applied

mesh.conecCell{ end+1, 1 } = [ 0 1 3 2 ] ; % node 2
mesh.conecCell{ end+1, 1 } = [ 0 1 3 3 ] ; % node 3
mesh.conecCell{ end+1, 1 } = [ 0 1 3 4 ] ; % node 4
mesh.conecCell{ end+1, 1 } = [ 0 1 3 6 ] ; % node 6
mesh.conecCell{ end+1, 1 } = [ 0 1 3 7 ] ; % node 7
mesh.conecCell{ end+1, 1 } = [ 0 1 3 8 ] ; % node 8

% frame elements
mesh.conecCell{ end+1, 1 } = [ 2 2 0   1 2] ;
mesh.conecCell{ end+1, 1 } = [ 2 2 0   2 3] ;
mesh.conecCell{ end+1, 1 } = [ 1 2 0   3 4] ;
mesh.conecCell{ end+1, 1 } = [ 1 2 0   4 5] ;
mesh.conecCell{ end+1, 1 } = [ 1 2 0   5 6] ;
mesh.conecCell{ end+1, 1 } = [ 1 2 0   6 7] ;
mesh.conecCell{ end+1, 1 } = [ 2 2 0   7 8] ;
mesh.conecCell{ end+1, 1 } = [ 2 2 0   8 9] ;

initialConds = {} ;

analysisSettings                    = {} ;
analysisSettings.methodName         = 'arcLength' ;
analysisSettings.deltaT             = 1 ;
analysisSettings.incremArcLen       = 1e-5*ones(1,1839) ;
analysisSettings.finalTime          = length(analysisSettings.incremArcLen) ;
analysisSettings.iniDeltaLamb       = 1 ;
analysisSettings.posVariableLoadBC  = 2 ;
analysisSettings.stopTolDeltau      = 1e-12 ;
analysisSettings.stopTolForces      = 1e-12 ;
analysisSettings.stopTolIts         = 10 ;
analysisSettings.ALdominantDOF      = [4*6+3 -1] ;

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

% Create arrow
annotation(figure(1),'arrow',[0.416071428571428 0.471428571428571],...
    [0.707142857142858 0.616666666666667]);

% Create textarrow
annotation(figure(1),'textarrow',[0.758928571428571 0.732142857142857],...
    [0.86904761904762 0.780952380952382]);

% Create textbox
annotation(figure(1),'textbox',...
    [0.470642857142857 0.571428571428575 0.157928571428572 0.0595238095238096],...
    'String','FIRST HINGE',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',12,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

% Create textbox
annotation(figure(1),'textbox',...
    [0.670642857142854 0.721428571428574 0.15614285714286 0.054761904761909],...
    'String',{'SECOND','HINGE'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',12,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

% Create textbox
annotation(figure(1),'textbox',...
    [0.817857142857138 0.726190476190483 0.130357142857147 0.0595238095238095],...
    'String',{'THIRD HINGE'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',12,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

% Create arrow
annotation(figure(1),'arrow',[0.873214285714286 0.841071428571429],...
    [0.888095238095238 0.780952380952381]);

figure('Name','Darvall-Mendis Frame / Plasticity (Hinge Moment)','NumberTitle','off') ;
hold on, grid on

plot(abs(displacements), abs(tMn_numericONSAS), '-x', 'linewidth', lw, 'markersize', ms, "Color", "#D95319") ;

labx = xlabel('Generalized displacements in free node (m, rad)') ;
laby = ylabel('Hinge Moment') ;

legend('ONSAS tMn [y]', 'location', 'Southeast') ;

set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Darvall-Mendis Frame / Plasticity (Hinge Moment)') ;

print('-f1','../../../Tesis/tex/imagenes/DarvallMendisFrameLoadFactors.pdf','-dpdf') ;
print('-f2','../../../Tesis/tex/imagenes/DarvallMendisFrameHingeMoment.pdf','-dpdf') ;