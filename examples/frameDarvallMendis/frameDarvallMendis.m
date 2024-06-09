% =========================================================================

% Darvall-Mendis Frame
% Elements with embeded discontinuity

% Elastic-Plastic-Softening Analysis of Plane Frames
% Peter LePoer Darvall and Priyantha Anumddha Mendis
% Dept. of Civ. Engrg., Monash Univ., Victoria, Australia
% Journal of Structural Engineering, Vol. Ill, No. 4, April, 1985.

% =========================================================================

close all ; clear ;
addpath( genpath( [ pwd '/../..src'] ) ) ;

% assumed XY plane

% -------------------------------------------
% scalar parameters
% material
EI  = 77650 ;       % KN.m^2
kh1 = 29400 ;       % KN.m^2
kh2 = 273 ;
Ks  = -18000 ;      % KN.m

nu = 0.3 ;          % Poisson's ratio

% geometry
l  = 2.5 ;              % m
ty = 0.3 ;              % width cross section
tz = 0.4 ;              % height cross section
Inertia = tz*ty^3/12 ;  % m^4

E = EI/Inertia ;        % KN/m^2 [KPa]

A  = ty*tz ;            % m^2
Mc = 37.9 ;             % KN.m
My = 268 ;
Mu = 374 ;

% /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\
% ONSAS (NUMBER OF ELEMENTS 10)

% at the beginning..., there was no softening hinge
% soft_hinge_boolean = false ;

% number of finite elements
num_elem = 10 ;

materials             = struct() ;
materials.modelName   = 'plastic-2Dframe' ;
materials.modelParams = [ E Mc My Mu kh1 kh2 Ks nu ] ;

elements             = struct() ;
elements(1).elemType = 'node' ;

elements(2).elemType = 'frame' ;
elements(2).elemCrossSecParams = {'generic' ; [A 1 Inertia Inertia] } ;

boundaryConds                  = {} ;
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;

boundaryConds(2).imposDispDofs = [ 2 4 5] ;
boundaryConds(2).imposDispVals = [ 0 0 0 ] ;
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsBaseVals = [ 0 0 -1 0 0 0 ] ;
boundaryConds(2).loadsTimeFact = @(t) t ;

boundaryConds(3).imposDispDofs = [ 2 4 5] ;
boundaryConds(3).imposDispVals = [ 0 0 0 ] ;

% The coordinates of the nodes of the mesh are given by the matrix:
mesh = {} ;
xs = linspace(0,l,num_elem+1);
mesh.nodesCoords = [ xs' zeros(num_elem + 1, 2) ] ;

mesh.conecCell = {} ;

mesh.conecCell{ 1, 1 } = [ 0 1 1 1 ] ; % node

if num_elem>1

    for k=2:num_elem
        mesh.conecCell{ end+1, 1 } = [ 0 1 3 k ] ;
    end

end

for k=1:num_elem

    mesh.conecCell{ end+1, 1 } = [ 1 2 0 k k+1 ] ;

end

mesh.conecCell{ end+1, 1 } = [ 0 1 2 num_elem+1 ] ; % loaded node

initialConds = {} ;

analysisSettings                    = {} ;
analysisSettings.methodName         = 'arcLength' ;
analysisSettings.deltaT             = 1 ;
analysisSettings.incremArcLen       = [1e-4*ones(1,10) 1e-5*ones(1,3000) 1e-4*ones(1,1000)] ;
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

girosUltimoNodo_10 = matUs((num_elem+1)*6,:) ;
descensosUltimoNodo_10 = matUs((num_elem+1)*6-3,:) ;
factorescarga_10 = loadFactorsMat(:,2) ;

%{

% GRAPHICS

lw = 2 ; ms = 1 ; plotfontsize = 14 ;

figure('Name','Cantilever Beam / Plasticity (load factors)','NumberTitle','off') ;
hold on, grid on

plot(abs(girosUltimoNodo), factorescarga, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#EDB120") ;
plot(abs(descensosUltimoNodo), factorescarga, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#0072BD") ;

plot(abs(girosUltimoNodo_2), factorescarga_2, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#D95319") ;
plot(abs(descensosUltimoNodo_2), factorescarga_2, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#4DBEEE") ;

plot(abs(girosUltimoNodo_5), factorescarga_5, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#7E2F8E") ;
plot(abs(descensosUltimoNodo_5), factorescarga_5, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#77AC30") ;

plot(abs(girosUltimoNodo_10), factorescarga_10, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#0072BD") ;
plot(abs(descensosUltimoNodo_10), factorescarga_10, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#D95319") ;


labx = xlabel('Generalized displacements in free node (m, rad)') ;
laby = ylabel('Forces') ;

legend('ONSAS (1 elem) [\theta]', 'ONSAS (1 elem) [y]', 'ONSAS (2 elem) [\theta]', 'ONSAS (2 elem) [y]', 'ONSAS (5 elem) [\theta]', 'ONSAS (5 elem) [y]', 'ONSAS (10 elem) [\theta]', 'ONSAS (10 elem) [y]', 'location', 'Southeast') ;

set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Cantilever Beam / Plasticity (load factors)') ;

figure('Name','Cantilever Beam / Plasticity (validation)','NumberTitle','off') ;
hold on, grid on

plot(abs(girosUltimoNodo), abs(Mn1_validation), '-x' , 'linewidth', lw, 'markersize', ms, "Color", "#EDB120") ;
plot(abs(descensosUltimoNodo), abs(Mn1_validation), '-x' , 'linewidth', lw, 'markersize', ms, "Color", "#0072BD") ;

plot(abs(matdes(6,1:length(load_factors)-1)), Mn,'-x' , 'linewidth', lw, 'markersize', ms, "Color", "#D95319") ;
plot(abs(matdes(4,1:length(load_factors)-1)), Mn, '-x' , 'linewidth', lw, 'markersize', ms, "Color", "#77AC30") ;

labx = xlabel('Generalized displacements in free node (m, rad)') ;
laby = ylabel('Bulk Moment at the first integration point (KN.m)') ;
legend('Semi Analytic (1 elem) [\theta]', 'Semi Analytic (1 elem) [y]', 'ALGOL (1 elem) [\theta]', 'ALGOL (1 elem) [y]', 'location', 'Southeast') ;

set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Cantilever Beam / Plasticity (validation)') ;

print('-f1','../../../Tesis/tex/imagenes/Load_factors.pdf','-dpdf') ;
print('-f2','../../../Tesis/tex/imagenes/Validation.pdf','-dpdf') ;

%}