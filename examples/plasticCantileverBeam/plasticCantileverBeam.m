
% numerical example
% cantilever beam loaded with a vertical force at the free end
% 
close all, if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear all, end
addpath( genpath( [ pwd '/../../src'] ) );

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
Mu = 376 ;

% /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\
% ONSAS (NUMBER OF ELEMENTS 1)

% at the beginning..., there was no softening hinge
soft_hinge_boolean = false ;

% number of finite elements
num_elem = 1 ;

global historic_parameters
historic_parameters = [] ;

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
mesh.nodesCoords = [  xs' zeros(num_elem + 1, 2) ] ;

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
analysisSettings.incremArcLen       = [1e-3*ones(1,460) 1e-4*ones(1,1000) 1e-5*ones(1,500)] ;
analysisSettings.finalTime          = length(analysisSettings.incremArcLen) ;
analysisSettings.iniDeltaLamb       = 1 ;
analysisSettings.posVariableLoadBC  = 2 ;
analysisSettings.stopTolDeltau      = 1e-8 ;
analysisSettings.stopTolForces      = 1e-8 ;
analysisSettings.stopTolIts         = 15 ;
analysisSettings.ALdominantDOF      = [ (num_elem+1)*6-3  -1] ;

otherParams              = struct() ;
otherParams.problemName  = 'plastic_2dframe' ;
# otherParams.plots_format = 'vtk' ;

[ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

[matUs, loadFactorsMat, cellFint, cellStress ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;

girosUltimoNodo = matUs((num_elem+1)*6,:) ;
descensosUltimoNodo = matUs((num_elem+1)*6-3,:) ;
factorescarga = loadFactorsMat(:,2) ;

moments_hist = zeros(4,length(cellFint));
for i =1:length(cellFint)
    aux = cellFint{1,i};
    moments_hist(:,i) = aux(1:4) ;
end

Mn1_numeric = moments_hist(1,:);

% /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\
% semi-analitic validation / ONSAS with the function moments_plus_internal_variables

xd = 0 ;
alpha = 0 ;
xin1val = zeros(1,length(matUs(1,:))) ;
kappa_plas_n = zeros(1,length(matUs(1,:))) ;
kappa_plas_n1 = zeros(1,3) ;
xin11val = zeros(1,3) ;

Mn1_semianalytic = zeros(1,length(matUs(1,:))) ;

for i = 1:length(matUs(1,:))

    v1 = matUs(3,i) ;
    v2 = matUs(9,i) ;

    theta1 = matUs(6,i) ;
    theta2 = matUs(12,i) ;

    kappa_plas_n(i) = kappa_plas_n1(1) ;
    xin1val(i) = xin11val(1) ;

    [kappa_plas_n1, xin11val, Mn1] = moments_plus_internal_variables(v1, v2, theta1, theta2 , xd, alpha, xin1val(i), kappa_plas_n(i), Mc, My, kh1, kh2, E, Inertia, l) ;

    Mn1_semianalytic(i) = Mn1(1) ;

end

% Plots

lw = 2 ; ms = 1 ; plotfontsize = 14 ;

figure('Name','Cantilever Beam / Plasticity (load factors)','NumberTitle','off') ;
hold on, grid on

plot(abs(girosUltimoNodo), factorescarga, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#EDB120") ;
plot(abs(descensosUltimoNodo), factorescarga, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#0072BD") ;

labx = xlabel('Generalized displacements in free node (m, rad)') ;
laby = ylabel('Forces') ;
legend('ONSAS (1 elem) [\theta]', 'ONSAS (1 elem) [y]', 'location', 'Southeast') ;
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Cantilever Beam / Plasticity (load factors)') ;


figure('Name','Cantilever Beam / Plasticity (validation)','NumberTitle','off') ;
hold on, grid on

plot(abs(girosUltimoNodo), abs(Mn1_semianalytic), '-x' , 'linewidth', lw, 'markersize', ms, "Color", "#EDB120") ;
plot(abs(descensosUltimoNodo), abs(Mn1_semianalytic), '-x' , 'linewidth', lw, 'markersize', ms, "Color", "#0072BD") ;

plot(abs(girosUltimoNodo), abs(moments_hist(1,:)), 'r-x' , 'linewidth', lw, 'markersize', ms ) ;
plot(abs(descensosUltimoNodo), abs(moments_hist(1,:)), 'r-x' , 'linewidth', lw, 'markersize', ms) ;

plot(abs(girosUltimoNodo), abs(moments_hist(2,:)), 'g-x' , 'linewidth', lw, 'markersize', ms ) ;
plot(abs(descensosUltimoNodo), abs(moments_hist(2,:)), 'g-x' , 'linewidth', lw, 'markersize', ms) ;

plot(abs(girosUltimoNodo), abs(moments_hist(3,:)), 'c-x' , 'linewidth', lw, 'markersize', ms ) ;
plot(abs(descensosUltimoNodo), abs(moments_hist(3,:)), 'c-x' , 'linewidth', lw, 'markersize', ms) ;

plot(abs(girosUltimoNodo), abs(moments_hist(4,:)), 'k-x' , 'linewidth', lw, 'markersize', ms ) ;
plot(abs(descensosUltimoNodo), abs(moments_hist(4,:)), 'k-x' , 'linewidth', lw, 'markersize', ms) ;

labx = xlabel('Generalized displacements in free node (m, rad)') ;
laby = ylabel('Bulk Moment at the first integration point (KN.m)') ;
legend('Semi-Analytic [\theta]', 'Semi-Analytic [y]', 'ONSAS [\theta]', 'ONSAS [y]', 'location', 'Southeast') ;
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Cantilever Beam / Plasticity (validation)') ;
