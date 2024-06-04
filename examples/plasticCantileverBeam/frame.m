% ----------------------------------------------------------------------------------
% numerical example
% cantilever beam loaded with a vertical force at the free end
% assumed XY plane
% ----------------------------------------------------------------------------------
%
close all,
if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear all, end
addpath( genpath( [ pwd '/../../src'] ) );
%
% ----------------------------------------------------------------------------------
% problem parameters
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


nelem = 10; nnodes = nelem+1;
% ----------------------------------------------------------------------------------
% /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\
% ONSAS solution with 1 element

% material
materials             = struct() ;
materials.modelName   = 'plastic-2Dframe' ;
materials.modelParams = [ E Mc My Mu kh1 kh2 Ks nu ] ;

% material
elements             = struct() ;
elements(1).elemType = 'node' ;
%
elements(2).elemType = 'frame' ;
elements(2).elemCrossSecParams = {'generic' ; [A 1 Inertia Inertia] } ;

boundaryConds                  = struct() ;
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
mesh = struct() ;
mesh.nodesCoords = [ linspace(0,l,nnodes)' zeros(nnodes,2)] ;

mesh.conecCell = {} ;
mesh.conecCell{ 1,     1 } = [ 0 1 1  1   ] ; % node
mesh.conecCell{ end+1, 1 } = [ 0 1 2  nnodes   ] ; % loaded node
for i=2:nnodes-1
    mesh.conecCell{ end+1, 1 } = [ 0 1 3  i   ] ; % intermediate node
end

for i=1:nelem
    mesh.conecCell{ end+1, 1 } = [ 1 2 0  i i+1 ] ; % 
end


initialConds = struct() ;

analysisSettings                    = struct() ;
analysisSettings.methodName         = 'arcLength' ;
analysisSettings.deltaT             = 1 ;

% the softening hinge is activated  1 elem
# analysisSettings.incremArcLen       = [1e-3*ones(1,832) eps*ones(1,1) 1e-3*ones(1,100)] ;

# analysisSettings.incremArcLen       = [1e-3*ones(1,432) ] ;

% second stage of plastic hardening
analysisSettings.incremArcLen       = [.2e-3*ones(1,320) .8e-4*ones(1,1000) .1e-4*ones(1,500)  ] ;

% first stage of plastic hardening
#  analysisSettings.incremArcLen       = [.3e-3*ones(1,2) ] ;

analysisSettings.finalTime          = length(analysisSettings.incremArcLen) ;
analysisSettings.iniDeltaLamb       = 1 ;
analysisSettings.posVariableLoadBC  = 2 ;
analysisSettings.stopTolDeltau      = 1e-14 ;
analysisSettings.stopTolForces      = 1e-8 ;
analysisSettings.stopTolIts         = 30 ;
analysisSettings.ALdominantDOF      = [nnodes*6-3 -1] ;

otherParams              = struct() ;
otherParams.problemName  = 'plastic_2dframe' ;
# otherParams.plots_format = 'vtk' ;

[ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

[matUs, loadFactorsMat, modelSolutions ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;

girosUltimoNodo     = matUs((nnodes)*6,:) ;
descensosUltimoNodo = matUs((nnodes)*6-3,:) ;
factorescargaONSAS  = loadFactorsMat(:,2) ;

moments_hist = zeros(4,length(modelSolutions)) ;

moments_integrados_izq = zeros(nelem,length(modelSolutions)) ;
moments_integrados_der = zeros(nelem,length(modelSolutions)) ;

for i =1:length(modelSolutions)
    aux = modelSolutions{i}.localInternalForces(1) ;
    moments_hist(:,i) = [ aux.Mz; aux.Mz2; aux.Mz3; aux.tM ] ;
end
Mn1_numericONSAS = moments_hist(1,:) ;
Mn2_numericONSAS = moments_hist(2,:) ;
Mn3_numericONSAS = moments_hist(3,:) ;
tMn_numericONSAS = moments_hist(4,:) ;

for i =1:length(modelSolutions)
  for j=1:nelem
    aux = modelSolutions{i}.localInternalForces(j) ;
    moments_integrados_izq(j,i) = aux.Mz_integrado_izq ;
    moments_integrados_der(j,i) = aux.Mz_integrado_der ;
  end
end

stop
% ----------------------------------------------------------------------------------
% /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\
% semi-analitic validation / ONSAS with the function moments_plus_internal_variables

xd = 0 ;
alpha = 0 ;
xin1val = zeros(1,length(matUs(1,:))) ;
xin2val = zeros(1,length(matUs(1,:))) ;
kappa_plas_n = cell(1,length(matUs(1,:))) ;
kappa_plas_n1 = zeros(1,3) ;
xin11val = zeros(1,3) ;
xin21val = zeros(1,3) ;

Mn1_semianalytic = zeros(1,length(matUs(1,:))) ;

for i = 1:length(matUs(1,:))

    v1 = matUs(3,i) ;
    v2 = matUs(9,i) ;

    theta1 = matUs(6,i) ;
    theta2 = matUs(12,i) ;

    kappa_plas_n{i} = kappa_plas_n1 ;
    xin1val(i) = xin11val(1) ;
    
    [kappa_plas_n1, xin11val, xin21val, alfan1, Mn1] = softHinge1DOF_semiAnalyticSol(v1, v2, theta1, theta2 , xd, alpha, xin1val(i), xin2val(i), kappa_plas_n{i}, Mc, My, Mu, kh1, kh2, Ks, E, Inertia, l) ;

    Mn1_semianalytic(i) = Mn1(1) ;

end

% ----------------------------------------------------------------------------------
% numerical solution ALGOL

[Mn_numer, Fn, matdes] = softHinge1DOF_numericSol(l, A, E, Inertia, Mc, My, Mu, kh1, kh2, Ks) ;

desp_numer = matdes(4,:) ;

% ----------------------------------------------------------------------------------
% Plots

lw = 1.6 ; ms = 1 ; plotfontsize = 14 ;

figure('Name','Cantilever Beam / Plasticity (validation)','NumberTitle','off') ;
hold on, grid on

step = 1 ;

plot(-descensosUltimoNodo, -Mn1_numericONSAS, '-x', 'linewidth', lw, 'markersize', ms*4, "Color", "#0072BD") ;


plot(desp_numer(1:length(Mn_numer)), Mn_numer(1:length(Mn_numer)), '-x' , 'linewidth', lw, 'markersize', ms, "Color", "#A2142F") ;

plot(-descensosUltimoNodo, -Mn2_numericONSAS, '-x', 'linewidth', lw, 'markersize', ms*4, "Color", "#0072BD") ;

% plot(-descensosUltimoNodo(1:step:length(Mn1_semianalytic)), -Mn1_semianalytic(1:step:length(Mn1_semianalytic)), 'b^', 'linewidth', lw, 'markersize', ms*4, "Color", "#EDB120") ;

labx = xlabel('Generalized displacements in free node (m, rad)') ;
laby = ylabel('Bulk Moments at the integration points (KN.m)') ;

legend('ONSAS Mp1 [y]', 'ALGOL Mp1 [y]','onsas mp2', 'location', 'Southeast') ;

% 'Semi-Analytic Mp1 [y]', 

set(gca, 'linewidth', 1, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Cantilever Beam / Plasticity (validation)') ;

figure
plot(-descensosUltimoNodo , factorescargaONSAS )
