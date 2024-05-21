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
Mu = 376 ;

% ----------------------------------------------------------------------------------
% /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\
% ONSAS solution with 1 element

% at the beginning..., there was no softening hinge
soft_hinge_boolean = false ;

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

% The coordinates of the nodes of the mesh are given by the matrix:
mesh = struct() ;
mesh.nodesCoords = [ 0 0 0 ; l 0 0 ] ;

mesh.conecCell = {} ;
mesh.conecCell{ 1,     1 } = [ 0 1 1  1   ] ; % node
mesh.conecCell{ end+1, 1 } = [ 0 1 2  2   ] ; % node
mesh.conecCell{ end+1, 1 } = [ 1 2 0  1 2 ] ; % frame

initialConds = struct() ;

analysisSettings                    = struct() ;
analysisSettings.methodName         = 'arcLength' ;
analysisSettings.deltaT             = 1 ;
% analysisSettings.incremArcLen       = [1e-3*ones(1,460) 1e-4*ones(1,1000) 1e-5*ones(1,500)] ;

% analysisSettings.incremArcLen       = [1e-3*ones(1,30)  ] ;
% analysisSettings.incremArcLen       = [.2e-3*ones(1,2)   ] ;
analysisSettings.incremArcLen       = [.2e-3*ones(1,20) 1e-3*ones(1,22)  ] ;
% analysisSettings.incremArcLen       = [1e-3*ones(1,428)  ] ;
% analysisSettings.incremArcLen       = [1e-3*ones(1,425) 1e-4*ones(1,40) ] ;

analysisSettings.finalTime          = length(analysisSettings.incremArcLen) ;
analysisSettings.iniDeltaLamb       = 1 ;
analysisSettings.posVariableLoadBC  = 2 ;
analysisSettings.stopTolDeltau      = 1e-8 ;
analysisSettings.stopTolForces      = 1e-8 ;
analysisSettings.stopTolIts         = 15 ;
analysisSettings.ALdominantDOF      = [ 2*6-3  -1] ;

otherParams              = struct() ;
otherParams.problemName  = 'plastic_2dframe' ;
% otherParams.plots_format = 'vtk' ;

[ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

[matUs, loadFactorsMat, cellFint, cellStress ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;

girosUltimoNodo     = matUs((1+1)*6,:) ;
descensosUltimoNodo = matUs((1+1)*6-3,:) ;
factorescargaONSAS  = loadFactorsMat(:,2) ;

moments_hist = zeros(4,length(cellFint)) ;
for i =1:length(cellFint)
    aux = cellFint{1,i} ;
    moments_hist(:,i) = aux(1:4) ;
end
Mn1_numericONSAS = moments_hist(1,:) ;
Mn2_numericONSAS = moments_hist(2,:) ;
Mn3_numericONSAS = moments_hist(3,:) ;
tMn_numericONSAS = moments_hist(4,:) ;

% ----------------------------------------------------------------------------------
% /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\
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

    [kappa_plas_n1, xin11val, xin21val, alfan1, Mn1] = softHinge1DOF_semiAnaSol(v1, v2, theta1, theta2 , xd, alpha, xin1val(i), xin2val(i), kappa_plas_n{i}, Mc, My, Mu, kh1, kh2, Ks, E, Inertia, l) ;

    Mn1_semianalytic(i) = Mn1(1) ;

end

% ----------------------------------------------------------------------------------
% numerical solution B

[ Mn_numer, Fn, matdes ] = softHinge1DOF_numericSol(l, A, E, Inertia, Mc, My, Mu, kh1, kh2, Ks) ;

desp_numer = matdes(4,:) ;


% ----------------------------------------------------------------------------------
% Plots

lw = 1.4 ; ms = 1 ; plotfontsize = 14 ;

% figure('Name','Cantilever Beam / Plasticity (load factors)','NumberTitle','off') ;
% hold on, grid on

% plot(abs(girosUltimoNodo),     factorescarga, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#EDB120") ;
% plot(abs(descensosUltimoNodo), factorescarga, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#0072BD") ;

% labx = xlabel('Generalized displacements in free node (m, rad)') ;
% laby = ylabel('Forces') ;
% legend('ONSAS (1 elem) [\theta]', 'ONSAS (1 elem) [y]', 'location', 'Southeast') ;
% set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
% set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
% title('Cantilever Beam / Plasticity (load factors)') ;

figure('Name','Cantilever Beam / Plasticity (validation)','NumberTitle','off') ;
hold on, grid on

paso = 2 ;

% plot(-girosUltimoNodo(1:paso:length(Mn1_semianalytic)), -Mn1_semianalytic(1:paso:length(Mn1_semianalytic)), '^', 'linewidth', lw, 'markersize', ms*10, "Color", "#77AC30") ;
plot(-descensosUltimoNodo(1:paso:length(Mn1_semianalytic)), -Mn1_semianalytic(1:paso:length(Mn1_semianalytic)), '^', 'linewidth', lw, 'markersize', ms*10, "Color", "#EDB120") ;

% plot(-girosUltimoNodo, -Mn1_numericONSAS, '-x' , 'linewidth', lw, 'markersize', ms*5, "Color", "#EDB120") ;
plot(-descensosUltimoNodo, -Mn1_numericONSAS, '-x' , 'linewidth', lw, 'markersize', ms*5, "Color", "#0072BD") ;

plot(-descensosUltimoNodo, -Mn2_numericONSAS, 'k-s' , 'linewidth', lw, 'markersize', ms*5) ;
plot(-descensosUltimoNodo, -Mn3_numericONSAS, 'r-o' , 'linewidth', lw, 'markersize', ms*5) ;

plot(-descensosUltimoNodo, -tMn_numericONSAS, 'g-^' , 'linewidth', lw, 'markersize', ms*10) ;

plot(-descensosUltimoNodo, factorescargaONSAS*l, 'c-v' , 'linewidth', lw, 'markersize', ms*10) ;

final = round(length(Mn_numer)/8) ;
plot(desp_numer(1:final), Mn_numer(1:final), 'r-+' , 'linewidth', lw, 'markersize', ms*10) ;

%{

plot(abs(girosUltimoNodo), abs(moments_hist(2,:)), '-x' , 'linewidth', lw, 'markersize', ms, "color", "#A2142F") ;
plot(abs(descensosUltimoNodo), abs(moments_hist(2,:)), '-x' , 'linewidth', lw, 'markersize', ms, "Color", "#7E2F8E") ;

plot(abs(girosUltimoNodo), abs(moments_hist(3,:)), '-x' , 'linewidth', lw, 'markersize', ms, "Color", "#77AC30") ;
plot(abs(descensosUltimoNodo), abs(moments_hist(3,:)), '-x' , 'linewidth', lw, 'markersize', ms, "Color", "#4DBEEE") ;

plot(abs(girosUltimoNodo), abs(moments_hist(4,:)), '-x' , 'linewidth', lw, 'markersize', ms, "Color", "#D95319") ;
plot(abs(descensosUltimoNodo), abs(moments_hist(4,:)), '-x' , 'linewidth', lw, 'markersize', ms, "Color", "#7E2F8E") ;

%}

labx = xlabel('Generalized displacements in free node (m, rad)') ;
laby = ylabel('Bulk Moments at the integration points (KN.m)') ;

% legend('Semi-Analytic Mp1 [\theta]', 'Semi-Analytic Mp1 [y]', 'ONSAS Mp1 [\theta]', 'ONSAS Mp1 [y]', 'location', 'Southeast') ;
legend('Semi-Analytic Mp1 [y]', 'ONSAS Mp1 [y]', 'ONSAS Mp2 [y]', 'ONSAS Mp3 [y]', 'tMn', 'factorcarga onsas', 'location', 'Southeast') ;

%  'ONSAS Mp2 [\theta]', 'ONSAS Mp2 [y]', 'ONSAS Mp3 [\theta]', 'ONSAS Mp3 [y]'

set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Cantilever Beam / Plasticity (validation)') ;