% Copyright 2024, ONSAS Authors (see documentation)
%
% This file is part of ONSAS.
%
% ONSAS is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% ONSAS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.

% =========================================================================

% Euler-Bernoulli element with embeded discontinuity

% Numerical modeling of softening hinges in thin Euler–Bernoulli beams
% Francisco Armero, David Ehrlich / University of California, Berkeley

% Embedded discontinuity finite element formulation
% For failure analysis of planar reinforced concrete beams and frames
% Miha Jukić, Boštjan Brank / University of Ljubljana
% Adnan Ibrahimbegović / Ecole normale supérieure de Cachan

% =========================================================================

% numerical example
% cantilever beam loaded with a vertical force at the free end

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

[matUs, loadFactorsMat, modelSolutions ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;

girosUltimoNodo_10 = matUs((num_elem+1)*6,:) ;
descensosUltimoNodo_10 = matUs((num_elem+1)*6-3,:) ;
factorescarga_10 = loadFactorsMat(:,2) ;

% /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\
% ONSAS (NUMBER OF ELEMENTS 1)

% at the beginning..., there was no softening hinge
% soft_hinge_boolean = false ;

% number of finite elements
num_elem = 1 ;

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
analysisSettings.incremArcLen       = 1e-3*ones(1,1000) ;
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

[matUs, loadFactorsMat, modelSolutions ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;

girosUltimoNodo = matUs((num_elem+1)*6,:) ;
descensosUltimoNodo = matUs((num_elem+1)*6-3,:) ;
factorescarga = loadFactorsMat(:,2) ;

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

Mn1_validation = zeros(1,length(matUs(1,:))) ;

for i = 1:length(matUs(1,:))

    v1 = matUs(3,i) ;
    v2 = matUs(9,i) ;

    theta1 = matUs(6,i) ;
    theta2 = matUs(12,i) ;

    kappa_plas_n{i} = kappa_plas_n1 ;
    xin1val(i) = xin11val(1) ;
    
    [kappa_plas_n1, xin11val, xin21val, alfan1, Mn1] = softHinge1DOF_semiAnalyticSol(v1, v2, theta1, theta2 , xd, alpha, xin1val(i), xin2val(i), kappa_plas_n{i}, Mc, My, Mu, kh1, kh2, Ks, E, Inertia, l) ;

    Mn1_validation(i) = Mn1(1) ;

end

% /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\
% Algorithm without ONSAS (NUMBER OF ELEMENTS 1)

addpath( genpath( [ pwd '/../../src'] ) ) ;

% Mc, My, Mu / from the moment-curvature diagram
% kh1, kh2   / hardening modules
% Ks         / from the moment-rotation jump diagram

l = 2.5 ;           % m
A = 0.4*0.3 ;       % m^2
E = 30000000 ;      % KN/m^2 KPa
EI = 77650 ;        % KN.m^2
Iy = EI/E ;         % m^4
Mc = 37.9 ;         % KN.m
My = 268 ;
Mu = 374 ;
kh1 = 29400 ;       % KN.m^2
kh2 = 272 ;
Ks = -18000 ;       % KN.m

freedofs = [2 4 6]; % u2 v2 theta2

% at the beginning..., there was no softening hinge
soft_hinge_boolean = false ;

% Gauss-Lobatto Quadrature with 3 integration points [a (a+b)/2 b]
npi = 3 ;
xpi = [0 l/2 l] ;
wpi = [1/3 4/3 1/3] * l * 0.5 ;

nu   = 0.3 ;
tol1 = 1e-8;
tol2 = 1e-8 ;
tolk = 15 ;

% initial values
dn   = [0 0 0 0 0 0]' ;
Fint = [0 0 0 0 0 0]' ;
tM   = 0 ;

kpn  = zeros(npi,1) ;
xin1 = zeros(npi,1) ;
xin2 = 0 ;

kpn1  = zeros(npi,1) ;
xin11 = zeros(npi,1) ;
xin21 = 0 ;

khat1 = zeros(npi,1) ;

M1 = zeros(npi,1) ;
Fn = zeros(npi,1) ;

alfan = 0 ;

xd = 0 ;
xdi = 1 ;

Final_force = 183 ; % value of the final force

load_case = [0 0 0 1 0 0]' ; % load applied in vertical direction (Y)
load_factors = 0:Final_force ;

% --- element params ---
elemParams = [l A Iy] ;

% --- elastoplastic params ---
elastoplasticParams = [E Mc My Mu kh1 kh2 Ks] ;

matdes = zeros (6, Final_force+1) ;

matdes(:,1) = dn ;

gxin = zeros(Final_force, 1) ;
gxin2 = zeros(Final_force, 1) ;
gkpn = zeros(Final_force, 1) ;

Mn = zeros(Final_force, 1) ;
TM = zeros(Final_force, 1) ;

Alf = zeros(Final_force, 1) ;

for ind = 2:length(load_factors)

    curr_load_factor = load_factors(ind) ;

    Fext = load_case * curr_load_factor ;

    dnk = matdes(:,ind-1) ;

    % iteration vars
    converged_boolean = false ;
    
    % set iterations zero
    k = 0 ;

    gxin(ind-1,1)   = xin1(1)   ;
    gxin2(ind-1,1)  = xin2      ;
    gkpn(ind-1,1)   = kpn(1)    ;
    Mn(ind-1,1)     = M1(1)     ;
    TM(ind-1,1)     = tM        ;
    Fn(ind-1,1)     = Fint(4)   ;
    Alf(ind-1,1)    = alfan     ;

    while converged_boolean == false && k < tolk

        k = k + 1 ;

        [soft_hinge_boolean, Fint, M1, Kelement, kpn1, xin11, xin21, alfan1, xd, xdi, tM] = framePlastic(soft_hinge_boolean, dnk, kpn, xin1, xin2, alfan, xd, xdi, tM, elemParams, elastoplasticParams) ;

        residualForce = Fext - Fint ;

        Krelement = Kelement(freedofs,freedofs) ;

        residualForceRed = residualForce(freedofs) ;

        % system of equilibrium equations

        deltadred = Krelement\residualForceRed ;

        % /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\

        deltad = zeros(6,1) ;
        deltad(freedofs) = deltadred ;

        dnk1 = dnk + deltad ;

        dnk     = dnk1   ;

        kpn     = kpn1   ;
        xin1    = xin11  ;

        xin2    = xin21  ;
        alfan   = alfan1 ;

        norm1 = norm(deltadred) ;
        norm2 = norm(residualForceRed) ;

        converged_boolean = norm1 < tol1 || norm2 < tol2 ;

    end

    matdes(:,ind) = dnk1 ;

end

% GRAPHICS

lw = 2 ; ms = 1 ; plotfontsize = 14 ;

figure('Name','Cantilever Beam / Plasticity (load factors)','NumberTitle','off') ;
hold on, grid on

plot(abs(girosUltimoNodo), factorescarga, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#EDB120") ;
plot(abs(descensosUltimoNodo), factorescarga, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#0072BD") ;

% plot(abs(girosUltimoNodo_2), factorescarga_2, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#0072BD") ;
% plot(abs(descensosUltimoNodo_2), factorescarga_2, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#77AC30") ;

% plot(abs(girosUltimoNodo_5), factorescarga_5, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#D95319") ;
% plot(abs(descensosUltimoNodo_5), factorescarga_5, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#77AC30") ;

plot(abs(girosUltimoNodo_10), factorescarga_10, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#D95319") ;
plot(abs(descensosUltimoNodo_10), factorescarga_10, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#0072BD") ;


labx = xlabel('Generalized displacements in free node (m, rad)') ;
laby = ylabel('Forces') ;

legend('ONSAS (1 elem) [\theta]', 'ONSAS (1 elem) [y]','ONSAS (10 elem) [\theta]', 'ONSAS (10 elem) [y]', 'location', 'Southeast') ;

set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Cantilever Beam / Plasticity (load factors)') ;

%{

figure('Name','Cantilever Beam / Plasticity (validation)','NumberTitle','off') ;
hold on, grid on

% plot(abs(girosUltimoNodo), abs(Mn1_validation), '-x' , 'linewidth', lw, 'markersize', ms, "Color", "#EDB120") ;
% plot(abs(descensosUltimoNodo), abs(Mn1_validation), '-x' , 'linewidth', lw, 'markersize', ms, "Color", "#0072BD") ;

plot(abs(matdes(6,1:length(load_factors)-1)), Mn,'-x' , 'linewidth', lw, 'markersize', ms, "Color", "#D95319") ;
plot(abs(matdes(4,1:length(load_factors)-1)), Mn, '-x' , 'linewidth', lw, 'markersize', ms, "Color", "#77AC30") ;

labx = xlabel('Generalized displacements in free node (m, rad)') ;
laby = ylabel('Bulk Moment at the first integration point (KN.m)') ;
legend('ALGOL (1 elem) [\theta]', 'ALGOL (1 elem) [y]', 'location', 'Southeast') ;

% 'Semi Analytic (1 elem) [\theta]', 'Semi Analytic (1 elem) [y]',

set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Cantilever Beam / Plasticity (validation)') ;

%}

print('-f1','../../../Tesis/tex/imagenes/Load_factors.pdf','-dpdf') ;
% print('-f2','../../../Tesis/tex/imagenes/Validation.pdf','-dpdf') ;