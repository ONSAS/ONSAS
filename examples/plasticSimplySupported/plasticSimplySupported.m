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

close all ; clear ;
addpath( genpath( [ pwd '/../../src'] ) ) ;

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
l1  = 1.0 ;              % m
l2  = 1.0 ;              % m
ty = 0.3 ;              % width cross section
tz = 0.4 ;              % height cross section
Inertia = tz*ty^3/12 ;  % m^4

E = EI/Inertia ;        % KN/m^2 [KPa]

A  = ty*tz ;            % m^2
Mc = 37.9 ;             % KN.m
My = 268 ;
Mu = 374 ;

materials             = struct() ;
materials.modelName   = 'plastic-2Dframe' ;
materials.modelParams = [ E Mc My Mu kh1 kh2 Ks nu ] ;

Pc = Mc * (l1+l2)/(l1*l2)

deltac = Pc * (l2*((l1+l2)^2-l2^2)^(1.5) )/(9 *sqrt(3)*EI)

elements             = struct() ;
elements(1).elemType = 'node' ;

elements(2).elemType = 'frame' ;
elements(2).elemCrossSecParams = {'generic' ; [A 1 Inertia Inertia] } ;

boundaryConds                  = {} ;
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 ] ;

boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsBaseVals = [ 0 0 -1 0 0 0 ] ;
boundaryConds(2).loadsTimeFact = @(t) t ;
boundaryConds(2).imposDispDofs = [ 2 4 5] ;
boundaryConds(2).imposDispVals = [ 0 0 0 ] ;

% The coordinates of the nodes of the mesh are given by the matrix:
mesh = {} ;
num_elem = 2 ;
xs = linspace(0,l1+l2,num_elem+1);
mesh.nodesCoords = [ 0 0 0 ;  l1 0 0; l1+l2 0 0 ] ;

mesh.conecCell = {} ;

mesh.conecCell{ 1    , 1 } = [ 0 1 1  1 ] ; % node with fixed end support
mesh.conecCell{ end+1, 1 } = [ 0 1 1  3 ] ; % node with fixed end support
mesh.conecCell{ end+1, 1 } = [ 0 1 2  2 ] ; % loaded node

for k=1:num_elem
    mesh.conecCell{ end+1, 1 } = [ 1 2 0 k k+1 ] ;
end

initialConds = {} ;

analysisSettings                    = {} ;
analysisSettings.methodName         = 'arcLength' ;
analysisSettings.deltaT             = 1 ;
# analysisSettings.incremArcLen       = [1e-4*ones(1,10) 1e-5*ones(1,3000) 1e-4*ones(1,1000)] ;
analysisSettings.incremArcLen       = [ deltac/10*ones(1,8) -deltac/10*ones(1,8)  deltac/10*ones(1,30) -deltac/10*ones(1,14) deltac/10*ones(1,180) +deltac/1000*ones(1,10) ] ;
analysisSettings.finalTime          = length(analysisSettings.incremArcLen) ;
analysisSettings.iniDeltaLamb       = 1 ;
analysisSettings.posVariableLoadBC  = 2 ;
analysisSettings.stopTolDeltau      = 1e-14 ;
analysisSettings.stopTolForces      = 1e-8 ;
analysisSettings.stopTolIts         = 30 ;
analysisSettings.ALdominantDOF      = [6+3 -1] ;

otherParams              = struct() ;
otherParams.problemName  = 'plastic_2dframe' ;
% otherParams.plots_format = 'vtk' ;

[ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

[matUs, loadFactorsMat, ~ ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;


girosUltimoNodo      = matUs((num_elem+1)*6,:) ;
descensosNodoCargado = matUs(6+3,:) ;
factorescarga        = loadFactorsMat(:,2) ;

u_elem1 = matUs(   [3 6 3+6 6+6],end);
u_elem2 = matUs( 6+[3 6 3+6 6+6],end);

[xs1, deformada1] = forma(u_elem1,.01,l1/2,l1);
[xs2, deformada2] = forma(u_elem2,0,0,l2);
xs2=xs2+l1;
% ----------------------------------------------------------------------------------

% Plots

lw = 2 ; ms = 1 ; plotfontsize = 14 ;

figure('Name','giros','NumberTitle','off') ;
grid on
plot( girosUltimoNodo , factorescarga, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#EDB120") ;

figure('Name','desplazamientos/load fact','NumberTitle','off') ;
grid on
plot( descensosNodoCargado, factorescarga, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#EDB120") ;

figure('Name','desplazamientos/tiempo','NumberTitle','off') ;
grid on
plot( descensosNodoCargado, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#EDB120") ;

figure('Name','load factors solo','NumberTitle','off') ;
grid on
plot( factorescarga, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#EDB120") ;

figure('Name','deformada','NumberTitle','off') ;
grid on, hold on
plot( xs1, deformada1, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#EDB120") ;
plot( xs2, deformada2, '-s', 'linewidth', lw, 'markersize', ms, "Color", "#EDB120") ;




stop
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

% plot(abs(rotations), abs(Mn1_validation), '-x' , 'linewidth', lw, 'markersize', ms, "Color", "#EDB120") ;
plot(abs(displacements), abs(Mn1_validation), '-x' , 'linewidth', lw, 'markersize', ms*6, "Color", "#0072BD") ;

plot(desp_numer(1:length(Mn_numer)), Mn_numer(1:length(Mn_numer)), '-x' , 'linewidth', lw, 'markersize', ms, "Color", "#7E2F8E") ;

% plot(abs(matdes(6,1:length(load_factors)-1)), Mn,'-x' , 'linewidth', lw, 'markersize', ms, "Color", "#D95319") ;
% plot(abs(matdes(4,1:length(load_factors)-1)), Mn, '-x' , 'linewidth', lw, 'markersize', ms, "Color", "#77AC30") ;

labx = xlabel('Generalized displacements in free node (m, rad)') ;
laby = ylabel('Bulk Moment at the first integration point (KN.m)') ;
legend('Semi Analytic (1 elem) [y]', 'ALGOL (1 elem) [y]', 'location', 'Southeast') ;

set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Cantilever Beam / Plasticity (validation)') ;

print('-f1','../../../Tesis/tex/imagenes/Load_factors.png','-dpng') ;
print('-f2','../../../Tesis/tex/imagenes/Validation.png','-dpng') ;