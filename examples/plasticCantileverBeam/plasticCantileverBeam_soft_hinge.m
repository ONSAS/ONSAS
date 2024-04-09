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

close all ;
clear ;
addpath( genpath( [ pwd '/../../src'] ) ) ;

% assumed XY plane

% -------------------------------------------
% scalar parameters
% material
E = 30000000 ;      % KN/m^2 KPa
kh1 = 29400 ;       % KN.m^2
kh2 = 273 ;
Ks = -18000 ;       % KN.m

nu = 0.3 ;          % Poisson's ratio

% geometry
l = 2.5 ;           % m
ty = 0.3 ;          % width cross section
tz = 0.4 ;          % height cross section

A = .4*.3 ;         % m^2
EI = 77650 ;        % KN.m^2
Inertia = EI/E ;    % m^4
Mc = 37.9 ;         % KN.m
My = 268 ;
Mu = 376 ;

% at the beginning..., there was no softening hinge
soft_hinge_boolean = false ;

% number of finite elements
num_elem = 1 ;

% /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ 

global historic_parameters

global arcLengthFlag % 1: cylindrical 2: jirasek
arcLengthFlag = 2 ;

global dominantDofs
dominantDofs = (num_elem+1)*6-3 ;

global scalingProjection
scalingProjection = -1 ;

global sizecmatrix
sizecmatrix = 6*(num_elem+1) ;

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
analysisSettings.incremArcLen       = [1e-3*ones(1,847) eps*ones(1,1)] ;
analysisSettings.finalTime          = length(analysisSettings.incremArcLen) ;
analysisSettings.iniDeltaLamb       = 1 ;
analysisSettings.posVariableLoadBC  = 2 ;
analysisSettings.stopTolDeltau      = 1e-8 ;
analysisSettings.stopTolForces      = 1e-8 ;
analysisSettings.stopTolIts         = 15 ;

otherParams              = struct() ;
otherParams.problemName  = 'plastic_2dframe' ;
otherParams.plots_format = 'vtk' ;

[matUs, loadFactorsMat, internalforces] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

girosUltimoNodo = matUs((num_elem+1)*6,:) ;
descensosUltimoNodo = matUs((num_elem+1)*6-3,:) ;
factorescarga = loadFactorsMat(:,2) ;

% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

lw = 1.5 ; ms = 1 ; plotfontsize = 14 ;

figure('Name','Cantilever Beam / Plasticity (load factors)','NumberTitle','off') ;
hold on, grid on
plot(abs(girosUltimoNodo), factorescarga,'-x' , 'linewidth', lw, 'markersize', ms, "Color", "#EDB120") ;
plot(abs(descensosUltimoNodo), factorescarga, '-x' , 'linewidth', lw, 'markersize', ms, "Color", "#0072BD") ;

labx = xlabel('Generalized displacements in free node (m, rad)') ; 
laby = ylabel('Forces') ;
legend('ONSAS (1 elem) \theta', 'ONSAS (1 elem) y', 'location','Southeast') ;
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Cantilever Beam / Plasticity (load factors)') ;