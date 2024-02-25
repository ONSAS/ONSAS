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

close all, if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear, end
addpath( genpath( [ pwd '/../../src'] ) ) ;

% assumed XY plane

% -------------------------------------------
% scalar parameters
% material
E = 30000000 ;      % KN/m^2 KPa
nu   = 0.3 ;
kh1 = 29400 ;       % KN.m^2
kh2 = 272 ;
Ks = -18000 ;       % KN.m

% geometry
l = 2.5 ;           % m
ty =0.3 ;           % width cross section
tz =0.4 ;           % height cross section

A = .4*.3 ;         % m^2
EI = 77650 ;        % KN.m^2
Inercia = EI/E ;    % m^4
Mc = 37.9 ;         % KN.m
My = 268 ;
Mu = 374 ;

% at the beginning..., there was no softening hinge
soft_hinge_boolean = false ;

% number of finite elements
num_elem = 1 ;

% -------------------------------------------

global historico_params

historico_params=[];

materials             = struct() ;
materials.modelName   = 'plastic-2Dframe' ;
materials.modelParams = [ E Mc My Mu kh1 kh2 Ks nu ] ;

%materials.modelName  = 'elastic-rotEngStr' ;
%md and in the field `modelParams` a vector with the parameters of the Engineering Strain model is set
%materials.modelParams = [ E nu ] ;

elements             = struct() ;
elements(1).elemType = 'node'  ;
%mdframe elements for modelling the blades
elements(2).elemType = 'frame' ;
elements(2).elemCrossSecParams = {'generic' ; [A 1 Inercia Inercia] };

boundaryConds                  = {} ;
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
%
boundaryConds(2).imposDispDofs = [ 2 4 5] ;
boundaryConds(2).imposDispVals = [ 0 0 0 ] ;
boundaryConds(2).loadsCoordSys = 'global'         ;
boundaryConds(2).loadsBaseVals = [ 0 0 -1 0 0 0 ] ;
boundaryConds(2).loadsTimeFact = @(t) t     ;
%md
boundaryConds(3).imposDispDofs = [ 2 4 5] ;
boundaryConds(3).imposDispVals = [ 0 0 0 ] ;

%mdThe coordinates of the nodes of the mesh are given by the matrix:
mesh             = {} ;
xs = linspace(0,l,num_elem+1);
mesh.nodesCoords =  [  xs' zeros(num_elem+1, 2) ] ;

mesh.conecCell = {} ;

mesh.conecCell{ 1, 1 } = [ 0 1 1   1   ] ; % nodo

if num_elem>1
  for k=2:num_elem
    mesh.conecCell{ end+1, 1 } = [ 0 1 3   k   ] ;
  end
end

for k=1:num_elem
  mesh.conecCell{ end+1, 1 } = [ 1 2 0   k k+1   ] ;
end
mesh.conecCell{ end+1, 1 } = [ 0 1 2   num_elem+1   ] ; % nodo cargado

initialConds = {} ;

%~ analysisSettings               = {} ;
%~ analysisSettings.methodName    = 'newtonRaphson' ;
%~ analysisSettings.deltaT        =   5  ;
%~ analysisSettings.finalTime     =   25 ;
%~ analysisSettings.stopTolDeltau =   1e-8 ;
%~ analysisSettings.stopTolForces =   1e-8 ;
%  analysisSettings.stopTolIts    =   15   ;
%md
analysisSettings                    = {}            ;
analysisSettings.methodName         = 'arcLength'   ;
analysisSettings.deltaT             = 1             ;
analysisSettings.incremArcLen       = 1e-4          ;
analysisSettings.finalTime          = 1000          ;
analysisSettings.iniDeltaLamb       = 1             ;
analysisSettings.posVariableLoadBC  = 2             ;
analysisSettings.stopTolDeltau      = 1e-8          ;
analysisSettings.stopTolForces      = 1e-8          ;
analysisSettings.stopTolIts         = 15            ;

otherParams              = struct();
otherParams.problemName  = 'plastic_2dframe';
otherParams.plots_format = 'vtk' ;

[matUs, loadFactorsMat ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

girosUltimoNodo = matUs((num_elem+1)*6,:);
descensosUltimoNodo = matUs((num_elem+1)*6-3,:);
factorescarga = loadFactorsMat(:,2) ;

lw = 2.5 ; ms = 0.5 ; plotfontsize = 16 ;

figure('Name','Cantilever Beam / Plasticity','NumberTitle','off');
hold on, grid on
plot(abs(girosUltimoNodo), factorescarga,'b-x' , 'linewidth', lw, 'markersize', ms, "Color", "#EDB120") ;
plot(abs(descensosUltimoNodo), factorescarga, 'k-o' , 'linewidth', lw, 'markersize', ms, "Color", "#0072BD") ;
labx = xlabel('Generalized displacements in free node (m, rad)'); 
laby = ylabel('Lambda') ;
legend('Degree of Freedom y','Degree of Freedom \theta','location','Southeast') ;
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Cantilever Beam / Plasticity') ;