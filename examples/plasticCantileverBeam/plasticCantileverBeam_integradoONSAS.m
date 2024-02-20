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

% assumed XZ plane

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
ty = .3; % width cross section
tz = .4; % height cross section

A = .4*.3 ;       % m^2
EI = 77650 ;        % KN.m^2
Inercia = EI/E ;         % m^4
Mc = 37.9 ;         % KN.m
My = 268 ;
Mu = 374 ;

% -------------------------------------------


materials             = struct() ;
materials.modelName   = 'plastic-2Dframe' ;
materials.modelParams = [ E Mc My Mu kh1 kh2 Ks nu ]        ;

%materials.modelName  = 'elastic-rotEngStr' ;
%md and in the field `modelParams` a vector with the parameters of the Engineering Strain model is set
%materials.modelParams = [ E nu ] ;


disp('hola')
elements             = struct() ;
elements(1).elemType = 'node'  ;
%mdframe elements for modelling the blades
elements(2).elemType = 'frame' ;
elements(2).elemCrossSecParams = {'generic' ; [A 1 Inercia Inercia] };


boundaryConds                  = {} ;
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
%
boundaryConds(2).loadsCoordSys = 'global'         ;
boundaryConds(2).loadsBaseVals = [ 0 0 -1 0 0 0 ] ;
boundaryConds(2).loadsTimeFact = @(t) t     ;
%md

%mdThe coordinates of the nodes of the mesh are given by the matrix:
mesh             = {} ;
mesh.nodesCoords = [   0  0   0 ; ...
                       l  0   0 ] ;

mesh.conecCell = {} ;

mesh.conecCell{ 1, 1 } = [ 0 1 1   1   ] ; % nodo
mesh.conecCell{ 2, 1 } = [ 0 1 2   2   ] ; % nodo

mesh.conecCell{ 3, 1 } = [ 1 2 0   1 2   ] ;

initialConds = {} ;


analysisSettings               = {};
analysisSettings.methodName    = 'newtonRaphson' ;

analysisSettings.deltaT        =   1  ;
analysisSettings.finalTime     =   16    ;
analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   15   ;
%md

otherParams              = struct();
otherParams.problemName  = 'plastic_2dframe';
otherParams.plots_format = 'vtk' ;

disp('hola')

[matUs, loadFactorsMat ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;


disp('chau')
