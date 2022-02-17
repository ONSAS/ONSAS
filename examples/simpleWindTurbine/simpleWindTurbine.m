%md# Wind turbine example
close all, clear all
addpath( genpath( [ pwd '/../../src'] ) );

% scalar parameters
E = 70e9 ;  nu = 0.3 ; rho = 500 ; G = E / (2 * (1+nu)) ;
l = 5 ; d = 0.3;

% materials
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu ]        ;
materials.density         = rho             ;

% elements
elements(1).elemType = 'node'  ;
elements(2).elemType = 'frame' ;
elements(2).elemTypeGeometry = [2 d d] ;

% boundaryConds
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;

% initial Conditions
initialConds = struct() ;

% mesh parameters
mesh.nodesCoords = [ 0              0                0  ; ...
                     0  l*cos(      0 )  l*sin(      0 ) ; ...
                     0  l*cos( 2*pi/3 )  l*sin( 2*pi/3 ) ; ...
                     0  l*cos( 4*pi/3 )  l*sin( 4*pi/3 ) ] ;

mesh.conecCell         = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 0   1   ] ;
mesh.conecCell{ 2, 1 } = [ 1 2 0 0   1 2 ] ;
mesh.conecCell{ 3, 1 } = [ 1 2 0 0   1 3 ] ;
mesh.conecCell{ 4, 1 } = [ 1 2 0 0   1 4 ] ;



% Static case

% analysisSettings
analysisSettings.methodName    = 'newtonRaphson'                    ;
analysisSettings.finalTime      =   1                                ;
analysisSettings.deltaT        =   analysisSettings.finalTime / 2    ;
analysisSettings.stopTolDeltau =   1e-6                             ;
analysisSettings.stopTolForces =   1e-6                             ;
analysisSettings.stopTolIts    =   30                               ;
analysisSettings.booleanSelfWeight = false                           ;

% otherParams
%----------------------------
otherParams.problemName      = 'simpleWindTurbine' ;
otherParams.plotsFormat      = 'vtk' ;


[ matUsStatic, ~ ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
