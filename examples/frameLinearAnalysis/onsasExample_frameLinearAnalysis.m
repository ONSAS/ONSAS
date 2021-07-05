%md Linear Beam Element example
%md
%md linear analysis
close all, clear all        ;
addpath( genpath( [ pwd '/../../src' ] ) ) ; % add ONSAS directory to path

% ----------------------------------------------------------------------
% scalar auxiliar parameters
E = 210e6 ; L = 5 ; nu = 0.3 ;  rho = 0 ; b = 0.3 ; P = -5 ;

% ----------------------------------------------------------------------
% MEBI parameters: Material-Element-BoundaryConditions-InitialConditions

% Materials
materials.hyperElasModel = {'linearElastic'} ;
materials.hyperElasParams = { [E, nu] } ;

% Elements
elements.elemType  = { 'node', 'frame' } ;
elements.elemTypeGeometry  = { [], [2, b, b] } ;
elements.elemTypeParams = { [], 1 };

% BoundaryConditions
% Loads
boundaryConds.loadsCoordSys = { [] ; 'global' } ;
boundaryConds.loadsTimeFact = { [] ; @(t) t } ;
boundaryConds.loadsBaseVals = { [] ; [ 0 0 0 0 P 0] } ;
% Supports
boundaryConds.imposDispDofs = { [ 1 2 3 4 5 6 ] ; [] } ;
boundaryConds.imposDispVals = { [ 0 0 0 0 0 0 ] ; 0 } ;

% InitialConditions
initialConds = struct() ;


% Mesh nodes

mesh.nodesCoords = ...
					[ 0		0	0		; ...
						L 	0 0 	; ...
						2*L 0 0 	] ;

% Conec Cell

mesh.conecCell = { } ;
% Node elements
mesh.conecCell{1,1} = [ 0 1 1 0 1 ] ;
mesh.conecCell{2,1} = [ 0 1 1 0 3 ] ;
mesh.conecCell{3,1} = [ 0 1 2 0 2 ] ;
% Frame elements
mesh.conecCell{4,1} = [ 1 2 0 0 1 2 ] ;
mesh.conecCell{5,1} = [ 1 2 0 0 2 3 ] ;

% Analysis settings
analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.deltaT        = 1 ;
analysisSettings.finalTime     = 1 ;
analysisSettings.stopTolDeltau = 1e-6 ;
analysisSettings.stopTolForces = 1e-6 ;
analysisSettings.stopTolIts    = 1 ;

otherParams.problemName = 'linearBeamElement_test' ;

A = b^2;
I = b^4/12 ;

axial = E*A/L *2;
bending = E*I/L^3 *4*L^2*2 ;

flecha = P*(2*L)^3/(192*E*I)


[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

verifBoolean = abs( flecha - matUs(11,2) ) < 1e-4
% ----------------------------------------------------------------------
