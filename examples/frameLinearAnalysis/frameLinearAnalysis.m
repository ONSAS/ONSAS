%md# A Frame Linear Analysis Example
%md
%md## Previous definitions
close all, clear all;
addpath( genpath( [ pwd '/../../src' ] ) ) ; % add ONSAS directory to path
%md
%md scalar auxiliar parameters
E  = 210e6 ; nu = 0.3 ; %
ty = 0.3   ; tz = 0.6 ; % cross-section widths
L1 = 2     ; L2 = 1.0 ; %
P = 1e3    ;            % applied nodal load
%md
%md## MEBI parameters: Material-Element-BoundaryConditions-InitialConditions
%md
%md### Materials
materials(1).hyperElasModel = 'linearElastic' ;
materials(1).hyperElasParams = [ E, nu] ;
%md### Elements
elements(1).elemType  = 'node'  ;
elements(2).elemType  = 'frame' ;
elements(2).elemCrossSecParams = { 'rectangle'; [ ty tz ] };
elements(2).massMatType     =  'consistent' ; 
%md### BoundaryConditions
% Supports
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
% Loads
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 -P 0] ;
%md### InitialConditions
%md empty struct
initialConds = struct() ;
%md
%md## Mesh
%md Mesh nodes
mesh.nodesCoords = ...
					[ 0 	0	 0	; ...
						L1  0  0 	; ...
						L1  L2 0 	] ;
%md
%md Conec Cell
mesh.conecCell = { } ;
%md nodes
mesh.conecCell{1, 1 } = [ 0 1 1 0  1 ] ;
mesh.conecCell{2, 1 } = [ 0 1 2 0  3 ] ;
%md and frame elements
mesh.conecCell{3, 1 } = [ 1 2 0 0  1 2 ] ;
mesh.conecCell{4, 1 } = [ 1 2 0 0  2 3 ] ;
%md
%md Analysis settings
analysisSettings = struct() ;

otherParams.problemName = 'linearBeamElement' ;
otherParams.plotsFormat = 'vtk' ;

[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

mu = E /( 2*(1+nu) ) ;
a = .5 * max( elements(2).elemCrossSecParams{2} ) ;
b = .5 * min( elements(2).elemCrossSecParams{2} ) ;
J = a * b^3 * ( 16/3 - 3.36 * b/a * ( 1 - b^4 / ( 12*a^4 ) ) ) ;
Iyy = ty * tz^3 / 12.0 ;

Mt = P*L2 ;

analyThetax1 = - L1 * Mt / ( mu * J ) ;
numerThetax1 = matUs( 8, 2 ) ;

analyDefl = -P*L1^3/ ( 3*E*Iyy ) +analyThetax1*L2 - P*L2^3/ ( 3*E*Iyy ) ;
numerDefl = matUs(6*2+5,2) ;

verifBoolean = ( abs( numerThetax1 - analyThetax1 ) < 1e-8 ) && ...
               ( abs( numerDefl     - analyDefl     ) < 1e-8 ) ;
