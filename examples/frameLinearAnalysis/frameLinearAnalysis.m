% md# A Frame Linear Analysis Example
% md
% md## Previous definitions
close all, if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear all, end
addpath( genpath( [ pwd '/../../src' ] ) ) ; % add ONSAS directory to path
% md
% md scalar auxiliar parameters
E  = 210e9 ; nu = 0.3 ; %
ty = 0.1   ; tz = 0.2 ; % cross-section widths
L1 = 2     ; L2 = 1.5 ; %
Pz = 2e3   ; Py = 1e3 ; % applied nodal loads
% md
% md## MEB parameters: Material-Element-BoundaryConditions
% md
% md### Materials
materials = struct();
materials(1).modelName = 'elastic-linear'  ;
materials(1).modelParams = [ E, nu] ;
% md### Elements
elements = struct();
elements(1).elemType  = 'node'  ;
elements(2).elemType  = 'frame' ;
elements(2).elemCrossSecParams = { 'rectangle'; [ ty tz ] };
elements(2).massMatType     =  'consistent' ; 
% md### BoundaryConditions
% Supports
boundaryConds = struct();
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
% Loads
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsBaseVals = [ 0 0 Py 0 -Pz 0] ;
% md
% md## Mesh
% md Mesh nodes
mesh = struct();
mesh.nodesCoords = [ 0   0  0 ; ...
					 L1  0  0 ; ...
					 L1  L2 0 ] ;
% md
% md Conec Cell
mesh.conecCell = { } ;
% md nodes
mesh.conecCell{1, 1 } = [ 0 1 1   1 ] ;
mesh.conecCell{2, 1 } = [ 0 1 2   3 ] ;
% md and frame elements
mesh.conecCell{3, 1 } = [ 1 2 0   1 2 ] ;
mesh.conecCell{4, 1 } = [ 1 2 0   2 3 ] ;
% md
% md### InitialConditions
% md empty struct
initialConds = struct() ;
%
% md Analysis settings
analysisSettings = struct() ;
%
otherParams = struct() ;
otherParams.problemName = 'frameLinearAnalysis' ;
otherParams.plots_format = 'vtk' ;
%
[ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%
% mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;
% md
% md the report is generated
outputReport( modelProperties.outputDir, modelProperties.problemName )

mu = E /( 2*(1+nu) ) ;
a = .5 * max( elements(2).elemCrossSecParams{2} ) ;
b = .5 * min( elements(2).elemCrossSecParams{2} ) ;
J = a * b^3 * ( 16/3 - 3.36 * b/a * ( 1 - b^4 / ( 12*a^4 ) ) ) ;
Iyy = ty * tz^3 / 12.0 ;

Mt = Pz*L2 ;

analyThetax1 = - L1 * Mt / ( mu * J ) ;
numerThetax1 = matUs( 8, 2 ) ;

analyDeflz = -Pz*L1^3/ ( 3*E*Iyy ) +analyThetax1*L2 - Pz*L2^3/ ( 3*E*Iyy ) ;
numerDeflz = matUs(6*2+5,2) ;

verifDisps = ( abs( numerThetax1 - analyThetax1 ) < 1e-8 ) && ...
               ( abs( numerDeflz     - analyDeflz     ) < 1e-8 ) ;

analyFintElem1 = [ -L1*Pz  ] ;  
numerFintElem1 = modelSolutions{end}.localInternalForces(1).My ;

analyFintElem2 = [ -Py ] ;
numerFintElem2 = modelSolutions{end}.localInternalForces(2).Nx ;

verifFints = ( norm( numerFintElem1 - analyFintElem1 ) < ( 1e-8 * norm(analyFintElem1) ) ) && ...
             ( norm( numerFintElem2 - analyFintElem2 ) < ( 1e-8 * norm(analyFintElem2) ) ) ;
   
verifBoolean = verifDisps && verifFints