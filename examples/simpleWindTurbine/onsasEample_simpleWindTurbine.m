%md# Wind turbine example
close all, clear all
addpath( genpath( [ pwd '/../../src'] ) );

% scalar parameters
E = 70e9 ;  nu = 0.3 ; rho = 500 ; G = E / (2 * (1+nu)) ;
l = 5 ; d = 0.3;
numElementsBlade = 1 ;

% materials
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu ]         ;
materials.density         = rho              ;

% elements
%----------------------------
elements(1).elemType         = 'node'  ;
% formulation type and number of gauss integration points
formulCase = 2 ;
numGaussPoints = 2 ;
% first type of blade
elements(2).elemType         = 'frame' ;
elements(2).elemTypeGeometry = [2 d d] ;
elements(2).elemTypeAero     = [0 0 d numGaussPoints formulCase ] ;
elements(2).userDragCoef     = 'dragCoefS809'   ;
elements(2).userLiftCoef     = 'liftCoefS809'   ;
% second type of blade
elements(3).elemType         = 'frame' ;
elements(3).elemTypeGeometry = [2 d d] ;
elements(3).elemTypeAero     = [0 sqrt(d) sqrt(d) numGaussPoints formulCase ] ;
elements(3).userDragCoef     = 'dragCoefS809'   ;
elements(3).userLiftCoef     = 'liftCoefS809'   ;
% boundaryConds
%----------------------------
% The elements are submitted to two different BC settings. The first BC corresponds to a welded condition (all 6 dofs set to zero)
boundaryConds(1).imposDispDofs = [ 1 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 ] ;
%
% initial Conditions
%----------------------------
% homogeneous initial conditions are considered, then an empty struct is set:
initialConds = struct() ;
%
% mesh parameters
%----------------------------
%The coordinates of the nodes of the mesh are given by the matrix:
localAxialBladeCords = ( 0:( numElementsBlade ) )'*l/ numElementsBlade ;
nodesLocalBladeZ = [ zeros( numElementsBlade +1,2) -localAxialBladeCords ] ;
nodesLocalBlade120 = [ zeros( numElementsBlade +1,1), +sin( deg2rad(120) )*localAxialBladeCords,  -cos( deg2rad(120) )*(localAxialBladeCords) ] ;
nodesLocalBlade120(1,:) = [] ;
nodesLocalBlade240 = [ zeros( numElementsBlade +1,1) -sin( deg2rad(60) )*localAxialBladeCords +cos( deg2rad(60) )*(localAxialBladeCords) ] ;
nodesLocalBlade240(1,:) = [] ;
%The final nodes coordinates matrix is:
% the mesh conecitvity cell is
mesh.nodesCoords = [   nodesLocalBladeZ;  nodesLocalBlade240; nodesLocalBlade120; ] ;
mesh.conecCell = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1 ] ;
%The conec cell is assamble as:
for i=1:numElementsBlade
    conecElemMatrix(i,:) = [ 1 2 0 0  i i+1 ] ;
    mesh.conecCell{ i+1,1 } = [ 1 2 0 0  i i+1 ] ;
  if i == 1 
    conecElemMatrix( i + numElementsBlade  , : ) = [ 1 3 0 0  1   numElementsBlade + 1 + i ] ;
    conecElemMatrix( i + 2*numElementsBlade, : ) = [ 1 2 0 0  1 2*numElementsBlade + 1 + i ] ;
    mesh.conecCell{ i  + numElementsBlade + 1 , 1 } = [ 1 3 0 0  1     numElementsBlade + 1 + i ] ;
    mesh.conecCell{ i  + 2*numElementsBlade +1, 1 } = [ 1 2 0 0  1   2*numElementsBlade + 1 + i ] ;
  elseif i > 1
    conecElemMatrix( i + numElementsBlade,   : ) = [ 1 2 0 0    numElementsBlade + i    numElementsBlade + i + 1  ] ;
    conecElemMatrix( i + 2*numElementsBlade, : ) = [ 1 3 0 0  2*numElementsBlade + i  2*numElementsBlade + i + 1  ] ;
    mesh.conecCell{ i  +   numElementsBlade +1 , 1 } = [ 1 2 0 0    numElementsBlade + i   numElementsBlade + i + 1  ] ;
    mesh.conecCell{ i  + 2*numElementsBlade +1 , 1 } = [ 1 2 0 0  2*numElementsBlade + i 2*numElementsBlade + i + 1  ] ;
  end
end
% analysisSettings
%----------------------------
analysisSettings.methodName    = 'newtonRaphson'                    ;
analysisSettings.finalTime      =   1                               ;
analysisSettings.deltaT        =   analysisSettings.finalTime / 2   ;
analysisSettings.stopTolDeltau =   1e-6                             ;
analysisSettings.stopTolForces =   1e-6                             ;
analysisSettings.stopTolIts    =   30                               ;
analysisSettings.booleanSelfWeight = false                          ;
analysisSettings.userWindVel = 'windVelDynamic';


% otherParams
%----------------------------
otherParams.problemName      = 'onsasExample_simpleWindTurbine' ;
otherParams.plotsFormat      = 'vtk' ;

% Execute ONSAS
%----------------------------
[ matUs, ~ ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
