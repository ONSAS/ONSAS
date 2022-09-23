%md# WOM VIV
%
% add ONSAS path
close all, clear all ; addpath( genpath( [ pwd '/../../src'] ) );
%tic ;
%
% declare global variables
%
global qvect; %VIV boolean is called inside hydroFrameForces 
global VIVBool; %VIV boolean is called inside hydroFrameForces
global constantLiftDir; %constantLiftDir is called inside hydroFrameForces
global uniformUdot; %constantLiftDir is called inside hydroFrameForces
VIVBool = true;
if VIVBool
  constantLiftDir = false; uniformUdot = false; 
end
%
% Load Parameters
parametersSet
%
% materials
%
materials.hyperElasParams = [ E nu ] ;
materials.density         = rho      ;
materials.hyperElasModel  = '1DrotEngStrain' ; 
%
% node
elements(1).elemType = 'node' ;
% hydro frame
elements(2).elemType = 'frame' ;
% cross section params
elements(2).elemCrossSecParams{1,1} = 'circle' ;
elements(2).elemCrossSecParams{2,1} = d          ;
% hydro cross-section props
numGaussPoints  = 4 ;
computeAeroTangMatrix = false ;
elements(2).aeroCoefs   = {nameDragFunc; nameLiftFunc; [] }   ;
%  chord vector and gauss points

elements(2).elemTypeAero = [0 0 -d numGaussPoints computeAeroTangMatrix ] ; % [chordVec1 chordVec2 chordVec3 numGauss  ]
% mass element formulation
elements(2).massMatType = 'consistent' ; 
%
% boundaryConds
%
%md The first and unique BC corresponds to a welded condition for a cantilever beam
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
%
% initialConds
%
%md homogeneous initial conditions are considered, then an empty struct is set:
initialConds = struct() ;
%
% mesh 
%
%mdThe coordinates of the mesh nodes are given by the matrix:
mesh.nodesCoords = [  zeros(numElements+1,1) (0:(numElements))'*l/numElements zeros(numElements+1,1) ] ;
%mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of nodes that compose the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
mesh.conecCell = { } ;
%md then the first welded node is defined with material (M) zero since nodes don't have material, the first element (E) type (the first entry of the `elements` struct), and (B) is the first entry of the the `boundaryConds` struct. For (I) no non-homogeneous initial condition is considered (then zero is used) and finally the node is assigned:
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1 ] ;
%md Next the frame elements MEBI parameters are set. The frame material is the first material of `materials` struct, then $1$ is assigned. The second entry of the `elements` struct correspond to the frame element employed, so $2$ is set. Finally no BC and no IC is required for this element, then $0$ is used.  Consecutive nodes build the element so then the `mesh.conecCell` is:
for i=1:numElements,
  mesh.conecCell{ i+1,1 } = [ 1 2 0 0  i i+1 ] ;
end
% Initialize qvect
qvect =  zeros(numElements*2,round(finalTime/dt)+1);
% fluid properties
analysisSettings.fluidProps = {rhoFluid; nuFluid; nameFuncVel} ;
%If drag reconfiguration then analysisSettings.geometricNonLinearAero  = true!! also it if is false then the lift direction will be constant
analysisSettings.geometricNonLinearAero = true;
% time parameters
analysisSettings.finalTime   = finalTime ;
analysisSettings.deltaT      = dt        ; 
% numerical parameters
analysisSettings.methodName     =   'newmark' ;  
analysisSettings.stopTolIts     =   20        ; 
analysisSettings.stopTolDeltau  =   1e-10      ; 
analysisSettings.stopTolForces  =   1e-5      ; 
% load settings
analysisSettings.booleanSelfWeight = false ;
% otherParams
otherParams.nodalDispDamping = cu ;
otherParams.problemName      = strcat('VIVTest') ;
%
% Run ONSAS
%
[ matUs, ~ ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ; 
% Extract numerical solution
uz = matUs(5:6:end, :) 
uz - uzsol
norm(uz - uzsol)
%save('testSolution', 'uzTest')
if length(uz) == length(uzsol)
    verifBoolean = norm(uz - uzsol) < 6e-08
else 
    verifBoolean = 0;
end
%-----------------------------------------------------------------------------------------------------
