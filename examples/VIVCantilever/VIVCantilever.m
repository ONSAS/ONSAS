% WOM VIV Example
%
% add ONSAS path
close all, if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear all, end
addpath( genpath( [ pwd '/../../src'] ) );
%
% ---------------------
% scalar parameters
% fluid
nuFluid = 1e-6;
rhoFluid = 1000;  
% solid
E = 5e10;
rho = 2*rhoFluid; 
nu = .3; 
%
nodalDamping = 0; 
dt = 0.0028 ; finalTime = 1000*dt; numElements = 3;
l = 1 ; d = 0.001;
St = 0.2; 

nameFuncVel = 'windVelUniform'; 
nameLiftFunc = 'liftCoefVIV'; 
nameDragFunc = 'dragCoef'; 
% ---------------------

% material (structural solid)
%
materials = struct();
materials.modelParams = [ E nu ] ;
materials.density         = rho      ;
materials.modelName  = 'elastic-rotEngStr' ; 

% elements
%
elements = struct();
elements(1).elemType = 'node' ;
elements(2).elemType = 'frame' ;
elements(2).elemCrossSecParams = {'circle' ; d } ;
elements(2).aeroCoefFunctions = { nameDragFunc, nameLiftFunc, []};
elements(2).chordVector = [0 0 -d];
elements(2).massMatType = 'consistent' ; 

% boundaryConds
%
% md The first and unique BC corresponds to a welded condition for a cantilever beam
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;

% initialConds
%
initialConds = struct() ;

% mesh 
%
mesh.nodesCoords = [  zeros(numElements+1,1) (0:(numElements))'*l/numElements zeros(numElements+1,1) ] ;
mesh.conecCell = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1   1 ] ;
for i=1:numElements,
  mesh.conecCell{ i+1,1 } = [ 1 2 0  i i+1 ] ;
end

% Numerical parameters
%
analysisSettings.methodName     = 'newmark'  ;  
analysisSettings.finalTime      = finalTime  ;
analysisSettings.deltaT         = dt         ; 
analysisSettings.stopTolIts     = 15         ; 
analysisSettings.stopTolDeltau  = 1e-10      ; 
analysisSettings.stopTolForces  = 1e-5       ; 

analysisSettings.fluidProps = {rhoFluid; nuFluid; nameFuncVel} ;
%If drag reconfiguration is considered then `analysisSettings.geometricNonLinearAero` should be set to `true`. On the other hand it this is set as `false` then the lift direction will be constant and given by the lift direction at the reference configuration.
analysisSettings.geometricNonLinearAero = true;
analysisSettings.booleanSelfWeight = false ;
analysisSettings.VIVBool        = true;

global qvect; %VIV boolean is called inside hydroFrameForces 
global uniformUdot; %constantLiftDir is called inside hydroFrameForces
uniformUdot = false; 

% Initialize qvect
qvect =  zeros(numElements*2,round(finalTime/dt)+1);
qvect(1:2:end,1) = (2*rand(numElements, 1)-1)*0.001 ;

% otherParams
otherParams = struct() ;
otherParams.nodalDispDamping = nodalDamping ;
otherParams.problemName      = 'vivCantilever' ;
otherParams.plots_format = 'vtk' ;

% Run ONSAS
%
[ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%
% mdAfter that the structs are used to perform the numerical time analysis
matUs = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;
%--------

% Fluid properties
cL0 = feval(nameLiftFunc, 0, 1); cD0 = feval(nameDragFunc, 0, 1); 
%
%
% --------
uzsol = 1.0e-07 *[ 0         0         0         0         0         0         0        0         0         0         0;
                   0   -0.0002   -0.0008   -0.0012   -0.0008    0.0004    0.0020        0.0035    0.0044    0.0042    0.0026;
                   0    0.0005    0.0023    0.0050    0.0082    0.0110    0.0127        0.0131    0.0131    0.0138    0.0153;
                   0    0.0011    0.0047    0.0111    0.0214    0.0356    0.0524        0.0719    0.0970    0.1284    0.1637];
%


% Extract numerical solution
uz = matUs(5:6:end, :);
verifBoolean = norm(uz - uzsol) < 6e-08 
