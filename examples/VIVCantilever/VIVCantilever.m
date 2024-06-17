%md# WOM VIV
%
% add ONSAS path
close all, if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear all, end
%
addpath( genpath( [ pwd '/../../src'] ) );
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
% Numerical parameters
%
dt = 0.0028 ; finalTime = 10*dt; numElements = 3;
%--------
% Fluid properties
%
nuFluid = 1e-6; rhoFluid = 1000;  
%
nameFuncVel = 'windVelUniform'; 
nameLiftFunc = 'liftCoefVIV'; 
nameDragFunc = 'dragCoef'; 
cL0 = feval(nameLiftFunc, 0, 1); cD0 = feval(nameDragFunc, 0, 1); 

% material
E = 5e10;
rho = 2*rhoFluid; 
nu = .3; 
cu = 0; 
%
l = 1 ; d = 0.001; I = pi * d^4 / 64 ;
%
St = 0.2; 
%--------
% With WOMV3
    % uzsol = 1.0e-07 *[  0    0         0         0         0         0         0            0         0         0         0;
    %                     0    0.0012    0.0058    0.0144    0.0264    0.0409    0.0569       0.0738    0.0907    0.1073    0.1232;
    %                     0    0.0011    0.0052    0.0126    0.0227    0.0347    0.0485       0.0642    0.0817    0.1010    0.1216;
    %                     0    0.0010    0.0050    0.0127    0.0229    0.0349    0.0491       0.0659    0.0852    0.1060    0.1275];
% With WOMV4
uzsol = 1.0e-07 *[ 0         0         0         0         0         0         0        0         0         0         0;
                   0   -0.0002   -0.0008   -0.0012   -0.0008    0.0004    0.0020        0.0035    0.0044    0.0042    0.0026;
                   0    0.0005    0.0023    0.0050    0.0082    0.0110    0.0127        0.0131    0.0131    0.0138    0.0153;
                   0    0.0011    0.0047    0.0111    0.0214    0.0356    0.0524        0.0719    0.0970    0.1284    0.1637];
%
% materials
%
materials.modelParams = [ E nu ] ;
materials.density         = rho      ;
materials.modelName  = 'elastic-rotEngStr' ; 
%
elements = struct();
% node
elements(1).elemType = 'node' ;
% hydro frame
elements(2).elemType = 'frame' ;
% cross section params
elements(2).elemCrossSecParams = {'circle' ; d } ;
elements(2).aeroCoefFunctions = {nameDragFunc, nameLiftFunc, []};
elements(2).chordVector = [0 0 -d];
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
%mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEB parameters, followed by the indexes of nodes that compose the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
mesh.conecCell = { } ;
%md then the first welded node is defined with material (M) zero since nodes don't have material, the first element (E) type (the first entry of the `elements` struct), and (B) is the first entry of the the `boundaryConds` struct. For (I) no non-homogeneous initial condition is considered (then zero is used) and finally the node is assigned:
mesh.conecCell{ 1, 1 } = [ 0 1 1   1 ] ;
%md Next the frame elements MEB parameters are set. The frame material is the first material of `materials` struct, then $1$ is assigned. The second entry of the `elements` struct correspond to the frame element employed, so $2$ is set. Finally no BC and no IC is required for this element, then $0$ is used.  Consecutive nodes build the element so then the `mesh.conecCell` is:
for i=1:numElements,
  mesh.conecCell{ i+1,1 } = [ 1 2 0  i i+1 ] ;
end
% Initialize qvect
qvect =  zeros(numElements*2,round(finalTime/dt)+1);
qvect(1:2:end,1) = (2*rand(numElements, 1)-1)*0.001 ;
% fluid properties
analysisSettings.fluidProps = {rhoFluid; nuFluid; nameFuncVel} ;
%If drag reconfiguration then `analysisSettings.geometricNonLinearAero` should be set to `true`. Also it if is `false` then the lift direction will be constant.
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
otherParams = struct() ;
otherParams.nodalDispDamping = cu ;
otherParams.problemName      = strcat('VIVTest') ;
%
% Run ONSAS
%
[ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%
%mdAfter that the structs are used to perform the numerical time analysis
matUs = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;

% Extract numerical solution
uz = matUs(5:6:end, :);
norm(uz - uzsol)      ;
%save('testSolution', 'uzTest')
if length(uz) == length(uzsol)
    verifBoolean = norm(uz - uzsol) < 6e-08
else 
    verifBoolean = 0;
end
%-----------------------------------------------------------------------------------------------------
