%md# Plate-element cantilever model
%md
close all, if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear all, end
addpath( genpath( [ pwd '/../../src'] ) );
%md
%md## Scalars
E = 200e9 ;
nu = 0.0 ;
tz = .05 ;
qx  = 1e3 ; % kN/m^2 
qz  = 1e3 ; % kN/m^2 
%md
Ly = .5;
Lx = 1 ;
%md## Numerical solution using plate elements
%md
%md### Materials
%md
materials                    = struct() ;
materials(1).modelName  = 'elastic-linear' ;
materials(1).modelParams = [ E nu ] ;
%md
%md### Elements
%md
elements             = struct() ;
elements(1).elemType = 'edge' ;
elements(1).elemCrossSecParams = tz         ;
elements(2).elemType = 'triangle-plate' ;
elements(2).elemCrossSecParams = {'thickness', tz } ;
%md
%md### Boundary conditions
%md
boundaryConds                  = struct() ;
boundaryConds(1).imposDispDofs =  [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals =  [ 0 0 0 0 0 0 ] ;
%
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsTimeFact = @(t) t  ;
boundaryConds(2).loadsBaseVals = [0 0 0 0 -qz 0 ] ;
%
boundaryConds(3).loadsCoordSys = 'global' ;
boundaryConds(3).loadsTimeFact = @(t) t  ;
boundaryConds(3).loadsBaseVals = [ qx 0 0 0 0 0 ] ;
%md
%md### mesh
%md
mesh = struct() ;
base_dir='';
if strcmp( getenv('TESTS_RUN'),'yes') && isfolder('examples'),
  base_dir=['.' filesep 'examples' filesep  'cantileverPlate' filesep];
end
[ mesh.nodesCoords, mesh.conecCell ] = meshFileReader( [base_dir 'geometry_cantileverPlate.msh'] ) ;
assert( max(mesh.nodesCoords(:,1)) == Lx && max(mesh.nodesCoords(:,2)) == Ly )

%md### Initial conditions
initialConds                  = struct() ;
%md
%md#### Analysis settings
analysisSettings               = struct() ;
analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.deltaT        =   1   ;
analysisSettings.finalTime     =   1  ;
analysisSettings.stopTolDeltau =   1e-10    ;
analysisSettings.stopTolForces =   1e-10    ;
analysisSettings.stopTolIts    =   10      ;
%md
%md#### OtherParams
%md The nodalDispDamping is added into the model using:
otherParams                  = struct() ;
%md The name of the problem is:
%md
otherParams.problemName  = 'cantileverPlate' ;
otherParams.plots_format = 'vtk' ;
%md
%md Execute ONSAS and save the results:
[ modelInitSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%
%mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions ] = ONSAS_solve( modelInitSol, modelProperties, BCsData ) ;
%md
%md## verification
nelem=size(modelProperties.Conec,1);
matSolic = getInternalForces( modelSolutions{end}.localInternalForces, 1:nelem, {'Mx','My','Mxy'} );
numer_maxMx = max(max(matSolic));
numer_wmax = min(matUs(5:6:end)) ;


%md
analy_maxMx = qz*Lx/2 ;
qlin = qz*Ly;  I = Ly*tz^3/12;
analy_wmax = -qlin*Lx^4/(8*E*I)  ;

elements(2).elemType           = 'triangle';
elements(2).elemTypeParams     = 2         ;
elements(2).elemCrossSecParams = tz         ;

[ modelInitSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%
%mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions ] = ONSAS_solve( modelInitSol, modelProperties, BCsData ) ;

numer_dxmax = max(matUs(1:6:end)) ;

analy_dxmax = qx*Lx/E ;

%md
verifBoolean = (abs( analy_wmax - numer_wmax   ) / abs(analy_wmax))  < 1e-3  ...
            && (abs( analy_maxMx - numer_maxMx ) / abs(analy_maxMx)) < 5e-3 ...
            && (abs( analy_dxmax - numer_dxmax ) / abs(analy_dxmax)) < 1e-3 ;
