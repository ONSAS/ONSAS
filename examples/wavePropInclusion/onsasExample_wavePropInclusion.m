% Example meshing InclusionWave using ONSAS

% command to generate msh file from Octave:
%   system('gmsh -2 inclusionCirc.geo')

clear all, close all

addpath(genpath('../../') )

E = 10e6     % 40 MPa
nu = 0.49 ;  %
p = 1.0 ;    %
thickness = 1 ;
density   = 2200 ;

materials.hyperElasModel  = {'linearElastic'; 'linearElastic'} ;
materials.hyperElasParams = { [ E nu ]; [ 10*E nu ] }      ;
materials.density         = { density ; density }      ;

elements.elemType = { 'node', 'edge', 'triangle' } ;
elements.elemTypeParams = { []; [] ; 2  } ;
elements.elemTypeGeometry = { []; thickness ; thickness } ;

boundaryConds.loadsCoordSys = {[]; []; 'global'  } ;
boundaryConds.loadsTimeFact = { []; []; @(t) (t>0)*(t<1e-3) } ;
boundaryConds.loadsBaseVals = { []; []; [ -p 0 0  0  0 0 ]  } ;
boundaryConds.imposDispDofs = { [1 3] ; [3] ; []  } ;
boundaryConds.imposDispVals = { [0 0] ; [0] ; []  } ;
%md

initialConds = struct();
%md

[ mesh.nodesCoords, mesh.conecCell ] = meshFileReader( 'wavePropInclusion.msh' ) ;

analysisSettings.methodName    = 'newmark' ;
analysisSettings.stopTolIts    = 30      ;
analysisSettings.stopTolDeltau = 1.0e-12 ;
analysisSettings.stopTolForces = 1.0e-12 ;
analysisSettings.finalTime      = 1e-3       ;
analysisSettings.deltaNM       = 0.5    ;
analysisSettings.alphaNM       = 0.25      ;
analysisSettings.deltaT        = 1e-4    ;
%md
%md
%md### Output parameters
otherParams.problemName = 'wavePropInclusion' ;
otherParams.plotsFormat = 'vtk' ;
%otherParams.spitMatrices = true ;
%md
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
