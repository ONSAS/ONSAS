%md# Lamb's wave problem
%md
%md In this case the ONSAS.m directory is loaded from an environment variable.
% add path
clear all, close all
addpath( genpath( getenv( 'ONSAS_PATH' ) ) )
% scalar parameters
E = 1.8773e10 ; nu = 0.25 ; thickness = 1 ; rho = 2200 ;
%md
global spitMatrices
spitMatrices = true
%md
%md## MEBI parameters
%md
%md### materials
%md since only one material is considered, the structs defined for the materials contain only one entr
materials.hyperElasModel  = {'linearElastic'} ;
materials.hyperElasParams = { [ E nu ] }      ;
materials.density         = { rho }           ;
%md
%md### elements
%md In this model two kinds of elements are used: tetrahedrons for the solid and triangles for introducing the external loads. Since two kinds of elements are used, the structs have length 2:
elements.elemType = { 'node', 'edge', 'triangle' } ;
%md since triangle and tetrahedron elements dont have specific parameters the struct entries contain empty vectors
elements.elemTypeParams = { []; [] ; 2  } ;
elements.elemTypeGeometry = { []; thickness ; thickness } ;
%md
%md#### boundaryConds
%md in this case four BCs are considered, one corresponding to a load and three to displacements.
boundaryConds.loadsCoordSys = {[]     ; []  ; 'global'  } ;
boundaryConds.loadsTimeFact = { []    ; []  ; @(t) 1e6*( 1*( t <= .15 ) - 3*( t <= .1 ) + 3*( t <= .05 ) ) } ;
boundaryConds.loadsBaseVals = { []    ; []  ; [ 0 0  1 0  0 0 ]  } ;
boundaryConds.imposDispDofs = { [1 3] ; [1] ; []  } ;
boundaryConds.imposDispVals = { [0 0] ; [0] ; []  } ;
%md
%md#### initialConds
%md since no initial non-homogeneous initial conditions are used, an empty struct is used .
initialConds = struct();
%md

%md### Mesh
%md An 8-node mesh is considered with its connectivity matrix
%md
%md```@raw html
%md<img src="https://raw.githubusercontent.com/ONSAS/ONSAS_docs/master/docs/src/solidCubeMeshHTML.svg" alt="structure diagram" width="500"/>
%md```
%md```@raw latex
%md\begin{center}
%md\def\svgwidth{0.6\textwidth}
%md\input{solidCubeMeshPDF.pdf_tex}
%md\end{center}
%md```
%md
[ mesh.nodesCoords, mesh.conecCell ] = meshFileReader( 'geometry_lambsProblem.msh' ) ;
%md
%md### Analysis parameters
%md
nAproxElem = length( mesh.conecCell )

CFL = 0.125 ;
cP  = 3200  ;
h   = sqrt( 2*3200^2 / nAproxElem ) * .5 ;
dt  = CFL * h / cP ;
%
analysisSettings.methodName    = 'newmark' ;
analysisSettings.stopTolIts    = 30      ;
analysisSettings.stopTolDeltau = 1.0e-12 ;
analysisSettings.stopTolForces = 1.0e-12 ;
analysisSettings.finalTime      = 1    ;
analysisSettings.alphaNM       = 0.25    ;
analysisSettings.deltaNM       = 0.5     ;
analysisSettings.deltaT        = dt      ;
%md
%md
%md### Output parameters
otherParams.problemName = 'lambsProblem' ;
otherParams.plotsFormat = 'vtk' ;
%md
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md