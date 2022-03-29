%md# Lamb's wave problem
%md
%md In this case the ONSAS.m directory is loaded from an environment variable.
% add path
clear all, close all
addpath( genpath( [ getenv( 'ONSAS_PATH' ) 'src/'] ) )
% scalar parameters
E = 1.8773e10 ; nu = 0.25 ; thickness = 1 ; rho = 2200 ;
%md
global spitMatrices, spitMatrices = true ;
%md
%md## MEBI parameters
%md
%md### materials
materials.hyperElasModel  = 'linearElastic' ;
materials.hyperElasParams =  [ E nu ]       ;
materials.density         =  rho            ;
%md
%md### elements
elements(1).elemType = 'node';
%
elements(2).elemType = 'edge';
elements(2).elemCrossSecParams = thickness ;
%
elements(3).elemType = 'triangle';
elements(3).elemTypeParams = 2 ;
elements(3).elemCrossSecParams = thickness ;
%md
%md#### boundaryConds
boundaryConds(1).imposDispDofs = [1 3] ;
boundaryConds(1).imposDispVals = [0 0] ;
%
boundaryConds(2).imposDispDofs = [1] ;
boundaryConds(2).imposDispVals = [0] ;
%
boundaryConds(3).loadsCoordSys = 'global' ;
boundaryConds(3).loadsTimeFact = @(t) 1e6*( 1*( t <= .15 ) - 3*( t <= .1 ) + 3*( t <= .05 ) ) ;
boundaryConds(3).loadsBaseVals = [ 0 0  1 0  0 0 ]  ;
%md
%md#### initialConds
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
analysisSettings.deltaT        = dt      ;
%md
%md
%md### Output parameters
otherParams.problemName = 'lambsProblem' ;
otherParams.plotsFormat = 'vtk' ;
%md
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md
