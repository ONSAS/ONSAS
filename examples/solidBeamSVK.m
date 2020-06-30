%% Example uniaxialSolid
% Linear elastic solid submitted to uniaxial loading. 
% Geometry given by $L_x$, $L_y$ and $L_z$, tension $p$ applied on 
% face $x=L_x$.
%
% Analytical solution to be compared with numerical:
% $$ u_x(x=L_x,y,z) = \frac{p L_x}{E} $$
%%

% uncomment to delete variables and close windows
clear all, close all

%% General data
dirOnsas = [ pwd '/..' ] ;
problemName = 'solidBeamSVK' ;

addpath( [ dirOnsas '/sources/' ] );



%% Structural properties

% an 8-node mesh is considered with its connectivity matrix
% Nodes and Conectivity matrix from .dxf file
[ nodesMat, conecMat, physicalNames ] = msh4Reader('solidBeamSVK.msh') ;
[ nodesMat, conecMat ] = esmacParser( nodesMat, conecMat, physicalNames ) ;

          
% Material and geometry properties
E = 200e9 ; nu = 0.3 ;
  
hyperElasParams = cell(1,1) ;  
hyperElasParams{1} = [ 6 E nu ] ;

secGeomProps = [ 0 0 0 0 ] ;

% Displacement boundary conditions and springs
suppsMat = [ inf 0  inf 	0  inf 	0 ] ;

loadsMat = [1   0 0 0 0 -1 0 ] ;

[Nodes, Conec, nodalVariableLoads, nodalConstantLoads, unifDisLoadL, unifDisLoadG, nodalSprings ] = inputFormatConversion ( nodesMat, conecMat, loadsMat, suppsMat ) ;


%% Analysis parameters
nonLinearAnalysisBoolean = 1 ; linearDeformedScaleFactor = 1.0 ;

stopTolIts       = 30     ;
stopTolDeltau    = 1.0e-8 ;
stopTolForces    = 1.0e-8 ;
targetLoadFactr  = 2e9    ;
nLoadSteps       = 20    ;

controlDofInfo = [ 7 5 -1 ] ;

numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactr nLoadSteps ] ; 


%% Output parameters
plotParamsVector = [ 3 ] ;
printflag = 2 ;

reportBoolean = 0;

%% ONSAS execution
% move to onsas directory and ONSAS execution
acdir = pwd ;
cd(dirOnsas);
ONSAS
cd(acdir) ;
  
