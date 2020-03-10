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
problemName = 'uniaxialSVKSolidManual' ;

addpath( [ dirOnsas '/sources/' ] );

%% Structural properties

% tension applied and x, y, z dimensions
p = 1 ; Lx = 1 ; Ly = 1 ; Lz = 1 ;


[ nodesMat, conecMat, physicalNames ] = msh4Reader('uniaxialSolid.msh') ;
[ nodesMat, conecMat ] = esmacParser( nodesMat, conecMat, physicalNames ) ;

suppsMat = [ inf 0  0 	0   0 	0 ; ...
             0 	 0  inf 0   0   0 ; ...
             0 	 0  0   0   inf 0 ] ;

% Loads matrix: 		Is defined by the corresponding load label. First entry is a boolean to assign load in Global or Local axis. (Recommendation: Global axis). 
%										Global axis -> 1, local axis -> 0. 
%										The structure of the matrix is: [ 1/0 Fx Mx Fy My Fz Mz ]

%~ loadsMat = [0   0 0 0 0 p 0 ] ;  % --- local loading ---
loadsMat = [1   p 0 0 0 0 0 ] ; % --- global loading ---

[Nodes, Conec, nodalVariableLoads, nodalConstantLoads, unifDisLoadL, unifDisLoadG, nodalSprings ] = inputFormatConversion ( nodesMat, conecMat, loadsMat, suppsMat ) ;

clear nodesMat conecMat loadsMat suppsMat
          
% Material and geometry properties
E = 1 ; nu = 0.3 ;
  
hyperElasParams = cell(1,1) ;  
hyperElasParams{1} = [ 6 E nu ] ;

secGeomProps = [ 0 0 0 0 ] ;


%% Loading parameters
%~ nodalForce = p * Ly * Lz / 6 ;
%~ nodalVariableLoads = [ (5:8)' nodalForce*[1 2 1 2]' zeros(4,5) ] ;

%% Analysis parameters
nonLinearAnalysisBoolean = 1 ; linearDeformedScaleFactor = 1.0 ;

stopTolIts       = 30     ;
stopTolDeltau    = 1.0e-12 ;
stopTolForces    = 1.0e-12 ;
targetLoadFactr  = 2    ;
nLoadSteps       = 3    ;

controlDofInfo = [ 7 1 1 ] ;

numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactr nLoadSteps ] ; 

% Analytic sol
analyticSolFlag = 2 ;
analyticCheckTolerance = 1e-8 ;
analyticFunc = @(w) E * 0.5 * ( (1 + w/Lx).^3 - (1+w/Lx) )

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
  
