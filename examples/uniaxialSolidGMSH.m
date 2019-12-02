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
dirOnsas = '..' ;
problemName = 'uniaxialSolidManual' ;

%% Structural properties

% Nodes and Conectivity matrix from .dxf file
[ nodesMat, conecMat ] = mshReader('uniaxialSolid.msh') ;

% Support matrix	: Is defined by the corresponding support label. I.e., in torre.dxf there is ony one label for supports, then
% 									the matrix will have only one row. The structure of the matrix is: [ ux thetax uy thetay uz thetaz ]
suppsMat = [ inf 0  0 	0   0 	0 ; ...
             0 	 0  inf 0   0   0 ; ...
             0 	 0  0   0   inf 0 ] ;

% Loads matrix: 		Is defined by the corresponding load label. First entry is a boolean to assign load in Global or Local axis. (Recommendation: Global axis). 
%										Global axis -> 1, local axis -> 0. 
%										The structure of the matrix is: [ 1/0 Fx Mx Fy My Fz Mz ]
p = -210e8 ; Lx = 0.5 ; Ly = 0.5 ; Lz = 0.5 ;

loadsMat = [0   0 0 0 0 p 0 ] ;

% Previously defined matrices to ONSAS format
[Nodes, Conec, nodalVariableLoads, nodalConstantLoads, unifDisLoadL, unifDisLoadG, nodalSprings ] = inputFormatConversion ( nodesMat, conecMat, loadsMat, suppsMat ) ;


% Constitutive properties: Structure: [ 1 E nu ]
E = 210e9 ; nu = 0.3 ;
hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [ 1    210e9    0.3 ] ;

% Sections
secGeomProps = [ 0 0 0 0 ] ;

% Analysis parameters
nonLinearAnalysisBoolean 	= 0 ;
dynamicAnalysisBoolean 		= 0 ;


% Analytic sol
analyticSolFlag = 3 ; analytSol = [ p*Lx/E ] ; analyticSolDofs = [ 6*(7-1)+1 ] ;
analyticCheckTolerance = 1e-8 ;

%% Output parameters
plotParamsVector = [ 3 ] ; printflag = 2 ;

%% ONSAS execution
% move to onsas directory and ONSAS execution

acdir = pwd ;
cd(dirOnsas);
ONSAS
cd(acdir) ;
  
