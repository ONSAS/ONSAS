%% cylinder solid
%
%
%%

% uncomment to delete variables and close windows
clear all, close all

addpath('../sources/');

%% General data
dirOnsas = '..' ;
problemName = 'cylinder' ;

% Pressure
p = 1e3 ;

% Mesh data
[ nodesMat, conecMat ] = mshReader('cylinder.msh') ;

% Material properties
E = 210e9 ; nu = 0.3 ;
hyperElasParams = cell(1,1) ;  
hyperElasParams{1} = [ 1 E nu ] ;

% Sections
secGeomProps = [ 0 0 0 0 ] ;

Rint = 0.2 ; Rext = 0.24 ;

suppsMat = [ inf 0  0 	0   0 	0 ; ...
             0 	 0  inf 0   0   0 ; ...
             0 	 0  0   0   inf 0 ] ;

loadsMat = [0   0 0 0 0 -p 0 ] ;


[Nodes, Conec, nodalVariableLoads, nodalConstantLoads, unifDisLoadL, unifDisLoadG, nodalSprings ] = inputFormatConversion ( nodesMat, conecMat, loadsMat, suppsMat ) ;

% Analysis parameters
nonLinearAnalysisBoolean = 0 ; linearDeformedScaleFactor = 10.0 ;

% Analytic solution
analyticSolFlag = 3 ; p = abs(p) ; r = Rext ;
a = ( (1+nu)*(1-2*nu)*Rint^2*p ) / ( E*(Rext^2-Rint^2) ) ;
b = ( (1+nu)*Rint^2*Rext^2*p )   / ( E*(Rext^2-Rint^2) ) ;
analytSol = a*r + b/r ; analyticSolDofs = [ 6*(6-1)+1 ] ;
analyticCheckTolerance = 1e-3 ;

%% Output parameters
plotParamsVector = [ 3 ] ; printflag = 2 ;

%% ONSAS execution
% move to onsas directory and ONSAS execution

acdir = pwd ;
cd(dirOnsas);
ONSAS
cd(acdir) ;
  
