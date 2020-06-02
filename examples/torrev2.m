% ONSAS version
inputONSASversion = '0.1.9';

% Problem name in order to write LaTex report 
problemName = 'torre' ; 

% Nodes and Conectivity matrix from .dxf file
[ nodesMat, conecMat ] = dxf2ONSAS('torre.dxf') ;


% Support matrix	: Is defined by the corresponding support label. I.e., in torre.dxf there is ony one label for supports, then
% 									the matrix will have only one row. The structure of the matrix is: [ ux thetax uy thetay uz thetaz ]
suppsMat = [ inf 0 inf 0 inf 0 ] ;

% Loads matrix: 		Is defined by the corresponding load label. First entry is a boolean to assign load in Global or Local axis. (Recommendation: Global axis). 
%										Global axis -> 1, local axis -> 0. 
%										The structure of the matrix is: [ 1/0 Fx Mx Fy My Fz Mz ]
loadsMat = [1 0 0 -4e5 0 -2e5 0 ] ;

% Previously defined matrices to ONSAS format
[Nodes, Conec, nodalVariableLoads, nodalConstantLoads, unifDisLoadL, unifDisLoadG, nodalSprings ] = inputFormatConversion ( nodesMat, conecMat, loadsMat, suppsMat )


% Constitutive properties: Structure: [ 1 E nu ]
hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [ 1    210e9    0.3 ] ;

% Geometrical properties: Structure: [ A Iy Iz J ]
secGeomProps = [    2.000e-03    0.000e+00    0.000e+00    0.000e+00 ] ;


% Analysis parameters
nonLinearAnalysisBoolean 	= 0 ;
dynamicAnalysisBoolean 		= 0 ;


% Output parameters
plotParamsVector 	= [2   2] ; 
reportBoolean 		= 1 ; 
printflag 				= 0 ; 
linearDeformedScaleFactor = 70 ; 
plotsViewAxis 						= [ -1.5 1 1 ] ; 
