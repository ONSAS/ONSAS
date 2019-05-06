% --------------------------------------------------------------------------------------------------
% IETFEM problem 1D barra
% --------------------------------------------------------------------------------------------------

inputONSASversion = '0.1.8';
problemName = 'IETFEM_1D_barra' ;

% Propiedades mecanicas
E = 210000000 ; % kN/m2
nu = 0.3 ;

hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 E nu] ;

% Propiedades geométricas
A1 = 0.001 ;
Iy1 = 1 ;
Iz1 = 1 ;
It1 = 1 ;

A2 = 0.002 ;
Iy2 = 1 ;
Iz2 = 1 ;
It2 = 1 ;

secGeomProps = [ A1 Iy1 Iz1 It1 ; ...
								 A2 Iy2 Iz2 It2 ] ; 

% Nodos
Nodes = [ 0 		0 0 ; ...
					3.28 	0 0 ; ...
					6.56 	0 0 ] ;
					
% Conectividad
Conec = [ 1 2 0 0 1 1 1 ; ...
					2 3 0 0 1 2 1 ] ;					

% Apoyos de la estructura 
nodalSprings = [ 1 inf inf inf inf inf inf ] ;

% Fuerzas
Fx = -1 ;
Mx = 0 ;
Fy = 0 ; 
My = 0 ;
Fz = 0 ;
Mz = 0 ;

nodalConstantLoads = [ 3 Fx Mx Fy My Fz Mz ] ;

% Opciones de analisis
nonLinearAnalysisBoolean = 0 ; 
dynamicAnalysisBoolean   = 0 ; 
LBAAnalyFlag             = 0 ;
dim = 1 ;

% Opciones de plot 

printflag = 0 ;
plotParamsVector = [2 1 1] ;
linearDeformedScaleFactor = 50000 ;

