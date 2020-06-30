% --------------------------------------------------------------------------------------------------
% Test problema: rectangular crosssection cantilever beam submitted to imposed displacement at the free end
% --------------------------------------------------------------------------------------------------


% Parametros constitutivos del material

E  = 2.1e5 ; % modulo de young kN/m2
nu = 0.3   ; % coeficiente de Poisson

% Parametros seccionales

l = 1   ; % largo de las barras
b = 0.1 ; % ancho de la seccion
h  = b * 1.5 ; % fixed-ratio heigth given by local z

% Fuerzas externas

Fx = 0 ;
Mx = 0 ;
Fy = 0 ;
My = 0 ;
Fz = 0 ;
Mz = 0 ;

inputONSASversion = '0.1.9';

problemName = 'prescribedDisps' ;

% Constitutive properties

hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 E nu] ;

% Geometrical properties

A  = b * h ;  
Iy = b * h^3 / 12 ;
Iz = h * b^3 / 12 ;
It = 0.196 * h * b^3 ;

secGeomProps = [ A Iy Iz It ] ;

% Boundaries of the structure

nodalSprings = [ 1  inf  inf  inf  inf  inf  inf ] ;

% Imposed displacements

Fz = -1 ; 

prescribedDisps = [ 2 5 Fz*l^3/(3*E*Iy) ] ;

% Nodes coordinate matrix

Nodes = [ 0 0 0 ; ...
          l 0 0 ] ;
          

% Conectivity matrix

Conec = [ 1 2 0 0 1 1 2 ] ;

% Loads applied

nodalConstantLoads = [ 2  Fx Mx Fy My Fz*0 Mz ] ;

% [ node nodaldof scalefactor(positive or negative) ]

controlDofInfo = [ 2 5 -1 ] ;

% Analysis parameters

nonLinearAnalysisBoolean = 0 ; 
dynamicAnalysisBoolean   = 0 ; 
LBAAnalyFlag             = 0 ;


printflag = 2 ;
plotParamsVector = [2 ];
tablesBoolean = 1 ;

linearDeformedScaleFactor = 1;

% Analytical solution 

analyticSolFlag = 3 ;

analytSol = [ -Fz*l^2/(2*E*Iy) ] ;  
analyticSolDofs = 10 ;
analyticCheckTolerance = 1e-4 ;
% --------------------------------------------------------------------------------------------------
