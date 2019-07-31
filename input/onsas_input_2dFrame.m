
E  = 30e9   ;
nu = 0.3     ;
Iy = 1943e-8 ;
Fz = -1350e3 ;
A  = 28.5e-4 ; 

inputONSASversion = '0.1.8';

problemName = '2dFrame' ;

% constitutive properties

hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 E nu] ;

% geometrical properties

secGeomProps = [ A Iy 1 1 ] ;

l1 = 3;
l2 = 6;
l3 = 3;

Nodes = [   0  0   0 ; ...
            0  0  l1 ; ... 
           l2  0  l1 ; ...
           l2  0  -l3+l1 ] ;

% in global system of coordinates
nodalSprings = [ 1  inf  inf  inf  inf  inf  inf ; ...
                 4  inf  inf  inf  inf  inf  inf ; ...
                 2    0  inf  inf    0    0  inf ; ...
                 3    0  inf  inf    0    0  inf ] ;

Conec = [ 1 2 0 0 1 1 2 ; ...
          2 3 0 0 1 1 2 ; ...
          3 4 0 0 1 1 2 ] ;

nodalConstantLoads   = [ 3  0 0 0 1 0 0 ] ;


% analysis parameters
nonLinearAnalysisBoolean = 0 ; 
dynamicAnalysisBoolean   = 0 ; 
LBAAnalyFlag             = 0 ;

% [ node nodaldof scalefactor(positive or negative) ]
controlDofInfo = [ 1 1 1 ] ;


printflag = 0 ;
tablesBoolean = 1;

plotParamsVector = [2];

plotsViewAxis = [ 0 -1 0] ;

linearDeformedScaleFactor = 1e5;

% analytical solution 
analyticSolFlag = 0 ;
