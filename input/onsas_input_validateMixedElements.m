% --------------------------------------------------------------------------------------------------
% Test problem: rectangular cross-section cantilever beam submitted to nodal external loads at the free end.
% --------------------------------------------------------------------------------------------------

E  = 210e9   ;
nu = 0.3     ;
l  = 3       ; % reference length
Iy = 1943e-8 ;
Fz = -1e6 ;
A  = 28.5e-4 ; 

inputONSASversion = '0.1.9';

problemName = 'validateMixedElements' ;

% constitutive properties

hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 E nu] ;

% geometrical properties

secGeomProps = [ A*1e4 Iy 1 1 ; ...
                 A*1e4 Iy*1e3 1 1 ] ;

Nodes = [   0  0  0 ; ...
            l  0  0 ; ... 
          2*l  0  0 ; ...
            l  0  l ] ;

% in global system of coordinates
nodalSprings = [ 1  inf  inf  inf    0  inf  inf ; ...
                 2    0  inf  inf    0    0  inf ; ...
                 3    0  inf  inf    0    0  inf ; ...
                 4  inf  inf  inf    0  inf  inf ] ;

Conec = [ 1 2 0 0 1 1 2 ; ...
          2 3 0 0 1 1 2 ; ...
          2 4 0 0 1 2 2 ] ;

nodalConstantLoads   = [ 3  0 0 0 0  Fz 0 ] ;

Releases = [3 1 0 1 0];
%~ Releases = [3 0 1 0 1];
%~ Releases = [3 1 1 1 1];

% analysis parameters
nonLinearAnalysisBoolean = 0 ; 
dynamicAnalysisBoolean   = 0 ; 
LBAAnalyFlag             = 0 ;

% [ node nodaldof scalefactor(positive or negative) ]
controlDofInfo = [ 3 5 -1 ] ;

printflag = 0 ;
tablesBoolean = 1;

plotParamsVector = [2 2];

plotsViewAxis = [ 0 -1 0] ;

linearDeformedScaleFactor = 1;

% analytical solution 
analyticSolFlag = 0 ;
G  = E / ( 2*(1+nu) )               ;
analytSol = [ -3.43990/(E*Iy)*Fz ]' ;
analyticSolDofs = [ 7 ]' ;
  %~ -1.90385
   %~ 0.93930
analyticCheckTolerance = 1e-4 ;  

% --------------------------------------------------------------------------------------------------
