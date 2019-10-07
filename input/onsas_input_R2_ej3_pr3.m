% --------------------------------------------------------------------------------------------------
% Test problem: Practico 3 ej 3 R2.
% --------------------------------------------------------------------------------------------------

E  = 210e6   ;
nu = 0.3    ;
Iy = 1450 / (100^4) ;
A  = sqrt( Iy )*1e9 ; 
%~ A	 = 27.9 / (100^2) ;

qz = -15 ;
H = 15 ;

inputONSASversion = '0.1.8';

problemName = 'R2_pr3_ej3' ;

% constitutive properties

hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 E nu] ;

% geometrical properties

secGeomProps = [ A Iy 1 1 ] ;

Nodes = [ 0 0 0 ; ...
					0 0 2.4 ; ...
					0 0 4.8 ; ...
					4 0 4.8 ; ...
					5.2 0 4.8 ; ...
					4 0 0 ] ;

% in global system of coordinates
nodalSprings = [ 1  	inf  inf  inf  inf  	inf  inf ; ...
                 6  	inf  inf  inf  	 0  	inf  inf ] ;


%%%%%%%%%%%%%%%%%%% opcion con elemento de biela %%%%%%%%%%%%%%%%%%%%%%%
Conec = [ 1 2 0 0 1 1 2 ; ...
          2 3 0 0 1 1 2 ; ...
          3 4 0 0 1 1 2 ; ...
          4 5 0 0 1 1 2 ; ...
          4 6 0 0 1 1 2 ; ...
          1 4 0 0 1 1 1 ] ;

%%%%%%%%%%%%%%%%%%% opcion con elemento de biela %%%%%%%%%%%%%%%%%%%%%%%
%~ Conec = [ 1 2 0 0 1 1 2 ; ...
          %~ 2 3 0 0 1 1 2 ; ...
          %~ 3 4 0 0 1 1 2 ; ...
          %~ 4 5 0 0 1 1 2 ; ...
          %~ 4 6 0 0 1 1 2 ; ...
          %~ 1 4 0 0 1 1 2 ] ;

%~ Releases = [ 6 1 1 0 0 ] ;           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

unifLoad = [ 3 1 0 0 qz ; ...
             4 1 0 0 qz ] ;

nodalConstantLoads = [ 2 H 0 0 0 0 0 ] ; 

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

linearDeformedScaleFactor = 1 ;  

% --------------------------------------------------------------------------------------------------
