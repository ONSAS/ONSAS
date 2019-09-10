
inputONSASversion = '0.1.8';

problemName = 'ej2' ;

%~ Constitutive properties

E   = 30e9 ; % N/m2
nu  = 0.3   ;
rho = 0     ;

hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 E nu rho] ;

%~ Geometrical properties

%~ PNI 18
a  = .4  ;
A  = a^2*1e8 ; % m2
Iz = a^4/12 ;
Iy = a^4/12 ;
It = 1 ; 

secGeomProps = [ A Iy Iz It ] ;

l = 6 ; % m

Nodes = [ 0       0     0               ; ...
          0       0     l*.5               ; ...
          l     0     l*.5               ; ...
          l 0   0                        ] ;

aux = [ 1 2 1 ; ...
        2 3 1 ; ...
        3 4 1 ] ; 

Conec = [ aux(:,1:2) zeros(size(aux,1),2) ones(size(aux,1),1) aux(:,3) ones(size(aux,1),1)*2 ] ;
            
nodalSprings = [ 1 inf inf inf inf inf inf ; ...
                 4 inf inf inf inf inf inf ] ;

selfWeightBoolean = 0 ;

unifLoad = 		 [ 1 1   50e3 0 0 ; ...
                 3 1 -50e3 0 0 ] ;

%~ Analysis parameters                       
nonLinearAnalysisBoolean = 0 ;
plotParamsVector = [ 2 ] ; 
printflag = 2 ;                       
linearDeformedScaleFactor = 10 ;
