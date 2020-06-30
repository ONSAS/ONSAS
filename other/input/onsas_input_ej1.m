
inputONSASversion = '0.1.9';

problemName = 'ej1' ;

%~ Constitutive properties

E   = 210e9 ; % N/m2
nu  = 0.3   ;
rho = 0     ;

hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 E nu rho] ;

%~ Geometrical properties

%~ PNI 18
fie = .15  ;
fii = .015 ;
A18 = 1 ; % m2
Iz18 = pi * ( fie^4 - fii^4) / 64 ; % m4 
Iy18 = pi * ( fie^4 - fii^4) / 64 ; % m4 
It18 = Iy18 * 2 ; 


secGeomProps = [ A18 Iy18 Iz18 It18 ] ;

%~ Nodes & Conectivity & Releases & Supports

l = 1 ; % m

Nodes = [ 0       0     0               ; ...
          l       0     0               ; ...
          2*l     0     0               ; ...
          l 0 -l                        ] ;

aux = [ 1 2 1 ; ...
        2 3 1 ; ...
        2 4 1 ] ; 

Conec = [ aux(:,1:2) zeros(size(aux,1),2) ones(size(aux,1),1) aux(:,3) ones(size(aux,1),1)*2 ] ;
            
nodalSprings = [ 1 inf inf inf inf inf inf ; ...
                 3 inf inf inf inf inf inf ; ...
                 4 inf inf inf inf inf inf ] ;

nodalConstantLoads   = [ 2  0 0    -200000 0 0 0 ] ;

selfWeightBoolean = 0 ;

%~ Analysis parameters                       
nonLinearAnalysisBoolean = 0 ;
plotParamsVector = [ 2 ] ; 
printflag = 0 ;                       
linearDeformedScaleFactor = 100 ;



