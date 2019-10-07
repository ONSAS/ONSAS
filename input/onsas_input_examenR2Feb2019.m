%~ Examen R2 Febrero 2019 - 14/02/2019

inputONSASversion = '0.1.9';

problemName = 'examenR2Feb2019' ;

%~ Constitutive properties

E = 210e9 ; % N/m2
nu = 0.3 ;
rho = 7850 * 9.81 ; % N/m3

hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 E nu rho] ;

%~ Geometrical properties

%~ PNI 18
A18 =  27.9 / 100^2 ; % m2
Iz18 = 1450 / 100^4 ; % m4 
Iy18 = 81.3 / 100^4 ; % m4
It18 = 1 ; 

%~ PNI 20
A20 =  33.5 / 100^2 ; % m2
Iz20 = 2140 / 100^4 ; % m4 
Iy20 =  117 / 100^4 ; % m4
It20 = 1 ; 

secGeomProps = [ A18 Iy18 Iz18 It18 ; ...
                 A20 Iy20 Iz20 It20 ] ;

%~ Nodes & Conectivity & Releases & Supports

l20_1 = 5.5 ; % m
l20_2 = 0.9 ; % m
l18 = 8 ; % m

Nodes = [ 0       0     0               ; ...
          0       0   l20_1            ; ...
          0       0   (l20_1 + l20_2)  ; ...
          l18/2   0   (l20_1 + l20_2)  ; ...
          l18     0   (l20_1 + l20_2)  ; ...
          l18/2   0   l20_1            ; ...
          l18     0   l20_1            ; ...
          l18     0     0               ] ;

aux = [ 1 2 2 ; ...
        2 3 2 ; ...
        3 4 1 ; ...
        4 5 1 ; ...
        2 6 1 ; ...
        6 7 1 ; ...
        5 7 2 ; ...
        7 8 2 ] ; 

Conec = [aux(:,1:2) zeros(size(aux,1),2) ones(size(aux,1),1) aux(:,3) ones(size(aux,1),1)*2 ] ;

Releases = [ 3 1 0 1 0 ; ...
             4 0 1 0 1 ; ...
             5 1 0 1 0 ; ...
             6 0 1 0 1 ] ;
             
nodalSprings = [ 1 inf inf inf inf inf inf ; ...
                 8 inf inf inf inf inf inf ] ;

%~ Loads

selfWeightBoolean = 1 ;


% Wind
a = 1 ; % m
b = 8 ; % m
w = 1000 ; % N/m2
q_wind_y = w * a / 2 ; % N/m

% pp

Pp = 800 ; % N
q_peso_z = - Pp / ( 2 * ( l18 ) ) ; % N/m
                   
unifLoad = 		 [ 3 0 0 q_wind_y q_peso_z ; ...
                 4 0 0 q_wind_y q_peso_z ; ...                   
                 5 0 0 q_wind_y q_peso_z ; ...                   
                 6 0 0 q_wind_y q_peso_z ] ;
                                   
                  
%~ Analysis parameters                       
nonLinearAnalysisBoolean = 0 ;
plotParamsVector = [ 2 ] ; 
printflag = 0 ;                       
linearDeformedScaleFactor = 10 ;



