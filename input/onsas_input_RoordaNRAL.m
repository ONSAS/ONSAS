% ----------------------------
% --- Roorda Frame Example ---
% ---    J. B. Bazzano     ---
% ----------------------------

% datos del ejemplo: en mm, N, y MPa.
% acero: E=200E3 MPa, nu=0.3.
% portico dimensiones: 1000mm x 1000mm
% apoyos simples. Estructura restringida al plano.
% Seccion barras: ancho=50mm, alto=10mm.
% Numero de elementos por barra: m

E  = 200e3 ;   nu = 0.3   ;
L= 1000   ;
a = 10    ;
b = 50    ;
m = 4     ;

inputONSASversion = '0.1.8';
problemName = 'RoordaNRAL' ;

hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 E nu] ;

sectPar = [ b a ]; 

A  = b*a      ; It = a*b^3/3 ;
Iy = b*a^3/12 ; Iz = a*b^3/12 ;

secGeomProps = [ A Iz Iy It ] ; % viga vertical [Iz Iy]

Nelem = 2*m ;  Np = Nelem +1;

dl = L/m ;
Nodes = zeros(Np,3);

% barra vertical
for i=1:(m+1)
  Nodes(i,1) = 0 ;            %x
  Nodes(i,2) = 0 + (i-1)*dl;  %y
  Nodes(i,3) = 0 ;            %z
end

% barra horizontal
for i=m+2:Np
  Nodes(i,1) = 0 + (i-(m+1))*dl;  %x
  Nodes(i,2) = 0 + L;         %y
  Nodes(i,3) = 0 ;            %z
end
                %node      ux   tx   uy   ty   z    tz
nodalSprings = [ 1         inf  inf  inf  inf  inf  0 ;
                 m+1        0   inf  0    inf  inf  0 ;
                 Np        inf  inf  inf  inf  inf  0 ] ;

Conec = [ (1:(Np-1))' (2:(Np))' zeros(Nelem,2)  (ones(Nelem,1)*[ 1 1 2]) ] ;

                      %nod    fx  mx  fy  my  fz  mz
nodalVariableLoads   = [ m+1    0   0   -1   0   0   .2 ] ; %mz>0 gira antih, imperf. sensitive

controlDofInfo = [ m+1 3 -1 ] ;

stopTolIts     = 30     ;
stopTolDeltau  = 1.0e-8 ;
stopTolForces  = 1.0e-8 ;
targetLoadFactr = 12.8e3 ;
nLoadSteps      = 100 ; incremArcLen     = 1    ;

plotParamsVector = [ 2 4 ];  plotsViewAxis = [ 0 0 1];  

printflag = 1;   reportBoolean = 1 ;

numericalMethodParams = [ 2 stopTolDeltau stopTolForces stopTolIts targetLoadFactr nLoadSteps  incremArcLen ] ;  % ArcLength
