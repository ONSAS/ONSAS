% ------------------------------------------------------------------------------ 
% ------      ONSAS example file: cantilever with nodal moment example    ------
% ------------------------------------------------------------------------------
clear all, close all

problemName = 'emparrillado_Mitad_Lineal' ;

% ------------------------------------------------------------------------------ 
% --- scalar parameters
l = 2 ;   width = .05 ;   P = 4e3 ;

Nodes = [ 0  0  l/2 ; ...
          l  0  l/2 ; ...
          l  0  0   ] ;

% ------------------------------------------------------------------------------ 
% --- MELCS parameters ---

E = 210e9 ;  nu = 0.3 ;  rho = 0 ;
materialsParams = { [ rho 1 E nu] } ;

% 1 node:  3: beam
elementsParams  = { 1; 3} ;

loadsParams   = {[ 1 1  0 0 -P/2 0 0 0 ]} ;

% --- cross section ---
A = width^2 ;   Iy = width^4/12 ;  Iz = Iy ;  It = 0.141 * width^4 ;

crossSecsParams = {[ 1 A It Iy Iz ]}      ;

springsParams      = {} ;
springsParams{1,1} = [ inf  inf  inf  inf  inf  inf ] ;
springsParams{2,1} = [ 0    inf  0    0    inf  0   ] ;

% ------------------------------------------------------------------------------ 
% --- connectivity ---

Conec = {} ;
Conec{1, 1} = [ 0 1 0 0 1   1   ] ; % fixed node
Conec{2, 1} = [ 0 1 1 0 2   3   ] ; % fixed/loaded node
Conec{3, 1} = [ 1 2 0 1 0   1 2 ] ; % beam
Conec{4, 1} = [ 1 2 0 1 0   2 3 ] ; % beam

plotParamsVector = [ 3 ] ; printFlag = 2 ;
storeBoolean     = 1 ;

% --- ONSAS execution ---
run( [ pwd '/../ONSAS.m' ] ) ;

matUsLineal = matUs

problemName = 'emparrillado_Mitad_NoLineal' ;

Conec = {} ;
Conec{1, 1} = [ 0 1 0 0 1   1   ] ; % fixed node
Conec{2, 1} = [ 0 1 1 0 2   3   ] ; % fixed/loaded node
Conec{3, 1} = [ 1 2 0 1 0   1 2 ] ; % beam
Conec{4, 1} = [ 1 2 0 1 0   2 3 ] ; % beam

stopTolIts      = 30     ;
stopTolDeltau   = 1.0e-10 ;
stopTolForces   = 1.0e-10 ;
targetLoadFactr = 20  ;
nLoadSteps      = 20 ;

plotParamsVector = [ 3 ] ;

numericalMethodParams = [ 1 ...
 stopTolDeltau stopTolForces stopTolIts ...
 targetLoadFactr nLoadSteps  ] ; 


% --- ONSAS execution ---
run( [ pwd '/../ONSAS.m' ] ) ;

problemName = 'emparrillado_Entero_NoLineal' ;

loadsParams   = {[ 1 1  0 0 -P 0 0 0 ]} ;

springsParams      = {} ;
springsParams{1,1} = [ inf  inf  inf  inf  inf  inf ] ;

Nodes = [ 0  0   l/2 ; ...
          l  0   l/2 ; ...
          l  0   0   ; ...
          l  0  -l/2 ; ...
          0  0  -l/2 ] ;
clear matNs

Conec = {} ;
Conec{1, 1} = [ 0 1 0 0 1   1   ] ; % fixed node
Conec{2, 1} = [ 0 1 1 0 0   3   ] ; % fixed/loaded node
Conec{3, 1} = [ 0 1 0 0 1   5   ] ; % fixed/loaded node
Conec{4, 1} = [ 1 2 0 1 0   1 2 ] ; % beam
Conec{5, 1} = [ 1 2 0 1 0   2 3 ] ; % beam
Conec{6, 1} = [ 1 2 0 1 0   3 4 ] ; % beam
Conec{7, 1} = [ 1 2 0 1 0   4 5 ] ; % beam

% --- ONSAS execution ---
run( [ pwd '/../ONSAS.m' ] ) ;
