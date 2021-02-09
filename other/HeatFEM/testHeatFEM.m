
% =========   caso de validacion 1 =========
% condiciones de contorno dirichlet con perfil de temperaturas inicial de https://onsas.github.io/ONSAS_docs/dev/tutorials/HeatDiffusion/heat/

close all, clear all
addpath( genpath( '../../src/'));

timeIncr  = 0.0001 ;
Tfinal    = 0.02 ;
rho       = 1. ;
cSpHe     = 1. ;
kCond     = 4 ;
Lx         = 1  ;
Ly         = .5  ;
Lz         = .5  ;
ndivs = [ 2 1 1  ];

hConv     = 10 ;

Tamb      = 20 ;
Tdiri     = 0  ;

%~ qInpLeft  = 0  ;
%~ qInpRight = 0  ;

initialTempFlag = 1 ;
anlyBoolean = 1 ;

%~ diriDofs = [ 1 nelem+1 ];
%~ robiDofs = [] ;

plotBoolean = true ;
nCurves = 4 ;

problemName = 'dirichlet' ;
  
HeatFEM( ...
  timeIncr, Tfinal, ...
  [rho, cSpHe, kCond], ...
  [ 2 Lx Ly Lz], ...
  ndivs, ...
  plotBoolean, nCurves, anlyBoolean, problemName );
return

Tdiri     = 1  ;
anlyBoolean = 0 ;
Tfinal = .5 ;

onsasExample_1DHeatTransfer( ...
  timeIncr, Tfinal, ...
  rho, cSpHe, kCond, ...
  L, Area, hConv, ...
  nelem, Tamb, Tdiri, ...
  initialTempFlag, ...
  diriDofs, robiDofs, qInpLeft, qInpRight, ...
  plotBoolean, nCurves, anlyBoolean, problemName );


% =========   caso de validacion 2 =========
%  diri no homogeneo ( ver fundamentos en https://www.math.uzh.ch/li/index.php?file&key1=25297 )

nCurves = 5 ;

qInpLeft  = 0  ;
qInpRight = 4  ;

Tdiri     = .5  ;
diriDofs = [ 1 ];

anlyBoolean = 0 ;

problemName = 'nonHomDiriAndNeumann' ;

onsasExample_1DHeatTransfer( ...
  timeIncr, Tfinal, ...
  rho, cSpHe, kCond, ...
  L, Area, hConv, ...
  nelem, Tamb, Tdiri, ...
  initialTempFlag, ...
  diriDofs, robiDofs, qInpLeft, qInpRight, ...
  plotBoolean, nCurves, anlyBoolean, problemName );



% =========   caso de validacion 3 =========



nCurves = 5 ;
Tfinal = 0.5;

qInpRight = 4  ;

Tdiri     = .5  ;
Tamb      = 5 ;
robiDofs  = [ 1 ] ;
anlyBoolean = 0 ;

problemName = 'neumAndNeumann' ;

onsasExample_1DHeatTransfer( ...
  timeIncr, Tfinal, ...
  rho, cSpHe, kCond, ...
  L, Area, hConv, ...
  nelem, Tamb, Tdiri, ...
  initialTempFlag, ...
  diriDofs, robiDofs, qInpLeft, qInpRight, ...
  plotBoolean, nCurves, anlyBoolean, problemName );

