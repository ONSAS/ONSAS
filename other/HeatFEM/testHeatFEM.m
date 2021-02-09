% =======================================
% sript for Heat transfer code validation
% =======================================


% =========   case 1 =========
% 1D problem
% boundary conditions: dirichlet on both ends
% initial temperature profile
% https://onsas.github.io/ONSAS_docs/dev/tutorials/HeatDiffusion/heat/

close all, clear all
addpath( genpath( '../../src/')); % add ONSAS src functions

timeIncr  = 0.0001 ;
Tfinal    = 0.02   ;
rho       = 1.     ;
cSpHe     = 1.     ;
kCond     = 4      ;
L         = 1      ;
Area      = 0.25   ; 
nelem     = 10     ;
Tdiri     = 0      ;

diriDofs = [ 1 nelem+1 ];
robiDofs = [] ;

nPlots = 4 ;

problemName = 'dirichlet' ;

initialTempFunc = 'myInitialTemp' ;


ndivs     = [ nelem  ];


boundaryCondParams = struct ( 'hConv'        , [], ...
                              'diriDofs'     , diriDofs, ...
                              'robiDofs'     , [], ...
                              'Tamb'         , [], ...
                              'qInpLeft'     , [], ...
                              'qInpRight'    , [], ...
                              'Tdiri'        , 0 );

[Ts, NodesCoord, times ] = HeatFEM( ...
  timeIncr, Tfinal, ...
  [rho, cSpHe, kCond], ...
  [ 1 L Area ], ...
  ndivs, ...
  boundaryCondParams, ...
  nPlots, problemName, initialTempFunc );


% plot history
figure, grid on, hold on
MS = 10 ; LW = 1.5 ;
indplot = round(nelem/2) ;

% numerical solution plot
plot( times, Ts(indplot,:), 'b',  'markersize', MS,'linewidth',LW )

% analytic solution computation
xsAnly = NodesCoord(indplot) ;
alpha = kCond / ( rho * cSpHe ) ;

tanali = exp(-(  pi)^2 * alpha * times ) * sin(pi * xsAnly) ...
                  + exp(-(3*pi)^2 * alpha * times ) * 0.5 * sin( 3 * pi * xsAnly ) ;

plot( times, tanali, 'r',  'markersize', MS,'linewidth',LW )

  xlabel('t'), ylabel('Temp')





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

