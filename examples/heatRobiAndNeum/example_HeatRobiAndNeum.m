

close all, clear all
addpath( genpath( '../../src/')); % add ONSAS src functions


timeIncr  = 0.001;
rho       = 1.     ;   cSpHe     = 1.     ;
kCond     = 4      ;   L         = 1      ;
Area      = 0.25   ;   nelem     = 20     ;

Tdiri     = 1  ; Tfinal = .1 ;
ndivs = [ nelem ] ;
diriDofs = [  ] ;
robiDofs = [ 1 ]
qInpRight = .5 ;

problemName = 'HeatRobiAndNeum' ;

ndivs     = [ nelem  ];
hConv = 10 ; Tamb = [];
qInpLeft = [];
Tdiri = 0 ;
nPlots = 4 ;

ambTempFuncName = 'myAmbTempFunc' ;
initFuncName =  'myInitialTemp' ;

[Ts, NodesCoord, times ] = HeatFEM( ...
  timeIncr, Tfinal, ...
  [rho, cSpHe, kCond], ...
  [ 1 L Area ], ...
  ndivs, ...
  hConv, diriDofs, robiDofs, ambTempFuncName, qInpLeft, qInpRight, Tdiri, ...
  nPlots, problemName, initFuncName, [], [], [], [] );





Lx        = 1      ; Ly        = .5     ;
Lz        = .5     ;
assert( Area == (Ly*Lz) ) ;

nPlots = inf ;

initialTempFunc = 'myInitialTemp' ;

ndivs     = [ nelem 2 2 ];

qInpLeft = []; qInpRight = [];

qInp = 0;

diriFaces = [  ];
neumFaces = [ qInp  2 3 4 5 6 ] ;
robiFaces = [ hConv 1 ] ;

[Ts3D, NodesCoord, times3D ] = HeatFEM( ...
  timeIncr, Tfinal, ...
  [rho, cSpHe, kCond], ...
  [ 2 Lx Ly Lz ], ...
  ndivs, ...
  hConv, diriDofs, robiDofs, ambTempFuncName, qInpLeft, qInpRight, Tdiri, ...
  nPlots, problemName, initFuncName, [], ...
  diriFaces, neumFaces, robiFaces  );





indplot = round( nelem/2 )+1;
xs = linspace( 0, Lx,nelem+1)' ;

% numerical solution plot
MS = 18 ; LW = 1.5 ;

figure, hold on, grid on
plot( times, Ts(indplot,:), 'b',  'markersize', MS,'linewidth',LW )
plot( times, Ts3D(indplot,:), 'g',  'markersize', MS,'linewidth',LW )
legend('1D','3D')

figure, hold on, grid on
plot( xs, Ts(:,50), 'b-o',  'markersize', MS,'linewidth',LW )
plot( xs, Ts3D(1:(nelem+1),50), 'g-s',  'markersize', MS,'linewidth',LW )
plot( xs, Ts3D((1:(nelem+1))+(nelem+1),50), 'r-x',  'markersize', MS,'linewidth',LW )
legend('1D','3D-line1','3D-line2')


