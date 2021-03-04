

close all, clear all
addpath( genpath( '../../src/')); % add ONSAS src functions

timeIncr  = 0.001  ;   Tfinal    = .1    ;
rho       = 1.     ;   cSpHe     = 1.    ;
kCond     = 4      ;   L         = 1     ;
Area      = 0.25   ;   nelem     = 20    ;

ndivs = [ nelem ] ;

diriDofs = [  ] ; robiDofs = [ 1 ]
qInpLeft = []; qInpRight = .5 ;
hConv = 10 ;  Tdiri = 0 ;

problemName = 'HeatRobiAndNeum' ;

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

ndivs     = [ nelem 2 2 ] ;

diriFaces = [  ];
neumFaces = [ 2 qInpRight  ;
              3 0 ;
              4 0 ;
              5 0 ;
              6 0 ] ;
robiFaces = [ 1 hConv ] ;

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

indTimePlot = min( [ size(Ts,2) 20 ] );

figure, hold on, grid on
plot( xs, Ts(:,indTimePlot), 'b-o',  'markersize', MS,'linewidth',LW )
plot( xs, Ts3D(1:(nelem+1), indTimePlot), 'g-s',  'markersize', MS,'linewidth',LW )
plot( xs, Ts3D((1:(nelem+1))+(nelem+1), indTimePlot), 'r-x',  'markersize', MS,'linewidth',LW )
legend('1D','3D-line1','3D-line2')

verifBoolean = ( norm( Ts(indplot,:) - Ts3D(indplot,:) ) / norm( Ts(indplot,:) ) ) < 1e-4 
