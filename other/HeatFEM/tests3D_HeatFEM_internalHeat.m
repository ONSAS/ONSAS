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
Tfinal    = 0.01   ;
rho       = 1.     ;
cSpHe     = 1.     ;
kCond     = 4      ;
Lx        = 1      ;
Ly        = .5     ;
Lz        = .5     ;
Tdiri     = 0      ;

diriDofs = [];
robiDofs = [] ;

nPlots = inf ;

problemName = 'dirichlet3D' ;

initialTempFunc = 'myInitialTemp' ;
internalHeatFunc = 'myInternalHeat' ;

ndivs     = [ 20 10 10 ];

hConv = 10 ;
Tamb = -5 ;
qInpLeft = [];
qInpRight = [];
Tdiri = 0 ;

qInp = 0;

diriFaces = [ Tdiri 1 ];
neumFaces = [ qInp 2 3 4 5 ] ;
robiFaces = [ hConv 6 ] ;
  
[Ts3D, NodesCoord, times ] = HeatFEM( ...
  timeIncr, Tfinal, ...
  [rho, cSpHe, kCond], ...
  [ 2 Lx Ly Lz ], ...
  ndivs, ...
  hConv, diriDofs, robiDofs, Tamb, qInpLeft, qInpRight, Tdiri, ...
  nPlots, problemName, initialTempFunc, internalHeatFunc, ...
  diriFaces, neumFaces, robiFaces  );



% evolution plot
figure, grid on, hold on
MS = 15 ; LW = 2 ;
indplot = round( ndivs(1)/2 ) ;

% numerical solution plot
plot( times(1:5:end), Ts3D(indplot,1:5:end), 'b-s',  'markersize', MS,'linewidth',LW )

